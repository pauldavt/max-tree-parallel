#include "opencv2/imgcodecs.hpp"
#include "../include/image/image.h"
#include "../include/maxtree/maxtree.h"
#include "../include/maxtree/tree_scan.h"
#include "../include/maxtree/reconstruct_image.h"

using namespace cv;
using index_t = uint32_t;
using attribute_t = uint32_t;

template <typename value_t>
void process(Mat const& img, attribute_t lambda, char const* out_file, int depth)
{
  check(img.cols > 0);
  check(img.rows > 0);
  check(img.isContinuous());

  info("processing loaded image");
  info("using " << pmt::thread_pool.max_threads() << " threads");

  size_t n = img.cols * img.rows;
  value_t *values = (value_t*)img.ptr<value_t>(0);

  // create an image data structure for the parallel max-tree algorithm
  // 2-D image with 4-neighbors
  using image_t = typename pmt::image<index_t, value_t, 2, 4>::type;
  image_t pmt_img(values, {index_t(img.cols), index_t(img.rows)});

  // allocate memory
  index_t* parents = new index_t[n];
  attribute_t *attributes = new attribute_t[n];

  {
    pmt::Timer t("max-tree construction");

    pmt::maxtree(pmt_img, parents);
  }

  {
    pmt::Timer t("attribute computation");

    // initial area attribute
    auto const &w = [](index_t i) ALWAYS_INL_L(attribute_t)
      {
        return 1U;
      };
    
    // merge (attributes of) pixel set b to pixel set a
    auto const &plus =
      [](attribute_t a, attribute_t b) ALWAYS_INL_L(attribute_t)
      {
        return a + b;
      };

    pmt::tree_scan(parents, n, attributes, w, plus);
  }

  {
    pmt::Timer t("image reconstruction");

    auto const &criterion = [=](index_t i) ALWAYS_INL_L(bool)
    {
      return attributes[i] >= lambda;
    };

    // direct filter: delete nodes in the max-tree that fail the criterion
    // nodes (pixels) take the value of the first ancestor that meets the criterion
    pmt::reconstruct_image(values, n, (value_t*)values, parents, criterion);
  }
  
  // write the filtered image
  imwrite(out_file, img);

  // free memory
  delete[] attributes;
  delete[] parents;
}

int main(int argc, char** argv)
{
  if (argc <= 3)
  {
    std::cout << "Usage: " << argv[0] << " <input image> <output image> <lambda> [thread count]\n";
    return 1;
  }

  if (argc == 5)
  {
    err("changing thread counts needs to be reimplemented");
  }

  attribute_t lambda = std::stoul(argv[3]);

  Mat img = cv::imread(argv[1], IMREAD_GRAYSCALE);

  if (img.data == NULL)
  {
    err("unable to read file");
  }

  switch(img.depth())
  {
    case CV_8U: 
      process<uint8_t>(img, lambda, argv[2], img.depth());
      break;
    case CV_16U:
      process<uint16_t>(img, lambda, argv[2], img.depth());
      break;
    case CV_16S:
      process<int16_t>(img, lambda, argv[2], img.depth());
      break;
    default:
      err("unable to process this image depth");
  }

  return 0;
}