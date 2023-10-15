
#include <cstdint>
#include <iostream>

#include "../include/image/image_blocks.h"
#include "../include/common.h"
#include "../include/misc/timer.h"
#include "../include/maxtree/maxtree_trie.h"
#include "../include/maxtree/maxtree.h"
#include "../include/misc/edge.h"
#include "../include/maxtree/check_equiv.h"

using index_t = uint32_t;

template <typename value_t>
void construct()
{
  index_t W = 300;
  index_t H = W;
  index_t D = W;
  index_t N = W * H * D;

  value_t* vals = new value_t[N];

  using image_t = typename pmt::image<index_t, value_t, 3>::type;  
  image_t img(vals, {W, H, D});
  check(img.dimensions().length() == N);

  using rng = typename pmt::rng<index_t>::type;

  size_t max_threads = pmt::thread_pool.max_threads();
  rng* rand = new rng[max_threads];

  info(max_threads << " threads");

  pmt::thread_pool.for_all(N, [=](index_t i, pmt::thread_nr_t t) {
    if (std::is_floating_point<value_t>::value)
    {
      vals[i] = pmt::random_fp(rand[t]);
    }
    else
    {
      vals[i] = rand[t]();
      //vals[i] = 0;
    }
  });

  index_t* parents = new index_t[N];

  {
    pmt::Timer t;
    pmt::maxtree(img, parents);

    printf("%f megapixel/s\n", N / 1e6 / t.stop());
  }

  index_t* parents2 = new index_t[N];

  {
    pmt::Timer t("seq trie");
    pmt::maxtree_trie(img, parents2);
  }

  {
    out("Checking if the max-tree is correct...");
    pmt::check_equiv(parents, N, parents2, vals);
  }

  delete[] parents2;
  delete[] parents;
  delete[] rand;
  delete[] vals;
}

int main(int argc, char** argv)
{
  construct<float>();

  return 0;
}