
#include <cstdint>
#include <iostream>

#include "../include/common.h"
#include "../include/misc/timer.h"
#include "../include/maxtree/connected_components.h"
#include "../include/misc/edge.h"
#include "../include/misc/random.h"

using index_t = uint32_t;
using edge_t = pmt::Edge<index_t>;

int main(int argc, char** argv)
{
  index_t W = 1024;
  index_t H = W;
  index_t D = 1;
  index_t N = W * H * D;

  index_t n_edges = (W - 1U) * H * D + (H - 1U) * W * D + (D - 1U) * W * H;

  edge_t* edges = new edge_t[n_edges];
  edge_t* aux = new edge_t[n_edges];

  index_t* values = new index_t[N];

  typename pmt::rng<index_t>::type rng;

  for (index_t i = 0; i < N; ++i)
  {
    values[i] = rng();
  }

  //values[0] = 0;

  index_t ctr = 0;

  for (index_t d = 0; d < D; ++d)
  for (index_t h = 0; h < H; ++h)
  {
    index_t i = d * H * W + h * W;

    for (index_t w = 0; w < W; ++w)
    {      
      if (w > 0)
      {
        edges[ctr++] = {i - 1, i};
      }

      if (h > 0)
      {
        edges[ctr++] = {i - W, i};
      }

      if (d > 0)
      {
        edges[ctr++] = {i - W*H, i};
      }

      ++i;
    }
  } 

  check(ctr == n_edges);

  index_t* roots = new index_t[N];

  for (index_t i = 0; i < N; ++i)
  {
    roots[i] = i;
  }

  {
    pmt::connected_components(edges, n_edges, aux, values, roots);
  }

  index_t min_root = 0;
  for (index_t i = 0; i < N; ++i)
  {
    if (values[i] < values[min_root])
    {
      min_root = i;
    }
  }
  index_t root = min_root;

  for (index_t i = 0; i < N; ++i)
  {    
    check(roots[i] == root);
    check(roots[i] == roots[roots[i]]);
  }

  out("Root test successful.");

  delete[] roots;
  delete[] values;
  delete[] edges;
  delete[] aux;

  return 0;  
}