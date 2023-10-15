#include "../include/misc/trie_queue.h"
#include "../include/misc/random.h"
#include "../include/misc/timer.h"

const size_t N = 1025 * 1027;
using Value = uint32_t;
using Index = uint32_t;

int main()
{
  Value* pixels = new Value[N];

  for (Index i = 0; i != N; ++i)
  {
    pixels[i] = i;
  }

  typename pmt::rng<Value>::type rng;

  double avg = 0;
  for (size_t i = 0; i != 50; ++i)
  {

    pmt::random_shuffle(pixels, N, rng);

    {
      pmt::Timer t;

      pmt::TrieQueue<Index> queue(N - 1);

      for (Index i = 0; i != N; ++i)
      {
        queue.insert(pixels[i]);
      }
        
      for (Index i = N; i--;)
      {
        Index val = queue.top();
        queue.remove();
        
        check(val == i);    
      }

      check(queue.empty());

      avg += (t.stop() - avg) / (i + 1.0);
    }
  }


  info("avg time = " << avg);  
  info("Success.");
  
  delete[] pixels;  
}