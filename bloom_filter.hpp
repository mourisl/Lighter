/*
 *********************************************************************
 *                                                                   *
 *                           Open Bloom Filter                       *
 *                                                                   *
 * Author: Arash Partow - 2000                                       *
 * URL: http://www.partow.net                                        *
 * URL: http://www.partow.net/programming/hashfunctions/index.html   *
 *                                                                   *
 * Copyright notice:                                                 *
 * Free use of the Open Bloom Filter Library is permitted under the  *
 * guidelines and in accordance with the most current version of the *
 * Common Public License.                                            *
 * http://www.opensource.org/licenses/cpl1.0.php                     *
 *                                                                   *
 *********************************************************************
*/


#ifndef INCLUDE_BLOOM_FILTER_HPP
#define INCLUDE_BLOOM_FILTER_HPP

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
// Add by Li
#include <pthread.h>
#include <memory.h>

#ifdef __AVX__
#include <immintrin.h>
#endif

#ifdef __SSE4_2__
#include <xmmintrin.h>
#endif

#define BLOCK_SIZE 256ull  // Should be multiple of 256's when using AVX 
#define NUM_OF_PATTERN 65536ull 

static const std::size_t bits_per_char = 0x08;    // 8 bits in 1 char(unsigned)
static const unsigned char bit_mask[bits_per_char] = {
                                                       0x01,  //00000001
                                                       0x02,  //00000010
                                                       0x04,  //00000100
                                                       0x08,  //00001000
                                                       0x10,  //00010000
                                                       0x20,  //00100000
                                                       0x40,  //01000000
                                                       0x80   //10000000
                                                     };

class bloom_parameters
{
public:

   bloom_parameters()
   : minimum_size(1),
     maximum_size(std::numeric_limits<unsigned long long int>::max()),
     minimum_number_of_hashes(1),
     maximum_number_of_hashes(std::numeric_limits<unsigned int>::max()),
     projected_element_count(10000),
     false_positive_probability(1.0 / projected_element_count),
     random_seed(0xA5A5A5A55A5A5A5AULL)
   {}

   bloom_parameters( unsigned long long int cnt, double fprate = 0.0001, unsigned long long int seed = 0xA5A5A5A55A5A5A5AULL ):
	   minimum_size(1),
	   maximum_size(std::numeric_limits<unsigned long long int>::max()),
	   minimum_number_of_hashes(1),
	   maximum_number_of_hashes(std::numeric_limits<unsigned int>::max() )
   {
   	projected_element_count = cnt ;
	false_positive_probability = fprate ;
	random_seed = seed ;
	compute_optimal_parameters() ;
   }

   virtual ~bloom_parameters()
   {}

   inline bool operator!()
   {
      return (minimum_size > maximum_size)      ||
             (minimum_number_of_hashes > maximum_number_of_hashes) ||
             (minimum_number_of_hashes < 1)     ||
             (0 == maximum_number_of_hashes)    ||
             (0 == projected_element_count)     ||
             (false_positive_probability < 0.0) ||
             (std::numeric_limits<double>::infinity() == std::abs(false_positive_probability)) ||
             (0 == random_seed)                 ||
             (0xFFFFFFFFFFFFFFFFULL == random_seed);
   }

   //Allowed min/max size of the bloom filter in bits
   unsigned long long int minimum_size;
   unsigned long long int maximum_size;

   //Allowed min/max number of hash functions
   unsigned int minimum_number_of_hashes;
   unsigned int maximum_number_of_hashes;

   //The approximate number of elements to be inserted
   //into the bloom filter, should be within one order
   //of magnitude. The default is 10000.
   unsigned long long int projected_element_count;

   //The approximate false positive probability expected
   //from the bloom filter. The default is the reciprocal
   //of the projected_element_count.
   double false_positive_probability;

   unsigned long long int random_seed;

   struct optimal_parameters_t
   {
      optimal_parameters_t()
      : number_of_hashes(0),
        table_size(0)
      {}

      unsigned int number_of_hashes;
      unsigned long long int table_size;
   };

   optimal_parameters_t optimal_parameters;

   virtual bool compute_optimal_parameters()
   {
      /*
        Note:
        The following will attempt to find the number of hash functions
        and minimum amount of storage bits required to construct a bloom
        filter consistent with the user defined false positive probability
        and estimated element insertion count.
      */

      if (!(*this))
         return false;
      double min_m = std::numeric_limits<double>::infinity();
      double min_k = 0.0;
      double curr_m = 0.0;
      double k = 2.0;

      while (k < 1000.0)
      {
         double numerator   = (- k * projected_element_count) ; 
         double denominator = std::log(1.0 - std::pow(false_positive_probability, 1.0 / k));
         //double denominator = std::log(1.0 - std::pow(10, -4.0 / k));
         curr_m = numerator / denominator + BLOCK_SIZE ;
         if (curr_m < min_m)
         {
            min_m = curr_m;
            min_k = k;
         }
      //printf( "%lf\n", curr_m ) ;
         k += 1.0;
      }
      //exit( 0 ) ;
      optimal_parameters_t& optp = optimal_parameters;

      optp.number_of_hashes = static_cast<unsigned int>(min_k);
      optp.table_size = static_cast<unsigned long long int>(min_m);
      optp.table_size += (((optp.table_size % bits_per_char) != 0) ? (bits_per_char - (optp.table_size % bits_per_char)) : 0);

	//++optp.number_of_hashes ;
      if (optp.number_of_hashes < minimum_number_of_hashes)
         optp.number_of_hashes = minimum_number_of_hashes;
      else if (optp.number_of_hashes > maximum_number_of_hashes)
         optp.number_of_hashes = maximum_number_of_hashes;

      if (optp.table_size < minimum_size)
         optp.table_size = minimum_size;
      else if (optp.table_size > maximum_size)
         optp.table_size = maximum_size;

      return true;
   }

};

class bloom_filter
{
protected:

   typedef std::size_t bloom_type;
   typedef unsigned char cell_type;

public:

   bloom_filter()
   : bit_table_(0),
     salt_count_(0),
     table_size_(0),
     raw_table_size_(0),
     projected_element_count_(0),
     inserted_element_count_(0),
     random_seed_(0),
     desired_false_positive_probability_(0.0),
     numOfThreads(1) 
   {}

   bloom_filter(const bloom_parameters& p)
   : bit_table_(0),
     projected_element_count_(p.projected_element_count),
     inserted_element_count_(0),
     random_seed_((p.random_seed * 0xA5A5A5A5) + 1),
     desired_false_positive_probability_(p.false_positive_probability),
     numOfThreads( 1 ) 
   {
      salt_count_ = p.optimal_parameters.number_of_hashes;
      table_size_ = p.optimal_parameters.table_size;
      generate_unique_salt();
      raw_table_size_ = table_size_ / bits_per_char;
      //printf( "before malloc %d %d\n", (int)raw_table_size_, (int)p.projected_element_count ) ;
      //bit_table_ = new cell_type[static_cast<std::size_t>(raw_table_size_)];
      bit_table_ = new cell_type[(raw_table_size_)];
      std::fill_n(bit_table_,raw_table_size_,0x00);
      /*blockSize = 1 ;
      while ( blockSize <= raw_table_size_ )
      	blockSize *= 2 ;
	blockSize /= 2 ;*/

	// Add by Li

	memset( patterns, 0, sizeof( patterns ) ) ;
	int randomArray[BLOCK_SIZE] ;
	std::size_t i, j ;
	for ( i = 0 ; i < BLOCK_SIZE ; ++i )
		randomArray[i] = i ;

	for ( i = 0 ; i < NUM_OF_PATTERN ; ++i )
	{
		// Pick which bits to set
		for ( j = 0 ; j < salt_count_ ; ++j )
		{
			int k = j + rand() % ( BLOCK_SIZE - j ) ;
			int tmp = randomArray[k] ;
			randomArray[k] = randomArray[j] ;
			randomArray[j] = tmp ;
		}
		// Set the patterns
		for ( j = 0 ; j < salt_count_ ; ++j )
		{
			int r =randomArray[j] ;
			patterns[i][ r / 64 ] |= ( 1ull << (uint64_t)( r % 64 ) ) ;
		}
	}
}

   bloom_filter(const bloom_filter& filter)
   {
      this->operator=(filter);
   }

   inline bool operator == (const bloom_filter& f) const
   {
      if (this != &f)
      {
         return
            (salt_count_                         == f.salt_count_)                         &&
            (table_size_                         == f.table_size_)                         &&
            (raw_table_size_                     == f.raw_table_size_)                     &&
            (projected_element_count_            == f.projected_element_count_)            &&
            (inserted_element_count_             == f.inserted_element_count_)             &&
            (random_seed_                        == f.random_seed_)                        &&
            (desired_false_positive_probability_ == f.desired_false_positive_probability_) &&
            (salt_                               == f.salt_)                               &&
            std::equal(f.bit_table_,f.bit_table_ + raw_table_size_,bit_table_);
      }
      else
         return true;
   }

   inline bool operator != (const bloom_filter& f) const
   {
      return !operator==(f);
   }

   inline bloom_filter& operator = (const bloom_filter& f)
   {
      if (this != &f)
      {
         salt_count_ = f.salt_count_;
         table_size_ = f.table_size_;
         raw_table_size_ = f.raw_table_size_;
         projected_element_count_ = f.projected_element_count_;
         inserted_element_count_ = f.inserted_element_count_;
         random_seed_ = f.random_seed_;
         desired_false_positive_probability_ = f.desired_false_positive_probability_;
         delete[] bit_table_;
         bit_table_ = new cell_type[static_cast<std::size_t>(raw_table_size_)];
         std::copy(f.bit_table_,f.bit_table_ + raw_table_size_,bit_table_);
         salt_ = f.salt_;
      }
      return *this;
   }

   virtual ~bloom_filter()
   {
      delete[] bit_table_;
   }

   inline bool operator!() const
   {
      return (0 == table_size_);
   }

   inline void clear()
   {
      std::fill_n(bit_table_,raw_table_size_,0x00);
      inserted_element_count_ = 0;
   }

   double occupancy()
   {
	   unsigned long long i, j ;
	   unsigned long long used = 0 ;
	   //printf( "(%llu)\n", table_size_ ) ;
	   for ( i = 0, j = 0 ; i <raw_table_size_ ; ++i, j += bits_per_char )
	   {
		int tmp = bit_table_[i] ;
		unsigned long long int t = 0 ;
		//if ( i >= table_size_ - 10 )
		//	printf( "(%llu)%llu %d\n", table_size_,i, tmp ) ;
		while ( tmp != 0 )
		{
			if ( j + t < table_size_ )
				used += ( tmp & 1 ) ;
			tmp = tmp >> 1 ;
			++t ;
		}
	   }
	   //printf( "%llu %llu\n", used, total ) ;
	   return (double)used/table_size_ ;
   }

   double GetActualFP()
   {
   	//return pow( occupancy(), salt_.size() ) ;
	   unsigned long long  i, j, k ;
	   int tagK ;
	   std::size_t used = 0 ;
	   unsigned long long denum = 0 ;
	   double sum = 0 ;
	   double blockFP[BLOCK_SIZE + 1] ;
	   tagK = 0;
	   k = 0;

	   // Precompute the false positive for each possible block's occupancy rate
	   for ( i = 0 ; i <= BLOCK_SIZE ; ++i )
	   	blockFP[i] = pow( (double)i / BLOCK_SIZE, salt_.size() ) ;
	   
	   //printf( "(%llu)\n", table_size_ ) ;
	   for ( i = 0, j = 0 ; i <raw_table_size_ && i < 1000000 ; ++i, j += bits_per_char )
	   {
		int tmp = bit_table_[i] ;
		unsigned long long int t = 0 ;
		//if ( i >= table_size_ - 10 )
		//	printf( "(%llu)%llu %d\n", table_size_,i, tmp ) ;
		while ( t < bits_per_char )
		{
			/*if ( used > BLOCK_SIZE || used < 0 )
			{
				printf( "BIG ERROR\n" ) ;
				printf( "%llu %d %llu\n", i, (int)bit_table_[i], used ) ;
			}*/
			if ( j + t < table_size_ )
				used += ( ( tmp >> ( bits_per_char - t - 1 ) ) & 1 ) ;
			else
				break ;

			if ( j + t == BLOCK_SIZE - 1 )
			{
				tagK = bits_per_char - 1 ;
				k = 0 ;
			}
			else if ( j + t >= BLOCK_SIZE )
			{
				used -= ( bit_table_[k] >> tagK ) & 1 ;
				//exit( 1 ) ;
				--tagK ;		
				if ( tagK < 0 )
				{
					++k ;
					tagK = bits_per_char - 1 ;
				}
			}
			//printf( "%u\n", (unsigned int)used ) ;
			if ( ( j & ( bits_per_char - 1 ) ) == 0 )
			{
				sum += blockFP[used] ; //pow( (double)used / BLOCK_SIZE, salt_.size() ) ;
				++denum ;
			}
			//tmp = tmp >> 1 ;
			++t ;
		}
	   }
	   //printf( "%llu %llu\n", used, total ) ;
	   //return (double)sum/( table_size_ - BLOCK_SIZE + 1 ) ;
	   //printf( "%lf\n", sum / denum ) ;
	   return (double)sum/denum ;
   }

   inline void insert(const unsigned char* key_begin, const std::size_t& length)
   {
	   //std::size_t bit_index = 0;
	   //std::size_t bit = 0;
	   // Implementation of pattern block bloom filter
	   std::size_t start = 0 ;
	   start = ( hash_ap( key_begin, length, salt_[0] ) ) % ( table_size_ - BLOCK_SIZE + 1 ) ;
	   
	   start = start / bits_per_char ; // Align the start position to char
	   int pid = ( hash_ap( key_begin, length, salt_[1] ) ) & ( NUM_OF_PATTERN - 1 ) ; // Which pattern to use
	   int lockId[2] = { -1, -1 } ; // The stride of lock is BLOCK_SIZE, so one block span at most two lock
      //for (std::size_t i = 0 ; i < salt_.size(); ++i)
	   if ( numOfThreads > 1 )
	   {
	   	lockId[0] = ( start / BLOCK_SIZE ) & lockMask ;
		lockId[1] = ( ( ( start + BLOCK_SIZE - 1 ) / BLOCK_SIZE) ) & lockMask ;
		if ( lockId[0] == lockId[1] )
			lockId[1] = -1 ;
		else if ( lockId[0] > lockId[1] )
		{
			int tmp = lockId[1] ;
			lockId[1] = lockId[0] ;
			lockId[0] = tmp ;
		}
		pthread_mutex_lock( &locks[ lockId[0] ] ) ;
		if ( lockId[1] != -1 )
			pthread_mutex_lock( &locks[ lockId[1] ] ) ;
	   }

	   for ( std::size_t j = 0 ; j < BLOCK_SIZE / 64 ; ++j )
	   {
	 	*( (uint64_t *)(  bit_table_ + start + 8 * j ) ) |= patterns[ pid ][ j ] ;
	   }
	   
	   if ( numOfThreads > 1 )
	   {
		if ( lockId[1] != -1 )
			pthread_mutex_unlock( &locks[ lockId[1] ] ) ;
		pthread_mutex_unlock( &locks[ lockId[0]] ) ;
	   }
	   ++inserted_element_count_;
   }

   template<typename T>
   inline void insert(const T& t)
   {
      // Note: T must be a C++ POD type.
      insert(reinterpret_cast<const unsigned char*>(&t),sizeof(T));
   }

   inline void insert(const std::string& key)
   {
      insert(reinterpret_cast<const unsigned char*>(key.c_str()),key.size());
   }

   inline void insert(const char* data, const std::size_t& length)
   {
      insert(reinterpret_cast<const unsigned char*>(data),length);
   }

   template<typename InputIterator>
   inline void insert(const InputIterator begin, const InputIterator end)
   {
      InputIterator itr = begin;
      while (end != itr)
      {
         insert(*(itr++));
      }
   }

   inline virtual bool contains(const unsigned char* key_begin, const std::size_t length) const
   {
      //std::size_t bit_index = 0;
      //std::size_t bit = 0;
      std::size_t start = 0 ;
      start = ( hash_ap( key_begin, length, salt_[0] ) ) % ( table_size_ - BLOCK_SIZE + 1 );
      start /= bits_per_char ;
      int pid = ( hash_ap( key_begin, length, salt_[1] ) ) & ( NUM_OF_PATTERN - 1 ) ;
#ifdef __AVX__ 
      for ( std::size_t j = 0 ; j < BLOCK_SIZE / 256 ; ++j )
      {
        __m256i b = _mm256_loadu_si256((__m256i*)(&bit_table_[start + 32 * j])) ;// bit_table
        __m256i p = _mm256_loadu_si256((__m256i*)(&patterns[pid][4 * j])) ; // pattern
        __m256i a = _mm256_and_si256(b, p) ;
        __m256i eq = _mm256_cmpeq_epi64(a, p) ;
        if (_mm256_movemask_pd((__m256d)eq) != 0x0F)
          return false ;
      }
#else
#ifdef __SSE4_2__
      for ( std::size_t j = 0 ; j < BLOCK_SIZE / 128 ; ++j )
      {
        __m128i b = _mm_loadu_si128((__m128i*)(&bit_table_[start + 16 * j])) ;// bit_table
        __m128i p = _mm_loadu_si128((__m128i*)(&patterns[pid][2 * j])) ; // pattern
        __m128i a = _mm_and_si128(b, p) ;
        __m128i eq = _mm_cmpeq_epi32(a, p) ;
        if (_mm_movemask_epi8(eq) != 0xFFFF)
          return false ;
      }
#else
      for ( std::size_t j = 0 ; j < BLOCK_SIZE / 64 ; ++j )
        if ( ( *(uint64_t *)( &bit_table_[start + 8 * j ] ) & patterns[ pid ][ j ] ) != patterns[ pid ][j] )
          return false ;
#endif
#endif
      return true;
   }

   template<typename T>
   inline bool contains(const T& t) const
   {
      return contains(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
   }

   inline bool contains( const uint64_t &key) const
   {
      return contains(reinterpret_cast<const unsigned char*>(&key),8);
   }

   inline bool contains(const std::string& key) const
   {
      return contains(reinterpret_cast<const unsigned char*>(key.c_str()),key.size());
   }

   inline bool contains(const char* data, const std::size_t& length) const
   {
      return contains(reinterpret_cast<const unsigned char*>(data),length);
   }

   template<typename InputIterator>
   inline InputIterator contains_all(const InputIterator begin, const InputIterator end) const
   {
      InputIterator itr = begin;
      while (end != itr)
      {
         if (!contains(*itr))
         {
            return itr;
         }
         ++itr;
      }
      return end;
   }

   template<typename InputIterator>
   inline InputIterator contains_none(const InputIterator begin, const InputIterator end) const
   {
      InputIterator itr = begin;
      while (end != itr)
      {
         if (contains(*itr))
         {
            return itr;
         }
         ++itr;
      }
      return end;
   }

   inline virtual unsigned long long int size() const
   {
      return table_size_;
   }

   inline std::size_t element_count() const
   {
      return inserted_element_count_;
   }

   inline double effective_fpp() const
   {
      /*
        Note:
        The effective false positive probability is calculated using the
        designated table size and hash function count in conjunction with
        the current number of inserted elements - not the user defined
        predicated/expected number of inserted elements.
      */
      return std::pow(1.0 - std::exp(-1.0 * salt_.size() * inserted_element_count_ / size()), 1.0 * salt_.size());
   }

   inline bloom_filter& operator &= (const bloom_filter& f)
   {
      /* intersection */
      if (
          (salt_count_  == f.salt_count_) &&
          (table_size_  == f.table_size_) &&
          (random_seed_ == f.random_seed_)
         )
      {
         for (std::size_t i = 0; i < raw_table_size_; ++i)
         {
            bit_table_[i] &= f.bit_table_[i];
         }
      }
      return *this;
   }

   inline bloom_filter& operator |= (const bloom_filter& f)
   {
      /* union */
      if (
          (salt_count_  == f.salt_count_) &&
          (table_size_  == f.table_size_) &&
          (random_seed_ == f.random_seed_)
         )
      {
         for (std::size_t i = 0; i < raw_table_size_; ++i)
         {
            bit_table_[i] |= f.bit_table_[i];
         }
      }
      return *this;
   }

   inline bloom_filter& operator ^= (const bloom_filter& f)
   {
      /* difference */
      if (
          (salt_count_  == f.salt_count_) &&
          (table_size_  == f.table_size_) &&
          (random_seed_ == f.random_seed_)
         )
      {
         for (std::size_t i = 0; i < raw_table_size_; ++i)
         {
            bit_table_[i] ^= f.bit_table_[i];
         }
      }
      return *this;
   }

   inline const cell_type* table() const
   {
      return bit_table_;
   }

   inline std::size_t hash_count()
   {
      return salt_.size();
   }

protected:

   inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const
   {
      bit_index = hash % table_size_;
      bit = bit_index % bits_per_char;
   }
   
   inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit, const std::size_t & start ) const
   {
      bit_index = start + ( hash & ( BLOCK_SIZE - 1 ) ) ;
      bit = bit_index % bits_per_char;
   }

   void generate_unique_salt()
   {
      /*
        Note:
        A distinct hash function need not be implementation-wise
        distinct. In the current implementation "seeding" a common
        hash function with different values seems to be adequate.
      */
      const unsigned int predef_salt_count = 128;
      static const bloom_type predef_salt[predef_salt_count] =
                                 {
                                    0xAAAAAAAA, 0x55555555, 0x33333333, 0xCCCCCCCC,
                                    0x66666666, 0x99999999, 0xB5B5B5B5, 0x4B4B4B4B,
                                    0xAA55AA55, 0x55335533, 0x33CC33CC, 0xCC66CC66,
                                    0x66996699, 0x99B599B5, 0xB54BB54B, 0x4BAA4BAA,
                                    0xAA33AA33, 0x55CC55CC, 0x33663366, 0xCC99CC99,
                                    0x66B566B5, 0x994B994B, 0xB5AAB5AA, 0xAAAAAA33,
                                    0x555555CC, 0x33333366, 0xCCCCCC99, 0x666666B5,
                                    0x9999994B, 0xB5B5B5AA, 0xFFFFFFFF, 0xFFFF0000,
                                    0xB823D5EB, 0xC1191CDF, 0xF623AEB3, 0xDB58499F,
                                    0xC8D42E70, 0xB173F616, 0xA91A5967, 0xDA427D63,
                                    0xB1E8A2EA, 0xF6C0D155, 0x4909FEA3, 0xA68CC6A7,
                                    0xC395E782, 0xA26057EB, 0x0CD5DA28, 0x467C5492,
                                    0xF15E6982, 0x61C6FAD3, 0x9615E352, 0x6E9E355A,
                                    0x689B563E, 0x0C9831A8, 0x6753C18B, 0xA622689B,
                                    0x8CA63C47, 0x42CC2884, 0x8E89919B, 0x6EDBD7D3,
                                    0x15B6796C, 0x1D6FDFE4, 0x63FF9092, 0xE7401432,
                                    0xEFFE9412, 0xAEAEDF79, 0x9F245A31, 0x83C136FC,
                                    0xC3DA4A8C, 0xA5112C8C, 0x5271F491, 0x9A948DAB,
                                    0xCEE59A8D, 0xB5F525AB, 0x59D13217, 0x24E7C331,
                                    0x697C2103, 0x84B0A460, 0x86156DA9, 0xAEF2AC68,
                                    0x23243DA5, 0x3F649643, 0x5FA495A8, 0x67710DF8,
                                    0x9A6C499E, 0xDCFB0227, 0x46A43433, 0x1832B07A,
                                    0xC46AFF3C, 0xB9C8FFF0, 0xC9500467, 0x34431BDF,
                                    0xB652432B, 0xE367F12B, 0x427F4C1B, 0x224C006E,
                                    0x2E7E5A89, 0x96F99AA5, 0x0BEB452A, 0x2FD87C39,
                                    0x74B2E1FB, 0x222EFD24, 0xF357F60C, 0x440FCB1E,
                                    0x8BBE030F, 0x6704DC29, 0x1144D12F, 0x948B1355,
                                    0x6D8FD7E9, 0x1C11A014, 0xADD1592F, 0xFB3C712E,
                                    0xFC77642F, 0xF9C4CE8C, 0x31312FB9, 0x08B0DD79,
                                    0x318FA6E7, 0xC040D23D, 0xC0589AA7, 0x0CA5C075,
                                    0xF874B172, 0x0CF914D5, 0x784D3280, 0x4E8CFEBC,
                                    0xC569F575, 0xCDB2A091, 0x2CC016B4, 0x5C5F4421
                                 };

      if (salt_count_ <= predef_salt_count)
      {
         std::copy(predef_salt,
                   predef_salt + salt_count_,
                   std::back_inserter(salt_));
          for (unsigned int i = 0; i < salt_.size(); ++i)
          {
            /*
              Note:
              This is done to integrate the user defined random seed,
              so as to allow for the generation of unique bloom filter
              instances.
            */
            salt_[i] = salt_[i] * salt_[(i + 3) % salt_.size()] + static_cast<bloom_type>(random_seed_);
          }
      }
      else
      {
         std::copy(predef_salt,predef_salt + predef_salt_count,std::back_inserter(salt_));
         srand(static_cast<unsigned int>(random_seed_));
         while (salt_.size() < salt_count_)
         {
            bloom_type current_salt = static_cast<bloom_type>(rand()) * static_cast<bloom_type>(rand());
            if (0 == current_salt) continue;
            if (salt_.end() == std::find(salt_.begin(), salt_.end(), current_salt))
            {
               salt_.push_back(current_salt);
            }
         }
      }
   }

   inline bloom_type hash_ap(const unsigned char* begin, std::size_t remaining_length, bloom_type hash) const
   {
      const unsigned char* itr = begin;
      unsigned int loop = 0;
      while (remaining_length >= 8)
      {
         const unsigned int& i1 = *(reinterpret_cast<const unsigned int*>(itr)); itr += sizeof(unsigned int);
         const unsigned int& i2 = *(reinterpret_cast<const unsigned int*>(itr)); itr += sizeof(unsigned int);
         hash ^= (hash <<  7) ^  i1 * (hash >> 3) ^
              (~((hash << 11) + (i2 ^ (hash >> 5))));
         remaining_length -= 8;
      }
      if (remaining_length)
      {
         if (remaining_length >= 4)
         {
            const unsigned int& i = *(reinterpret_cast<const unsigned int*>(itr));
            if (loop & 0x01)
               hash ^=    (hash <<  7) ^  i * (hash >> 3);
            else
               hash ^= (~((hash << 11) + (i ^ (hash >> 5))));
            ++loop;
            remaining_length -= 4;
            itr += sizeof(unsigned int);
         }
         if (remaining_length >= 2)
         {
            const unsigned short& i = *(reinterpret_cast<const unsigned short*>(itr));
            if (loop & 0x01)
               hash ^=    (hash <<  7) ^  i * (hash >> 3);
            else
               hash ^= (~((hash << 11) + (i ^ (hash >> 5))));
            ++loop;
            remaining_length -= 2;
            itr += sizeof(unsigned short);
         }
         if (remaining_length)
         {
            hash += ((*itr) ^ (hash * 0xA5A5A5A5)) + loop;
         }
      }
      return hash;
   }

   std::vector<bloom_type> salt_;
   unsigned char*          bit_table_;
   unsigned int            salt_count_;
   unsigned long long int  table_size_;
   unsigned long long int  raw_table_size_;
   unsigned long long int  projected_element_count_;
   unsigned int            inserted_element_count_;
   unsigned long long int  random_seed_;
   double                  desired_false_positive_probability_;

   // Add by Li. variables for multi-threads.
   int numOfThreads ;
   pthread_mutex_t *locks ;
   std::size_t lockMask ;
   
   uint64_t patterns[ NUM_OF_PATTERN ][ BLOCK_SIZE / 64 ] ;

public:
	void SetNumOfThreads( int in ) 
	{
		numOfThreads = in ;
		if ( in == 1 )
			return ;
		std::size_t i, k ;
		for ( k = 1 ; k <= (std::size_t)in ; k *= 2 )
			;
		k *= 8 ;
		lockMask = k - 1 ;
		locks = ( pthread_mutex_t * )malloc( sizeof( pthread_mutex_t ) * k ) ;
		for ( i = 0 ; i < k ; ++i )
			pthread_mutex_init( &locks[i], NULL ) ;
	}

	void Output( FILE *fp )
	{
		// Output the patterns
		fwrite( patterns, sizeof( patterns ), 1, fp ) ;

		// Output the bit table
		fwrite( bit_table_, sizeof( unsigned char ), raw_table_size_, fp ) ;
	}

	void Input( FILE *fp )
	{
		fread( patterns, sizeof( patterns ), 1, fp ) ;
		fread( bit_table_,sizeof( unsigned char ), raw_table_size_, fp ) ;
	}
};

inline bloom_filter operator & (const bloom_filter& a, const bloom_filter& b)
{
   bloom_filter result = a;
   result &= b;
   return result;
}

inline bloom_filter operator | (const bloom_filter& a, const bloom_filter& b)
{
   bloom_filter result = a;
   result |= b;
   return result;
}

inline bloom_filter operator ^ (const bloom_filter& a, const bloom_filter& b)
{
   bloom_filter result = a;
   result ^= b;
   return result;
}

class compressible_bloom_filter : public bloom_filter
{
public:

   compressible_bloom_filter(const bloom_parameters& p)
   : bloom_filter(p)
   {
      size_list.push_back(table_size_);
   }

   inline unsigned long long int size() const
   {
      return size_list.back();
   }

   inline bool compress(const double& percentage)
   {
      if ((0.0 >= percentage) || (percentage >= 100.0))
      {
         return false;
      }

      unsigned long long int original_table_size = size_list.back();
      unsigned long long int new_table_size = static_cast<unsigned long long int>((size_list.back() * (1.0 - (percentage / 100.0))));
      new_table_size -= (((new_table_size % bits_per_char) != 0) ? (new_table_size % bits_per_char) : 0);

      if ((bits_per_char > new_table_size) || (new_table_size >= original_table_size))
      {
         return false;
      }

      desired_false_positive_probability_ = effective_fpp();
      cell_type* tmp = new cell_type[static_cast<std::size_t>(new_table_size / bits_per_char)];
      std::copy(bit_table_, bit_table_ + (new_table_size / bits_per_char), tmp);
      cell_type* itr = bit_table_ + (new_table_size / bits_per_char);
      cell_type* end = bit_table_ + (original_table_size / bits_per_char);
      cell_type* itr_tmp = tmp;

      while (end != itr)
      {
         *(itr_tmp++) |= (*itr++);
      }

      delete[] bit_table_;
      bit_table_ = tmp;
      size_list.push_back(new_table_size);

      return true;
   }

private:

   inline void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const
   {
      bit_index = hash;
      for (std::size_t i = 0; i < size_list.size(); ++i)
      {
         bit_index %= size_list[i];
      }
      bit = bit_index % bits_per_char;
   }

   std::vector<unsigned long long int> size_list;
};

#endif


/*
  Note 1:
  If it can be guaranteed that bits_per_char will be of the form 2^n then
  the following optimization can be used:

  hash_table[bit_index >> n] |= bit_mask[bit_index & (bits_per_char - 1)];

  Note 2:
  For performance reasons where possible when allocating memory it should
  be aligned (aligned_alloc) according to the architecture being used.
*/
