//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Alejandro Cabrera 2011.
// Distributed under the Boost
// Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
// copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/bloom_filter for documentation.
//
//////////////////////////////////////////////////////////////////////////////

#ifndef BOOST_BLOOM_FILTER_COUNTING_BLOOM_FILTER_SERIALIZED_HPP
#define BOOST_BLOOM_FILTER_COUNTING_BLOOM_FILTER_SERIALIZED_HPP 1

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/config.hpp>
#include <boost/array.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/size.hpp>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_unsigned.hpp>

#include <boost/bloom_filter/detail/counting_apply_hash.hpp>
#include <boost/bloom_filter/hash/default.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#ifndef BOOST_NO_0X_HDR_INITIALIZER_LIST
#include <initializer_list>
#endif 

using namespace std;
using namespace boost::archive;

namespace boost {
  namespace bloom_filters {
    template <typename T,
	      size_t NumBins,
	      size_t BitsPerBin = 4,
	      class HashFunctions = mpl::vector<boost_hash<T> >,
	      typename Block = size_t>
    class serializeable_counting_bloom_filter {

      // Block needs to be an integral type
      BOOST_STATIC_ASSERT( boost::is_integral<Block>::value == true);

      // Block needs to be an unsigned type
      BOOST_STATIC_ASSERT( boost::is_unsigned<Block>::value == true);

      // BitsPerBin needs to be greater than 0
      BOOST_STATIC_ASSERT( BitsPerBin > 0);

      // it doesn't make sense to ever support using a BitsPerBin value larger
      // than the number of bits per Block. In that case, the user shouldn't
      // be using a Bloom filter to represent their data.
      BOOST_STATIC_ASSERT( (BitsPerBin < (sizeof(Block) * 8) ) );

      // because of the nature of this implementation, the Bloom filter
      // can have internal fragmentation if the calculation for 
      // bins_per_slot has a remainder. The severity of the  internal
      // fragmentation is equal to the remainder * the number of slots.
      // This check prevents internal fragmentation.
      // This also necessarily limits to bin sizes to one of:
      // [1,2,4,8,16,32(64-bit system only)] bits
      BOOST_STATIC_ASSERT( ((sizeof(Block) * 8) % BitsPerBin) == 0);

      // a slot is one element position in the array
      // a bin is a segment of a slot
      static const size_t slot_bits = sizeof(Block) * 8;
      static const size_t bin_bits = NumBins * BitsPerBin;
      static const size_t array_size = bin_bits / slot_bits + 1;

    public:
      typedef T value_type;
      typedef T key_type;
      typedef HashFunctions hash_function_type;
      typedef Block block_type;
      typedef serializeable_counting_bloom_filter<T, NumBins, BitsPerBin, 
				    HashFunctions, Block> this_type;

      typedef boost::array<Block, array_size> bucket_type;
      typedef typename bucket_type::iterator bucket_iterator;
      typedef typename bucket_type::const_iterator bucket_const_iterator;

    private:
      typedef detail::counting_apply_hash<mpl::size<HashFunctions>::value - 1,
					  this_type> apply_hash_type;

    public:
      //* constructors
      serializeable_counting_bloom_filter() 
      {
	this->clear();
      }

      
      template <typename _T, size_t _Bins, size_t _BitsPerBin,
                typename _HashFns, typename _Block>
      static serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin, _HashFns, _Block> load(const string& serialized) {
         size_t numbins, bitsperbin, bits;
         text_iarchive ia{serialized};

         ia.template register_type< serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin, _HashFns, _Block> >();

         serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin, _HashFns, _Block> toret;
         ia >> toret.bits;
         return toret;
      }

      string str() {
         stringstream ss;
         text_oarchive oa{ss};
         oa << this->bits;
         return ss.str();
      }

      template <typename InputIterator>
      serializeable_counting_bloom_filter(const InputIterator start, 
			    const InputIterator end) 
      {
	this->clear();

	for (InputIterator i = start; i != end; ++i)
	  this->insert(*i);
      }

#ifndef BOOST_NO_0X_HDR_INITIALIZER_LIST
      serializeable_counting_bloom_filter(const std::initializer_list<T>& ilist) 
      {
	this->clear();

	typedef typename std::initializer_list<T>::const_iterator citer;
	for (citer i = ilist.begin(), end = ilist.end(); i != end; ++i) {
	  this->insert(*i);
	}
      }
#endif

      //* meta functions
      static BOOST_CONSTEXPR size_t num_bins()
      {
	return NumBins;
      }

      static BOOST_CONSTEXPR size_t bits_per_bin()
      {
	return BitsPerBin;
      }

      static BOOST_CONSTEXPR size_t bins_per_slot()
      {
	return sizeof(block_type) * 8 / BitsPerBin;
      }

      static BOOST_CONSTEXPR size_t mask()
      {
	return static_cast<Block>(0 - 1) >> (slot_bits - BitsPerBin);
      }

      static BOOST_CONSTEXPR size_t bit_capacity() 
      {
        return NumBins * BitsPerBin;
      }

      static BOOST_CONSTEXPR size_t num_hash_functions() 
      {
        return mpl::size<HashFunctions>::value;
      }

      double false_positive_rate() const 
      {
        const double n = static_cast<double>(this->count());
        static const double k = static_cast<double>(num_hash_functions());
        static const double m = static_cast<double>(NumBins);
        static const double e =
	  2.718281828459045235360287471352662497757247093699959574966;
        return std::pow(1 - std::pow(e, -k * n / m), k);
      }

      //? returns the number of bins that have at least 1 bit set
      size_t count() const 
      {
	size_t ret = 0;

	for (bucket_const_iterator i = this->bits.begin(), 
	       end = this->bits.end(); 
	     i != end; ++i) {
	  for (size_t bin = 0; bin < this->bins_per_slot(); ++bin) {
	    const size_t offset_bits = bin * BitsPerBin;
	    const size_t target_bits = (*i >> offset_bits) & this->mask();

	    if (target_bits > 0)
	      ++ret;
	  }
	}

        return ret;
      }

      bool empty() const
      {
	return this->count() == 0;
      }

      const bucket_type&
      data() const
      {
	return this->bits;
      }

      //* core ops
      void insert(const T& t)
      {
	apply_hash_type::insert(t, 
				this->bits,
				this->num_bins());
      }

      template <typename InputIterator>
      void insert(const InputIterator start, const InputIterator end)
      {
	for (InputIterator i = start; i != end; ++i) {
	  this->insert(*i);
	}
      }

      void remove(const T& t)
      {
	apply_hash_type::remove(t, 
				this->bits,
				this->num_bins());
      }

      template <typename InputIterator>
      void remove(const InputIterator start, const InputIterator end)
      {
	for (InputIterator i = start; i != end; ++i) {
	  this->remove(*i);
	}
      }

      bool probably_contains(const T& t) const
      {
	return apply_hash_type::contains(t,
					 this->bits,
					 this->num_bins());
      }

      //* auxiliary ops
      void clear()
      {
	for (bucket_iterator i = bits.begin(), end = bits.end();
	     i != end; ++i) {
	  *i = 0;
	}
      }

      void swap(serializeable_counting_bloom_filter& other)
      {
	serializeable_counting_bloom_filter tmp = other;
	other = *this;
	*this = tmp;
      }

      //* equality comparison operators
      template <typename _T, size_t _Bins, size_t _BitsPerBin,
		typename _HashFns, typename _Block>
      friend bool
      operator==(const serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin,
					     _HashFns, _Block>& lhs,
		 const serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin,
					     _HashFns, _Block>& rhs);

      template <typename _T, size_t _Bins, size_t _BitsPerBin,
		typename _HashFns, typename _Block>
      friend bool
      operator!=(const serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin,
					     _HashFns, _Block>& lhs,
		 const serializeable_counting_bloom_filter<_T, _Bins, _BitsPerBin,
					     _HashFns, _Block>& rhs);


    private:
      bucket_type bits;
    };

    template<class T, size_t NumBins, size_t BitsPerBin, class HashFunctions,
	     typename Block>
    void
    swap(serializeable_counting_bloom_filter<T, NumBins, BitsPerBin, 
			       HashFunctions, Block>& lhs,
	 serializeable_counting_bloom_filter<T, NumBins, BitsPerBin,
			       HashFunctions, Block>& rhs)

    {
      lhs.swap(rhs);
    }

    template<class T, size_t NumBins, size_t BitsPerBin, class HashFunctions,
	     typename Block>
    bool
    operator==(const serializeable_counting_bloom_filter<T, NumBins, BitsPerBin, 
					   HashFunctions, Block>& lhs,
	       const serializeable_counting_bloom_filter<T, NumBins, BitsPerBin,
					   HashFunctions, Block>& rhs)
    {
      return (lhs.bits == rhs.bits);
    }

    template<class T, size_t NumBins, size_t BitsPerBin, class HashFunctions,
	     typename Block>
    bool
    operator!=(const serializeable_counting_bloom_filter<T, NumBins, BitsPerBin, 
					   HashFunctions, Block>& lhs,
	       const serializeable_counting_bloom_filter<T, NumBins, BitsPerBin,
					   HashFunctions, Block>& rhs)
    {
      return !(lhs == rhs);
    }

  } // namespace bloom_filter
} // namespace boost
#endif
