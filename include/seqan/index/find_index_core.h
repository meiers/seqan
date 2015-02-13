// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_FIND_INDEX_CORE_H
#define SEQAN_HEADER_INDEX_FIND_INDEX_CORE_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// different layouts of a suffix array or lcp table

	struct SortedList {};			// classical sorted list (suffix array, sorted list, ...)
	struct LeftCompleteTree {};		// flattened search tree root, left child, right child, left child's left child, left child's right child, ...

	template < unsigned BlockSize = 4096 >
	struct BTree {};				// b-tree compacts nodes and its children to blocks of BlockSize

	template < typename TString, typename TSpec >
    class SearchTreeIterator {};


	//////////////////////////////////////////////////////////////////////////////
	// class to return a suffix given a suffix start position
    //

    template <typename TText, typename TSAValue>
    struct SuffixFunctor :
        std::unary_function<TSAValue, typename Suffix<TText const>::Type>
    {
        TText const &text;

        SuffixFunctor(TText const &text) :
        text(text)
        {}

        typename SuffixFunctor::result_type
        operator() (TSAValue const &pos) const
        {
            return suffix(text, pos);
        }
    };

    template <typename TText, typename TSpec, typename TSAValue>
    struct SuffixFunctor <Index<TText, TSpec>, TSAValue> :
        public SuffixFunctor<typename Fibre<Index<TText, TSpec>, FibreText>::Type, TSAValue>
    {
        typedef SuffixFunctor<typename Fibre<Index<TText, TSpec>, FibreText>::Type, TSAValue> TBaseClass;

        SuffixFunctor(Index<TText, TSpec> const &index) : TBaseClass(indexText(index))
        {}
    };


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <typename TString>
	class SearchTreeIterator< TString, SortedList >
	{
	public:

		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		inline SearchTreeIterator(TString &string):
			first(begin(string, Standard())),
			count(length(string))
		{
			count2 = count / 2;
			_mid = first;
			goFurther(_mid, count2);
		}

        inline SearchTreeIterator(TIterator &first, TSize count):
            first(first),
            count(count)
		{
			count2 = count / 2;
			_mid = first;
			goFurther(_mid, count2);
		}

        inline const TValue& operator*() const {
			return *_mid;
		}

        inline const TValue* operator->() const {
			return &*_mid;
		}

		inline TSize mid() {
			return count2;
		}

		// descend left
		inline SearchTreeIterator & left()
		{
			count = count2;
			count2 /= 2;
			_mid = first;
			goFurther(_mid, count2);
			return *this;
		}

		// descend right
        inline SearchTreeIterator & right()
		{
			first = ++_mid, count -= count2 + 1;
			count2 = count / 2;
			goFurther(_mid, count2);
			return *this;
		}

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline bool atEnd() {
            return !count;
        }

		inline operator TIterator & () {
			return _mid;
		}

	private:
		TIterator	first, _mid;
		TSize		count, count2;
	};


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//
    // deprecated: Was used for the practically-slower enhanced lcp table search

	template <typename TString>
	class SearchTreeIterator< TString, LeftCompleteTree > 
	{
	public:
		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		TIterator	it;
		TSize		size;

		inline SearchTreeIterator(TString &string, TSize _size):
			it(begin(string, Standard())),
			size(_size)
		{
			_left = 0;
			_lSize = 1;
            _iSize = size;
			for (_xSize = 1; _xSize < size; _xSize <<= 1)
                continue;  // Create smallest 2^k that is >= size.
            if (!size) _xSize = 0;
		}

		inline SearchTreeIterator():
			it(),
            size(0),
            _xSize(0) {}

        inline const TValue& operator*() const {
			return *it;
		}

        inline const TValue* operator->() const {
			return &*it;
		}

		inline TSize mid() {
			return _xSize >> 1;
		}

		inline TSize leftSize() {
			return mid();
		}

		inline TSize rightSize() {
			return _iSize - mid();
		}

		// descend left
        inline SearchTreeIterator & left()
		{
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
            _descendLeft();
			_iSize = _xSize;	    // = mid();
            return *this;
        }

		// descend right
        inline SearchTreeIterator & right()
		{
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
			_iSize -= mid();
            SEQAN_ASSERT_NEQ(_iSize, 0u);    // _xSize/2 is less than _iSize by invariant

            // step down right
            _descendRight();

            // step down left until we have two sons or are a leaf
            while (_iSize <= (_xSize >> 1))
                _descendLeft();

            return *this;
        }

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline SearchTreeIterator& operator--() {
            --it;
            return *this;
        }

        inline SearchTreeIterator operator--(int) {
            SearchTreeIterator before = *this;
            --it;
            return before;
        }

        inline SearchTreeIterator & operator++() {
            ++it;
            return *this;
        }

        inline SearchTreeIterator operator++(int) {
            SearchTreeIterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(SearchTreeIterator const &I) {
            return (_xSize == I._xSize) && (_xSize == 0 || it == I.it);
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool atEnd() {
            return !_xSize;
        }

		inline operator TIterator & () {
			return it;
		}

	private:
		TSize _left;		// left iterator offset of current interval
		TSize _lSize;		// iterator elements of current level
		TSize _xSize;		// max interval size of current level
		TSize _iSize;		// current interval size

        inline void _descendLeft() 
		{
			goFurther(it, _left + _lSize);
			_left <<= 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }

        inline void _descendRight() 
		{
			goFurther(it, _left + 1 + _lSize);
			_left = (_left << 1) + 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search b-tree like a real tree
	//

	template <typename TString, unsigned BlockSize>
	class SearchTreeIterator< TString, BTree<BlockSize> > 
	{
	public:
		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		enum { BlockHeight = Log2Floor<BlockSize>::VALUE };
		enum { BlockElements = (1 << BlockHeight) - 1 };
		enum { BlockInnerElements = (1 << (BlockHeight - 1)) - 1 };
		enum { BlockLeafs = 1 << (BlockHeight - 1) };

		TIterator	it;
		TSize		size;

		inline SearchTreeIterator(TString &string):
			it(begin(string, Standard())),
			size(length(string))
		{
			//_left = 0;
			//_lSize = 0;
   //         _iSize = size;

			_heightHigh = 1;
			_stepSizeLow = 1;
			for(TSize _xSizeLow = 2; _xSizeLow <= size; _xSizeLow <<= 1) {
				if (_stepSizeLow == BlockLeafs) {
					_stepSizeLow = 1;
					++_heightHigh;
				} else
					_stepSizeLow <<= 1;
			}
			
			_stepSizeLow >>= 1;
			for(_xSizeHigh = 1; _xSizeHigh * BlockSize <= size; _xSizeHigh *= BlockSize) ;

			_leftLow = (_stepSizeLow << 1) - 1;		// point to the middle
			_leftHigh = 0;
			_lSizeHigh = 1;

			it += _leftLow;

			_updateElements();

			//if (!size) _xSizeLow = 0;
   //         if (!size) _xSizeHigh = 0;
		}

		inline SearchTreeIterator():
			it() {}

        inline const TValue& operator*() const {
			return *it;
		}

		inline const TValue* operator->() const {
			return &*it;
		}

		// descend left
        inline SearchTreeIterator & left() 
		{
            if (_heightHigh == 1 && !_stepSizeLow) {
                _heightHigh = 0;
                return *this;
            }
            _descendLeft();
            return *this;
        }

		// descend right
        inline SearchTreeIterator & right() 
		{
            if (_heightHigh == 1 && !_stepSizeLow) {
				++it;
                _heightHigh = 0;
                return *this;
            }

            // step down right
            _descendRight();

            // step down left until we have two sons or are a leaf
            while (_elements <= _leftLow)
                _descendLeft();

            return *this;
        }

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline SearchTreeIterator& operator--() {
            --it;
            return *this;
        }

        inline SearchTreeIterator operator--(int) {
            SearchTreeIterator before = *this;
            --it;
            return before;
        }

        inline SearchTreeIterator& operator++() {
            ++it;
            return *this;
        }

        inline SearchTreeIterator operator++(int) {
            SearchTreeIterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(SearchTreeIterator const &I) {
			return (it == I.it) || (atEnd() && I.atEnd());
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool atEnd() const {
            return !_heightHigh;
        }

		inline operator TIterator & () {
			return it;
		}

	private:
		TSize _heightHigh;	// height measured in BBlocks
		unsigned _elements;		// elements in current BBlock

		unsigned _leftLow;		// left iterator offset of current interval
		unsigned _stepSizeLow;	// left/right step size of current level

		TSize _leftHigh;		// left BBlock offset of current interval
		TSize _lSizeHigh;	// BBlocks of current level
		TSize _xSizeHigh;	// max BBlocks of current level

		inline void _descendLeft() 
		{
			if (_stepSizeLow) {
				it -= _stepSizeLow;
				_leftLow -= _stepSizeLow;
				_stepSizeLow >>= 1;
			} else
				if (--_heightHigh) {
					_descendBBlock(_leftLow);

					_leftLow = BlockInnerElements;		// point to the middle
					_stepSizeLow = BlockLeafs / 2;		// with step size of half an interval length
				}
        }

		inline void _descendRight() 
		{
			if (_stepSizeLow) {
				it += _stepSizeLow;
				_leftLow += _stepSizeLow;
				_stepSizeLow >>= 1;
			} else
				if (--_heightHigh) {
					_descendBBlock(_leftLow + 1);

					_leftLow = BlockInnerElements;		// point to the middle
					_stepSizeLow = BlockLeafs / 2;		// with step size of half an interval length
				}
        }

		inline void _updateElements() 
		{
			TSize firstElement = (1 + _leftHigh * BlockSize) * _xSizeHigh - 1;
			TSize lastElement = (1 + (_leftHigh + 1) * BlockSize) * _xSizeHigh - 2;

			if (lastElement >= size)
				_elements = (size - firstElement) / _xSizeHigh;
			else
				_elements = BlockElements;
		}

        inline void _descendBBlock(TSize _childIndex) 
		{
			// goFurther to the begin of the current BBlock and further to the destination BBlock
			goFurther(it, BlockSize * (_leftHigh * (BlockSize - 1) + _childIndex + _lSizeHigh) + BlockInnerElements - _leftLow);

			_leftHigh = _leftHigh * BlockSize + _childIndex;
			_xSizeHigh /= BlockSize;
			_lSizeHigh = (size / _xSizeHigh + BlockSize - 1) / BlockSize;

			_updateElements();
        }
    };

	template <
		typename TStringDereferer,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_lowerBoundSA(
        TStringDereferer &dereferer,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query,
        typename Size<TQuery>::Type parentRepLen = 0)
	{	// find first element not before query, using operator<
		typedef typename Size<TQuery>::Type                 TSize;
		typedef typename TStringDereferer::result_type      TString;
		typedef typename Iterator<TString, Standard>::Type	TStringIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;
        typedef typename Value<TString>::Type               TStringValue;

		TSize lcpLower = 0;
		TSize lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.atEnd(); )
		{	// divide and conquer, find half that contains answer

			TString		suf = dereferer(*treeIter);
			TStringIter	t = begin(suf, Standard());
			TStringIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TSize		lcp = _min(lcpLower, lcpUpper);

			goFurther(t, lcp + parentRepLen);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && ordEqual(*t, convert<TStringValue>(*q))) {
				++t;
				++q;
				++lcp;
			}
			
            // is text < query ?
			if (q != qEnd && (t == tEnd || ordLess(*t, convert<TStringValue>(*q)))) {
				treeIter.right();
				lcpLower = lcp;
			} else {
				treeIter.left();
				lcpUpper = lcp;
			}
		}
		return treeIter;
	}

	template <
		typename TStringDereferer,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_upperBoundSA(
        TStringDereferer &dereferer,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query,
        typename Size<TQuery>::Type parentRepLen = 0)
	{	// find first element that query is before, using operator<
		typedef typename Size<TQuery>::Type                 TSize;
		typedef typename TStringDereferer::result_type      TString;
		typedef typename Iterator<TString, Standard>::Type	TStringIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;
        typedef typename Value<TString>::Type               TStringValue;

		TSize lcpLower = 0;
		TSize lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.atEnd(); )
		{	// divide and conquer, find half that contains answer

			TString		suf = dereferer(*treeIter);
			TStringIter	t = begin(suf, Standard());
			TStringIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TSize		lcp = _min(lcpLower, lcpUpper);

			goFurther(t, lcp + parentRepLen);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && ordEqual(*t, convert<TStringValue>(*q))) {
				++t;
				++q;
				++lcp;
			}
			
            // is text <= query ?
			if (q == qEnd || t == tEnd || !ordGreater(*t, convert<TStringValue>(*q))) {
				treeIter.right();
				lcpLower = lcp;
			} else {
				treeIter.left();
				lcpUpper = lcp;
			}
		}
		return treeIter;
	}

	//////////////////////////////////////////////////////////////////////////////
    // binary search with mlr-heuristic
	template <
		typename TStringDereferer,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	_equalRangeSA(
		TStringDereferer &dereferer,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query,
        typename Size<TQuery>::Type parentRepLen = 0)
	{	// find range equivalent to query, using operator<
		typedef typename Size<TQuery>::Type                 TSize;
		typedef typename TStringDereferer::result_type      TString;
		typedef typename Iterator<TString, Standard>::Type	TStringIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;
        typedef typename Value<TString>::Type               TStringValue;

		TSize lcpLower = 0;
		TSize lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.atEnd(); )
		{	// divide and conquer, check midpoint

			TString		suf = dereferer(*treeIter);
			TStringIter	t = begin(suf, Standard());
			TStringIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TSize       lcp = _min(lcpLower, lcpUpper);

			goFurther(t, lcp + parentRepLen);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && ordEqual(*t, convert<TStringValue>(*q))) {
				++t;
				++q;
				++lcp;
			}
			
            // is text < query ?
			if (q != qEnd && (t == tEnd || ordLess(*t, convert<TStringValue>(*q))))
			{	// range begins above mid, loop
				treeIter.right();
				lcpLower = lcp;
			}
            // is text > query ?
			else if (q != qEnd && (t != tEnd && ordGreater(*t, convert<TStringValue>(*q))))
			{	// range in first half, loop
				treeIter.left();
				lcpUpper = lcp;
			} else
            // is text == query ?
			{	// range straddles mid, find each end and return
				return Pair<TSAIter> (
					_lowerBoundSA(dereferer, treeIter.leftChild(), query, parentRepLen),
					_upperBoundSA(dereferer, treeIter.rightChild(), query, parentRepLen)
				);
			}
		}
		return Pair<TSAIter> (treeIter, treeIter);	// empty range
	}


	//////////////////////////////////////////////////////////////////////////////
	// Finder wrappers (return iterators instead of positions)

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSANaiveIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSANaive(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	upperBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA const, Standard>::Type >
	equalRangeSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _equalRangeSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	upperBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA const, Standard>::Type >
	equalRangeSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _equalRangeSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}


	//////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find first element not before query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element not before query, using operator<
		return lowerBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find first element that query is before, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return upperBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA const, Standard>::Type > itPair =
			_equalRangeSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return equalRangeSA(text, sa, query, SortedList());
	}


	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find first element not before query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element not before query, using operator<
		return lowerBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find first element that query is before, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return upperBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA const, Standard>::Type > itPair =
			_equalRangeSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return equalRangeSA(text, sa, query, SortedList());
	}

}
#endif // #ifdef SEQAN_HEADER_INDEX_FIND_INDEX_CORE_H

