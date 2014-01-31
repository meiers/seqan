// ==========================================================================
//                       index_gapped_dislex_external.h
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_DISLEX_EXTERNAL_H_
#define CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_DISLEX_EXTERNAL_H_

namespace SEQAN_NAMESPACE_MAIN
{
    template <typename TShape, typename TSACA = Skew7>
    struct DislexExternal {};
// --------------------------------------------------------------------------
// Filter to transform positions to lengths                          [String]
// --------------------------------------------------------------------------

template <typename TSize>
struct _positionToLengthTransform
{
    TSize N;
    _positionToLengthTransform(TSize len): N(len)
    {}

    TSize operator()(TSize pos) const
    {
        return N - pos;
    }
};

// --------------------------------------------------------------------------
// Filter to transform positions to lengths                       [StringSet]
// --------------------------------------------------------------------------

template <typename TLimitString, typename TIdAndPosPair>
struct _positionToLengthTransformMulti
{
    typedef typename Value<TIdAndPosPair, 2>::Type TSize;
    TLimitString const & limits;

    _positionToLengthTransformMulti(TLimitString const & limitStr) : limits(limitStr)
    {}

    inline TIdAndPosPair operator() (TIdAndPosPair x) const
    {
        TSize N = limits[x.i1+1] - limits[x.i1];
        x.i2 = N - x.i2;
        return x;
    }
};


// --------------------------------------------------------------------------
// Comparator for naming tuples                                      [String]
// --------------------------------------------------------------------------

// TODO: _dislexTupleComp for bitvectors

/*
 * @signature _dislexTupleComp<TValue, TShape, TResult = int>
 *
 * @tparam TValue expects a Pair<TSize, Tuple> where the 1st parameter contains the
 *                <b>length</b> of the underlying suffix and the 2nd parameter the
 *                fixed-length sequence tuple (possibly bitpacked)
 * @tparam TShape expects a fixed CyclicShape (CyclicShape<FixedShape<...> >)
 *
 * Only for hardwired cyclic shapes!!
 *
 * @see _dislexTupleCompMulti
 */
template <typename TValue, typename TShape, typename TResult=int>
struct _dislexTupleComp : public std::binary_function<TValue, TValue, TResult>
{
    typedef typename Value<TValue, 1>::Type                 TSize;
    typedef typename Value<TValue, 2>::Type                 TTuple;
    typedef typename Value<TTuple>::Type                    TTupleValue;
    typedef typename StoredTupleValue_<TTupleValue>::Type TStoredValue;

    enum
    {
        _span = TShape::span,
        _weight = WEIGHT<TShape>::VALUE
    };
    TSize realLengths[2*_span];
    _positionToLengthTransform<TSize> posToLen;

    _dislexTupleComp(TShape const & shape, TSize strLen) : posToLen(strLen)
    {
        cyclicShapeToSuffixLengths(realLengths, shape);

        // extend the table to a size of 2*_span
        for (unsigned i=0; i<_span; ++i)
            realLengths[i+_span] = _weight + realLengths[i];
    }

    inline TResult operator() (const TValue &a, const TValue &b) const
    {
        const TStoredValue * sa = a.i2.i;
        const TStoredValue * sb = b.i2.i;

        TSize la = posToLen(a.i1);
        TSize lb = posToLen(b.i1);

        /*
         // the usual case:
         // both tuples have more than _weight real characters.
         // TODO: Tell compiler this is the more likely if-branch?
         if(la >= 2* _span && lb >= 2*_span)
         {
         // TODO: Unroll loop using Loop struct?
         for (TSize i = 0; i < _weight; ++i, ++sa, ++sb)
         {
         if (*sa == *sb) continue;
         return (*sa < *sb)? -1 : 1;
         }
         return 0;
         }
         */

        // find out the real lengths of the gapped strings
        TSize rla = (la < 2*_span ? realLengths[la] : 2*_weight);
        TSize rlb = (lb < 2*_span ? realLengths[lb] : 2*_weight);

        // compare the overlap of the first n bases
        TSize n = std::min(static_cast<TSize>(_weight), std::min(rla, rlb) );
        for (TSize i = 0; i < n; i++, ++sa, ++sb)
        {
            if (*sa == *sb) continue;
            return (*sa < *sb)? -1 : 1;
        }

        // if both strings have more than _weight chars, they are equal.
        if (rla > _weight && rlb > _weight) return 0;

        // if they differ in size, the shorter one is smaller.
        if (rla != rlb)
            return (rla < rlb ? -1 : 1);

        // In the same string, the length of the underlying suffix decides:
        if (la < lb) return -1;
        if (la > lb) return 1;

        // last case a.i1 == b.i1
        return 0;
    }
};

// --------------------------------------------------------------------------
// Comparator for naming tuples                                   [StringSet]
// --------------------------------------------------------------------------

// TODO: _dislexTupleCompMulti for bitvectors

/*
 * @signature _dislexTupleCompMulti<TValue, TShape, TResult = int>
 *
 * @tparam TValue expects a Pair<Pair<TSize, TSize>, Tuple> where the 1st parameter
 *                is a Pair of sequence ID and suffix <b>length</b> and the 2nd parameter
 *                the fixed-length sequence tuple (possibly bitpacked)
 * @tparam TShape expects a fixed CyclicShape (CyclicShape<FixedShape<...> >)
 *
 * @see _dislexTupleComp
 */
template <typename TValue, typename TShape, typename TLimitString, typename TResult=int>
struct _dislexTupleCompMulti  : public std::binary_function<TValue, TValue, TResult>
{
    typedef typename Value<TValue, 1>::Type                 TSetPos;
    typedef typename Value<TSetPos, 2>::Type                TSize;
    typedef typename Value<TValue, 2>::Type                 TTuple;
    typedef typename Value<TTuple>::Type                    TTupleValue;
    typedef typename StoredTupleValue_<TTupleValue>::Type TStoredValue;

    enum {
        _span = TShape::span,
        _weight = WEIGHT<TShape>::VALUE
    };
    TSize realLengths[2*_span];
    _positionToLengthTransformMulti<TLimitString, TSetPos> posToLen;

    _dislexTupleCompMulti(TShape const & shape, TLimitString const & limits) : posToLen(limits)
    {
        cyclicShapeToSuffixLengths(realLengths, shape);

        // extend the table to a size of 2*_span
        for (unsigned i=0; i<_span; ++i)
            realLengths[i+_span] = _weight + realLengths[i];
    }


    inline TResult operator() (const TValue &a, const TValue &b) const
    {
        const TStoredValue * sa = a.i2.i;
        const TStoredValue * sb = b.i2.i;

        TSetPos la = posToLen(a.i1);
        TSetPos lb = posToLen(b.i1);

        /*
         // the usual case:
         // both tuples have more than _weight real characters.

         // TODO: Tell compiler this is the more likely if-branch?
         if(la.i2 >= 2*_span && lb.i2 >= 2*_span)
         {
         // TODO: Unroll loop using Loop struct?
         for (TSize i = 0; i < _weight; ++i, ++sa, ++sb)
         {
         if (*sa == *sb) continue;
         return (*sa < *sb)? -1 : 1;
         }
         return 0;
         }
         */

        // find out the real lengths of the gapped strings
        TSize rla = (la.i2 < 2*_span ? realLengths[la.i2] : 2*_weight);
        TSize rlb = (lb.i2 < 2*_span ? realLengths[lb.i2] : 2*_weight);

        // compare the overlap of the first n bases
        TSize n = std::min(static_cast<TSize>(_weight), std::min(rla, rlb) );
        for (TSize i = 0; i < n; i++, ++sa, ++sb)
        {
            if (*sa == *sb) continue;
            return (*sa < *sb)? -1 : 1;
        }

        // if both strings have more than _weight chars, they are equal.
        if (rla > _weight && rlb > _weight) return 0;

        // if they differ in size, the shorter one is smaller.
        if (rla != rlb)
            return (rla < rlb ? -1 : 1);

        // if both have the same number of chars,
        // at first the string ID is relevant:
        if (la.i1 > lb.i1) return -1;
        if (la.i1 < lb.i1) return 1;

        // In the same string, the length of the underlying suffix decides:
        if (la.i2 < lb.i2) return -1;
        if (lb.i2 > lb.i2) return 1;

        // last case a.i1 == b.i1
        return 0;
    }
};


// --------------------------------------------------------------------------
// Mapping functor from text to lexText                              [String]
// --------------------------------------------------------------------------

// wrapper for the dislex Pipe
// takes a tuple <l, ACGACA> where p is the suffix position
// and returns L(N-l)
template <
typename TValue,
typename TResult = typename Value<TValue, 1>::Type>
struct _dislexMap :
public std::unary_function<TValue, TResult>
{
    _dislexTransform<TResult> formula;

    _dislexMap(TResult S_, TResult N_) : formula(S_, N_)
    {}

    inline TResult operator() (const TValue & x) const
    {
        return formula(x.i1);
    }
};


// --------------------------------------------------------------------------
// Mapping functor from text to lexText                           [StringSet]
// --------------------------------------------------------------------------

// dislex transformation used in the mapper pool
// takes a Pair( Pair(s,p), ACGATCG), where s is the seq id and p the suffix position,
// returns a global position L(s,p)=pos
template <
typename TValue,
typename Tlimits,
typename TResult = typename Value<typename Value<TValue, 1>::Type, 2>::Type >
struct _dislexMapMulti :
public std::unary_function<TValue, TResult>
{
    typedef typename Value<TValue, 1>::Type TPair;
    typedef typename Value<TPair, 2>::Type TSize;

    _dislexTransformMulti<TPair, Tlimits> formula;
    
    _dislexMapMulti(TResult S_, Tlimits const & stringSetLimits) : formula(S_, stringSetLimits)
    {}
    
    inline TResult operator() (const TValue & x) const
    {
        return formula(x.i1);
    }
};

// --------------------------------------------------------------------------
// Pipe Dislex                                                       [String]
// --------------------------------------------------------------------------

    // TODO(meiers): Define metafunctions, e.g. Value, or do I need them?

template <typename TInput, typename TShape, typename TSACA>
struct Pipe<TInput, DislexExternal<TShape, TSACA> >
{
    typedef Pipe<TInput, GappedTupler<TShape, false> >          TPipeTupler;
    typedef _dislexTupleComp<TypeOf_(TPipeTupler), TShape>      TTupleComparator;
    typedef Pool<TypeOf_(TPipeTupler), SorterSpec<
            SorterConfigSize<TTupleComparator,
            TSizeOf_(TPipeTupler) > > >                         TPoolSorter;

    typedef Pipe< TPoolSorter, Namer<TTupleComparator> >        TPipeNamer;
    typedef _dislexMap<TypeOf_(TPipeNamer) >                    TDislexMapper;
    typedef Pool< TypeOf_(TPipeNamer), MapperSpec<
            MapperConfigSize< TDislexMapper,
            TSizeOf_(TPipeNamer) > > >                          TPoolMapper;

    typedef Pipe< TPoolMapper, Filter<
            filterI2<TypeOf_(TPoolMapper)> > >                  TPipeFilterI2;
    typedef Pipe<TPipeFilterI2, TSACA>                          TPipeSACA;
    typedef _dislexReverseTransform<TypeOf_(TPipeSACA),
            TypeOf_(Pipe)>                                      TDislexReverse;
    typedef Pipe<TPipeSACA, Filter<TDislexReverse> >            TPipeReverseTransform;


    TPipeSACA pool;                 // last pool (skew); will be filled when calling process().
    TPipeReverseTransform in;       // final Pipe

    Pipe()
    {}

    Pipe(TInput & textIn) : in(pool, TDislexReverse(TShape::span, length(textIn)))
    {
        // fill pool right away
        process(textIn);
    }

    inline typename Value<Pipe>::Type const operator*() {
        return *in;
    }

    inline Pipe& operator++() {
        ++in;
        return *this;
    }

    template < typename TInput_ >
    bool process(TInput_ &textIn)
    {
        // 1. Generate Gapped Tuples
        TPipeTupler                                             tupler(textIn);

        // 2. Sort Tuples by the first few characters
        TTupleComparator                                        _comp(TShape(), length(textIn));
        TPoolSorter                                             sorter(tupler, _comp);
        sorter << tupler;

        // 3. Name tuples by their rank
        TPipeNamer                                              namer(sorter, _comp);

        // 4. Map text Positions to lexText positions
        TDislexMapper                                           _map(TShape::span, length(textIn));
        TPoolMapper                                             mapper(namer, _map);
        mapper << namer;

        // 5. Discard positions, keep rank
        TPipeFilterI2                                           filter(mapper);

        // 6. Run SACA on lex text
        pool << filter;

        // 7. Reverse Transform is done during the reading process
        return true;
    }
};

// --------------------------------------------------------------------------
// Pipe Dislex                                                    [StringSet]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, typename TSACA, typename TPair, typename TLimits>
struct Pipe<TInput, Multi<DislexExternal<TShape, TSACA>, TPair, TLimits> >
{
    typedef Pipe<TInput, Multi<GappedTupler<TShape, false>,
            TPair, TLimits> >                                   TPipeTupler;
    typedef _dislexTupleCompMulti<TypeOf_(TPipeTupler),
            TShape, TLimits>                                    TTupleComparator;
    typedef Pool<TypeOf_(TPipeTupler), SorterSpec<
            SorterConfigSize<TTupleComparator,
            TSizeOf_(TPipeTupler) > > >                         TPoolSorter;

    typedef Pipe< TPoolSorter, Namer<TTupleComparator> >        TPipeNamer;
    typedef _dislexMapMulti<TypeOf_(TPipeNamer), TLimits>       TDislexMapper;
    typedef Pool< TypeOf_(TPipeNamer), MapperSpec<
            MapperConfigSize< TDislexMapper,
            TSizeOf_(TPipeNamer) > > >                          TPoolMapper;

    typedef Pipe< TPoolMapper, Filter<
            filterI2<TypeOf_(TPoolMapper)> > >                  TPipeFilterI2;
    typedef Pipe<TPipeFilterI2, TSACA>                          TPipeSACA;
    typedef _dislexReverseTransformMulti<TypeOf_(TPipeSACA),
            TLimits, TypeOf_(Pipe)>                             TDislexReverse;
    typedef Pipe<TPipeSACA, Filter<TDislexReverse> >            TPipeReverseTransform;


    TLimits const & limits;         // StringSetLimits
    TPipeSACA pool;                 // last pool (skew); will be filled when calling process().
    TPipeReverseTransform in;       // final Pipe


    template <typename TLimits_>
    Pipe(TLimits_ const & _limits, SEQAN_CTOR_ENABLE_IF(IsSameType<TLimits, TLimits_>)) :
        limits(_limits)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TLimits_>
    Pipe(TInput& _textIn, TLimits_ const & _limits, SEQAN_CTOR_ENABLE_IF(IsSameType<TLimits, TLimits_>)) :
        limits(_limits), in(pool, TDislexReverse(TShape::span, limits))
    {
        // fill pool right away
        process(_textIn);
    }

    inline typename Value<Pipe>::Type const operator*() {
        return *in;
    }

    inline Pipe& operator++() {
        ++in;
        return *this;
    }

    template < typename TInput_ >
    bool process(TInput_ &textIn)
    {
        // 1. Generate Gapped Tuples
        TPipeTupler                                             tupler(textIn, limits);

        // 2. Sort Tuples by the first few characters
        TTupleComparator                                        _comp(TShape(), limits);
        TPoolSorter                                             sorter(tupler, _comp);

        // 3. Name tuples by their rank
        TPipeNamer                                              namer(sorter, _comp);

        // 4. Map text Positions to lexText positions
        TDislexMapper                                           _map(TShape::span, limits);
        TPoolMapper                                             mapper(namer, _map);
        mapper << namer;

        // 5. Discard positions, keep rank
        TPipeFilterI2                                           filter(mapper);

        // 6. Run SACA on lex text
        pool << filter;

        // 7. Reverse Transform is done during the reading process
        return true;
    }
};


// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Operator << for Pipes                                             [String]
// --------------------------------------------------------------------------

template < typename TInput, typename TShape, typename TSACA, typename TObject >
inline bool operator<<(Pipe< TInput, DislexExternal<TShape, TSACA> > &me, TObject &textIn)
{
    typedef Pipe< TInput, DislexExternal<TShape, TSACA> > TPipe;
    me.in = TPipe::TPipeReverseTransform(me.pool, TPipe::TDislexReverse(TShape::span, length(textIn)));
    return me.process(textIn);
}

// --------------------------------------------------------------------------
// Operator << for Pipes                                             [String]
// --------------------------------------------------------------------------

template < typename TInput, typename TShape, typename TSACA, typename TPair, typename TLimits, typename TObject >
inline bool operator<<(
    Pipe< TInput, Multi<DislexExternal<TShape, TSACA>, TPair, TLimits> > &me,
    TObject &textIn)
{
    typedef Pipe< TInput, Multi<DislexExternal<TShape, TSACA>, TPair, TLimits> > TPipe;
    me.limits = stringSetLimits(textIn);
    me.in = TPipe::TPipeReverseTransform(me.pool, TPipe::TDislexReverse(TShape::span, me.limits));
    return me.process(textIn);
}


// --------------------------------------------------------------------------
// function createGappedSuffixArray()                                [String]
// --------------------------------------------------------------------------

template < typename TSA, typename TText, typename TCyclicShape, typename TSACA>
inline void createGappedSuffixArray(
    TSA &SA, // must already be resized already
    TText const &s,
    TCyclicShape const & shape,
    ModCyclicShape<TCyclicShape> const &,
    DislexExternal<TCyclicShape, TSACA> const &)
{
    _createSuffixArrayPipelining(SA, s, DislexExternal<TCyclicShape, TSACA>());
}

// --------------------------------------------------------------------------
// function createGappedSuffixArray()                             [StringSet]
// --------------------------------------------------------------------------

template < typename TSA, typename TText, typename TSpec, typename TCyclicShape, typename TSACA>
inline void _createGappedSuffixArrayPipelining(
    TSA &SA, // must already be resized already
    StringSet<TText, TSpec> const &s,
    TCyclicShape const & shape,
    ModCyclicShape<TCyclicShape> const &,
    DislexExternal<TCyclicShape, TSACA> const &)
{
    // TODO
}


}
#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_DISLEX_EXTERNAL_H_