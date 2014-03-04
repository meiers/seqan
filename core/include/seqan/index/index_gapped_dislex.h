// ==========================================================================
//                          index_gapped_sa_dislex.h
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

#ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_
#define CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_

namespace SEQAN_NAMESPACE_MAIN
{



template <typename TString>
void printStr(TString const & str) { for(unsigned i=0; i< length(str); ++i)
    std::cout << str[i] << ","; std::cout << std::endl; }

template <typename T, typename A>
void printStrSet(StringSet<T,A> const & strSet) { for(unsigned i=0; i< length(strSet); ++i) {
    std::cout << "[" << i << "]"; printStr(strSet[i]); std::cout << std::endl; } }

template <typename TSACA = Skew7>
struct Dislex {};


// ==========================================================================
// DisLex Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction LexValue
// --------------------------------------------------------------------------

// default definition
template <typename TAlph, typename TMod>
struct LexValue {
    typedef unsigned int Type;
};

// TODO: smaller integer type if possible?


// --------------------------------------------------------------------------
// Metafunction LexText
// --------------------------------------------------------------------------

template <typename TText, typename TMod>
struct LexText {
    typedef typename Value<TText>::Type TAlph;
    typedef typename LexValue<TAlph, TMod>::Type TSize;

    typedef String<TSize> Type;
};

/*    template <typename TText, typename TSetSpec, typename TMod>
 struct LexText
 < StringSet<TText, TSetSpec>, TMod >
 {
 typedef StringSet<typename LexText<TText, TMod>::Type, Owner<ConcatDirect<> > > Type;
 };
 */



// ==========================================================================
// DisLex mapping functors
// ==========================================================================

// --------------------------------------------------------------------------
// Map: Text -> LexText                                              [String]
// --------------------------------------------------------------------------

// real dislex transformation:
// takes a position x and returns L(x)
template <
typename TInput,
typename TResult = TInput>
struct _dislexTransform :
public std::unary_function<TInput, TResult>
{
    TInput const S,N,B,C;

    _dislexTransform(TInput S_, TInput N_) : S(S_), N(N_), B(N_/S_), C(N_%S_)
    {}

    inline TResult operator() (TInput const & x) const
    {
        TInput r = x/S;
        TInput b = (N-x)%S;
        TInput ret = (S-1-b)*B + r;
        if (b <= C)
        ret += C-b;
        return static_cast<TResult>(ret);
    }
};


// --------------------------------------------------------------------------
// Map: Text -> LexText                                           [StringSet]
// --------------------------------------------------------------------------

// dislex transformation for multiple strings
// takes a Pair (seq,pos) and returns L(seq,pos) = pos
template <
typename TValue,
typename TString,
typename TResult = typename Value<TValue, 2>::Type>
struct _dislexTransformMulti :
public std::unary_function<TValue, TResult>
{
    TString const & limits;
    TResult const S;

    _dislexTransformMulti(TResult S_, TString const & stringSetLimits) : limits(stringSetLimits), S(S_)
    {}

    inline TResult operator() (const TValue & x) const
    {
        TResult seq = x.i1;
        TResult pos = x.i2;
        TResult N = limits[seq+1] - limits[seq];

        // TODO: Store the following values for each String?
        TResult B = N/S;
        TResult C = N%S;

        TResult r = (pos)/S;
        TResult b = (N-pos)%S;
        TResult ret = limits[seq] + (S-1-b)*B + r;
        if (b > C)  return ret;
        else        return ret + C-b;
    }
};


// --------------------------------------------------------------------------
// Map: LexText -> Text                                              [String]
// --------------------------------------------------------------------------

template <
    typename TInput,
    typename TResult = TInput>
struct _dislexReverseTransform :
public std::unary_function<TInput, TResult>
{
    const TInput  S,N,B,C;

    _dislexReverseTransform(TInput S_, TInput N_) : S(S_), N(N_), B(N_/S_), C(N_%S_)
    {}

    inline TResult operator() (const TInput & x) const
    {
        TInput b,r,ret;

        if (x < (S-C-1)*B) {
            b = S-1- x/B;
            r = x -(S-b-1)*B;
            ret = r*S + C - b + S ;
        } else {
            b = C -(x-(S-C-1)*B)/(B+1);
            r = x-(S-b-1)*B -C + b;
            ret = r*S + C - b;
        }
        return static_cast<TResult>(ret);
    }
};

// --------------------------------------------------------------------------
// Map: LexText -> Text                                           [StringSet]
// --------------------------------------------------------------------------

template <
    typename TInput,                    // global pos
    typename TLimits,                   // limits
    typename TResult = Pair<TInput> >   // local pos
struct _dislexReverseTransformMulti :
public std::unary_function<TInput, TResult>
{
    TLimits const & limits;
    TInput const S;

    _dislexReverseTransformMulti(TInput S_, TLimits const & stringSetLimits) : limits(stringSetLimits), S(S_)
    {}

    inline TResult operator() (const TInput & x) const
    {
        TResult ret;
        // binary search to find the corresponding sequence
        posLocalize(ret, x, limits);
        TInput N = limits[ret.i1+1] - limits[ret.i1];

        // TODO: Store the following values for each String?
        TInput B = N/S;
        TInput C = N%S;
        TInput b,r, i = ret.i2;
        if (i < (S-C-1)*B) {
            b = S-1- i/B;
            r = i -(S-b-1)*B;
            ret.i2 = r*S + C - b + S ;
        } else {
            b = C -(i-(S-C-1)*B)/(B+1);
            r = i-(S-b-1)*B -C + b;
            ret.i2 = r*S + C - b;
        }
        return ret;
    }
};




// ==========================================================================
// DisLex comparsion functors
// ==========================================================================

// NOTE:
//  - For a full suffix comparison use _SuffixLess (index_sa_gapped_qsort.h)
//  - For the comparison of suffixes of which we know are equal up to their
//    last character, use _ZeroBucketComparator (radix_inplace.h)
//  - Inside the _dislex() we need to compare k-mers. We use the
//    GappedSuffixQgramLess_ functors for that.
//  - The external algorithm saves a part of the sequence into a tuple and
//    then uses _dislexTupleComp/Multi to compare them. These have overloads
//    for bitpacked tuples.


// --------------------------------------------------------------------------
// struct GappedSuffixQgramLess_                                     [String]
// --------------------------------------------------------------------------

    template <typename TSAValue, typename TShape, typename TText, typename TResult=int>
    struct GappedSuffixQgramLess_;

    template <typename TSAValue, typename TShape, typename TText, typename TResult>
    struct GappedSuffixQgramLess_ : public std::binary_function<TSAValue, TSAValue, TResult>
    {
        typedef typename Size<TText>::Type                                          TSize;
        typedef ModifiedString<typename Suffix<TText const>::Type, ModCyclicShape<TShape> >      TSuffix;
        typedef typename Iterator<TSuffix,Standard>::Type                           TSuffIter;

        TText const &   _text;
        TShape const &  _shape;
        TSize const     _weight;

        GappedSuffixQgramLess_(TText const &text, TShape const & shape, TSize weight):
        _text(text), _shape(shape), _weight(weight)
        {}

        inline int operator() (TSAValue a, TSAValue b) const
        {
            if (a == b) return 0;

            TSuffix sa(suffix(_text, a), _shape);
            TSuffix sb(suffix(_text, b), _shape);

            TSuffIter saIt = begin(sa, Standard());
            TSuffIter sbIt = begin(sb, Standard());

            TSuffIter saEnd = end(sa, Standard());
            TSuffIter sbEnd = end(sb, Standard());

            TSize p = 0;

            for (; saIt < saEnd && sbIt < sbEnd && p < _weight; ++saIt, ++sbIt, ++p)
            {
                if (*saIt < *sbIt) return -1;
                if (*saIt > *sbIt) return 1;
            }

            if (!(saIt < saEnd)) return -1;
            if (!(sbIt < sbEnd)) return 1;

            return 0;
        }
    };

// --------------------------------------------------------------------------
// struct GappedSuffixQgramLess_                                  [StringSet]
// --------------------------------------------------------------------------

    template <typename TSAValue, typename TShape, typename TText, typename TSpec, typename TResult>
    struct GappedSuffixQgramLess_ <TSAValue, TShape, StringSet<TText, TSpec>, TResult> :
        public std::binary_function<TSAValue, TSAValue, TResult>
    {
        typedef StringSet<TText, TSpec>                         TSet;
        typedef typename Size<TText>::Type                      TSize;
        typedef ModifiedString<typename Suffix<TText const>::Type, ModCyclicShape<TShape> >      TSuffix;
        typedef typename Iterator<TSuffix,Standard>::Type       TSuffIter;

        TSet const &   _text;
        TShape const &  _shape;
        TSize const     _weight;

        GappedSuffixQgramLess_(TSet const &text, TShape const & shape, TSize weight):
        _text(text), _shape(shape), _weight(weight)
        {}

        inline int operator() (TSAValue a, TSAValue b) const
        {
            if (a == b) return 0;

            TSuffix sa(suffix(_text, a), _shape);
            TSuffix sb(suffix(_text, b), _shape);

            TSuffIter saIt = begin(sa, Standard());
            TSuffIter sbIt = begin(sb, Standard());

            TSuffIter saEnd = end(sa, Standard());
            TSuffIter sbEnd = end(sb, Standard());

            TSize p = 0;

            for (; saIt < saEnd && sbIt < sbEnd && p < _weight; ++saIt, ++sbIt, ++p)
            {
                if (*saIt < *sbIt) return -1;
                if (*saIt > *sbIt) return 1;
            }

            // if both suffixes are not yet empty, the they are equal
            if (saIt < saEnd && sbIt < sbEnd)
                return 0;

            // if both suffixes are empty, the seq. id and then the underlying suffix length decides
            if (!(saIt < saEnd) && !(sbIt < sbEnd))
            {
                if (getSeqNo(a) > getSeqNo(b)) return -1;
                if (getSeqNo(a) < getSeqNo(b)) return 1;

                if (getSeqOffset(a) > getSeqOffset(b)) return -1;
                if (getSeqOffset(a) < getSeqOffset(b)) return 1;

                // Does not occur
                SEQAN_ASSERT_EQ(true,false);
                return 0;
            }

            // only one suffix is empty
            if (!(saIt < saEnd)) return -1;
            if (!(sbIt < sbEnd)) return 1;

            // Does not occur
            SEQAN_ASSERT_EQ(true, false);
            return 0;
        }
    };

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// function _dislex()                                                [String]
// --------------------------------------------------------------------------

template <typename TLexText, typename TSA, typename TText, typename TCyclicShape>
inline void _dislex(
    TLexText & lexText,                         // random access (write)
    TSA const & partiallyOrderedSA,             // sequential scan
    TText const & origText,                     // random access
    TCyclicShape const & cyclic)
{
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Value<TLexText>::Type                  TRank;
    typedef ModifiedString<typename Suffix<TText const
    >::Type, ModCyclicShape<TCyclicShape> >      TModText;

    if (length(partiallyOrderedSA) < 1) return; // otherwise *sa will fail

    // dislex position calculator
    _dislexTransform<TSize> dislex(cyclic.span, length(origText));

    // q-gram comparator to determine rank
    GappedSuffixQgramLess_<TSAValue, TCyclicShape,TText> comp(origText, cyclic, weight(cyclic));

    resize(lexText, length(origText), Exact());

    TSAIter sa    = begin(partiallyOrderedSA, Standard());
    TSAIter saEnd = end(partiallyOrderedSA, Standard());
    TRank       rank    = 0;
    TSAValue    txtPos  = *sa++;

    // scan along the SA
    for(; sa < saEnd; txtPos = *sa++)
    {
        // write rank to position in lexText
        lexText[dislex(txtPos)] = rank;

        // compare two consecutive values: this is probably slow
        if(comp(txtPos, *sa))
            ++rank;

        //std::cout << txtPos << "  ...   " << *sa << " -> " << comp(txtPos, *sa) << std::endl;
        SEQAN_ASSERT_GEQ(0, comp(txtPos,*sa));
    }
    lexText[dislex(txtPos)] = rank;
}


// --------------------------------------------------------------------------
// function _dislex()                                             [StringSet]
// --------------------------------------------------------------------------

template < typename TLexText, typename TSA, typename TText, typename TTextSpec, typename TCyclicShape>
inline void _dislex(
    TLexText & lexText,                             // random access
    TSA const & partiallyOrderedSA,                 // sequential scan
    StringSet<TText, TTextSpec> const & origText,   // random access
    TCyclicShape const & cyclic)
{
    typedef typename Size<TSA>::Type                                    TSize;
    typedef typename Value<TSA>::Type                                   TSAValue;
    typedef typename Value<typename Concatenator<TLexText>::Type>::Type TRank;
    typedef typename Suffix<StringSet<TText, TTextSpec> const>::Type    TSuffix; // NOTE: This is != Suffix<TText>
    typedef ModifiedString<TSuffix, ModCyclicShape<TCyclicShape> >      TModText;
    typedef typename Iterator<TSA const, Standard>::Type                TSAIter;

    // position calculator
    typedef StringSet<TText, TTextSpec> const               TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type      TStringSetLimits;   // expected: String<unsigned>
    typedef Pair<typename Size<TText>::Type>                TSetPosition;       // expected: Pair<unsigned>

    if (length(partiallyOrderedSA) < 1) return; // otherwise *sa will fail
    
    TStringSetLimits xxx = stringSetLimits(origText);
    _dislexTransformMulti<TSetPosition, TStringSetLimits>
    dislex(cyclic.span, xxx);


    // q-gram comparator to determine rank
    GappedSuffixQgramLess_<TSAValue, TCyclicShape, StringSet<TText, TTextSpec> const> comp(origText, cyclic, static_cast<TSize>(weight(cyclic)));

    resize(lexText, lengthSum(origText));

    // scan along the SA
    TSAIter sa    = begin(partiallyOrderedSA, Standard());
    TSAIter saEnd = end(partiallyOrderedSA, Standard());

    TRank       rank = 0;
    TSAValue    txtPos = *sa++;


    for(; sa < saEnd; txtPos = *sa++)
    {
        // write rank to position in lexText
        lexText[dislex(txtPos)] = rank;

        // compare two consecutive values, this is probably slow
        if(comp(txtPos, *sa))
            ++rank;

        //std::cout << txtPos << "  ...   " << *sa << " -> " << comp(txtPos, *sa) << std::endl;
        SEQAN_ASSERT_GEQ(0, comp(txtPos,*sa));
    }
    lexText[dislex(txtPos)] = rank;
}


// --------------------------------------------------------------------------
// function _dislexReverse()                                         [String]
// --------------------------------------------------------------------------

template <typename TSA, typename TCyclicShape, typename TText, typename TLexSA>
void _dislexReverse(
    TSA & finalSA,                                  // random access
    TLexSA const & lexSA,                           // sequential scan
    TText const &,                                  // not needed
    TCyclicShape const & cyclic)
{
    typedef typename Iterator<TSA const, Standard>::Type    TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Size<TSA>::Type                        TSize;

    _dislexReverseTransform<TSize> dislexRev(cyclic.span, length(lexSA));

    TLexSAIter sa       = begin(lexSA, Standard());
    TLexSAIter saEnd    = end(lexSA, Standard());
    TSAIter insert      = begin(finalSA, Standard());

    for(; sa < saEnd; ++sa, ++insert)
    *insert = dislexRev (*sa);
}

// --------------------------------------------------------------------------
// function _dislexReverse()                                      [StringSet]
// --------------------------------------------------------------------------

template <typename TSA, typename TLexSA, typename TCyclicShape, typename TText, typename TTextSpec>
void _dislexReverse(
    TSA & finalSA,                                  // random access
    TLexSA const & lexSA,                           // sequential scan
    StringSet<TText, TTextSpec> const & origText,
    TCyclicShape const & cyclic)
{
    typedef typename Iterator<TLexSA const, Standard>::Type TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Size<TSA>::Type                        TSize;

    typedef StringSet<TText, TTextSpec> const               TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type      TStringSetLimits;   // expected: String<unsigned>

    TStringSetLimits xxx = stringSetLimits(origText);
    _dislexReverseTransformMulti<TSize, TStringSetLimits>
    dislexRev(cyclic.span, xxx);

    TLexSAIter sa       = begin(lexSA, Standard());
    TLexSAIter saEnd    = end(lexSA, Standard());
    TSAIter insert      = begin(finalSA, Standard());

    for(; sa < saEnd; ++sa, ++insert)
    *insert = dislexRev (*sa);
}

// --------------------------------------------------------------------------
// function createGappedSuffixArray()                                [Dislex]
// --------------------------------------------------------------------------
// Only defined for CyclicShapes

template < typename TSA, typename TText, typename TCyclicShape, typename TSACA>
inline void createGappedSuffixArray(
    TSA &SA, // must already be resized already
    TText const &s,
    TCyclicShape const & shape,
    ModCyclicShape<TCyclicShape> const &,
    Dislex<TSACA> const &)
{
    typedef typename LexText<TText, TCyclicShape>::Type         TLexText;

    // if alph too big, problem with counter array!
    SEQAN_ASSERT_GEQ(256u, valueSize<typename Value<TText>::Type>());

    double teim = sysTime();

    // insert positions into SA
    _initializeSA(SA, s);

    //if (length(SA) < 1) return;

    std::cout << "   |     init: " << sysTime() - teim << "s" << std::endl; teim = sysTime();

    // sort newSA according to the Shape
    inplaceRadixSort(SA, s, weight(shape)+1, shape, ModCyclicShape<TCyclicShape>());

    std::cout << "   | radix[" << (int)weight(shape) << "]: " << sysTime() - teim << "s" << std::endl; teim = sysTime();

    // disLexTransformation
    TLexText lexText;
    _dislex(lexText, SA, s, shape);

    std::cout << "   |   dislex: " << sysTime() - teim << "s" << std::endl; teim = sysTime();

    // Build Index using Skew7
    Index<TLexText, IndexSa<> > normalIndex(lexText);
    indexCreate(normalIndex, EsaSA(), TSACA());

    std::cout << "   |     saca: " << sysTime() - teim << "s (len = " << length(concat(lexText)) << ")" << std::endl; teim = sysTime();


    // reverse Transform of Index:
    _dislexReverse(SA, indexSA(normalIndex), s, shape);

    std::cout << "   |  reverse: " << sysTime() - teim << "s (len = " << length(indexSA(normalIndex)) << ")" << std::endl; teim = sysTime();

}












// ==========================================================================
// Garbage
// ==========================================================================

template <typename TSA,
typename TSuffixModifier,
typename TText, typename TTextSpec>
void dislexReverse___old(
     TSA & finalSA,                              // random access
     TSA const & lexSA,                          // sequential scan
     StringSet<TText, TTextSpec> const & text,   // only needed for lengths
     TSuffixModifier const & cyclic)
{
    typedef typename Iterator<TSA const, Standard>::Type    TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;

    TSize   S = cyclic.span;
    TSize   M = length(text);

    TLexSAIter sa       = begin(lexSA, Standard());
    TLexSAIter saEnd    = end(lexSA, Standard());
    TSAIter insert      = begin(finalSA, Standard());
    TSize blockNum, blockRank;

    // artificial border characters get sorted to the first S*M SA entries. Ignore them.
    for(sa += S*M; sa < saEnd; ++sa, ++insert)
    {
        TSize seq = getSeqNo(*sa);
        TSize   p = getSeqOffset(*sa);
        TSize   N = length(text[seq]);
        TSize   B = N / S + 1;      // block size
        TSize   E = N % S;          // number of blocks of size B+1

        blockNum  = p > E*(B+1) ? E + (p - E*(B+1))/B : p/(B+1);
        blockRank = p > E*(B+1) ? (p - E*(B+1)) % B : p % (B+1);
        *insert = TSAValue(seq, blockRank*S + blockNum);

        //std::cout << "lexSA[" << sa - begin(lexSA) << "] = " << *sa << "\t=> num: " << blockNum << ", rank: " << blockRank << std::endl;
        //if(sa - begin(lexSA) > 100) break;
    }
}





template < typename TLexText,
typename TSA,
typename TText,      typename TTextSpec,
typename TSuffixModifier>
inline void dislexTransform___old(
                                  StringSet<TLexText, Owner<ConcatDirect<> > > & lexText,        // random access
                                  TSA const & partiallyOrderedSA,                 // sequential scan
                                  StringSet<TText, TTextSpec> const & origText,   // random access
                                  TSuffixModifier const & cyclic)
{
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Value<typename Concatenator<TLexText>::Type>::Type TRank;
    typedef ModifiedString<typename Suffix<TText const
    >::Type, ModCyclicShape<TSuffixModifier> >          TModText;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;

    TSize   W = weight(cyclic);         // weight for comparison of consecutive entries
    TSize   S = cyclic.span;            // shape length
    TSize   M = length(origText);
    TSize   overall = 0;

    // allocate lexText
    String<TSize> lengths;
    resize(lengths,M);
    clear(lexText);
    for (TSize i = 0; i < M; ++i)
    {
        lengths[i] = length(origText[i]);
        overall += lengths[i] + S;
        appendValue(stringSetLimits(lexText), overall);
    }
    resize(concat(lexText), overall);


    // scan along the SA
    TSAIter sa    = begin(partiallyOrderedSA, Standard());
    TSAIter saEnd = end(partiallyOrderedSA, Standard());

    TRank       rank = S;        // rank 0..S-1 are reserved for delimiters
    TSAValue    txtPos = *sa++;
    TSAValue    lexPos;
    TSize seq   =  getSeqNo(txtPos);
    TSize pos    = getSeqOffset(txtPos);
    TSize N, B, E;
    for(; sa < saEnd; txtPos = *sa++)
    {
        N = lengths[seq];           // text length
        B = N/S + 1;                // block length
        E = N % S;                  // number of blocks of length B+1

        // determine insert position in lexText
        pos = getSeqOffset(txtPos);
        seq = getSeqNo(txtPos);
        lexPos = TSAValue(seq, (pos%S)*B + pos/S + (pos%S < E ? pos%S : E));
        concat(lexText)[ posGlobalize(lexPos, stringSetLimits(origText)) ] = rank;

        // compare two consecutive values, this is probably slow
        TModText modStr1(suffix(origText, txtPos), cyclic);
        TModText modStr2(suffix(origText, *sa), cyclic);

        // if k-mer changes, increase the rank by 1
        if(prefix(modStr1, std::min(length(modStr1),W))
           != prefix(modStr2, std::min(length(modStr2),W)))
        ++rank;

        SEQAN_ASSERT_LEQ(prefix(modStr1, std::min(length(modStr1),W)),
                         prefix(modStr2, std::min(length(modStr2),W)));
    }
    N = lengths[seq];
    B = N/S + 1;
    E = N % S;
    pos = getSeqOffset(txtPos);
    seq = getSeqNo(txtPos);
    lexPos = TSAValue(seq, (pos%S)*B + pos/S + (pos%S < E ? pos%S : E));
    concat(lexText)[ posGlobalize(lexPos, stringSetLimits(origText)) ] = rank;


    // manually insert delimiters between blocks:
    for (seq = 0; seq < M; ++seq)
    {
        N = lengths[seq];
        B = N/S + 1;
        E = N % S;
        TSize insertPos = 0;
        for(TSize p = 0; p < S; ++p)
        {
            insertPos += B + (p<E ? 1:0);
            concat(lexText)[ posGlobalize(insertPos-1, stringSetLimits(origText)) ] = S-p-1;
        }
        SEQAN_ASSERT_EQ(lengths[seq]+S, insertPos);
    }

    //        for(unsigned i=0; i< 100; ++i) std::cout << i << "\t" << lexText[0][i] << std::endl;
}


template <typename TSA,
typename TSuffixModifier,
typename TText>
void dislexReverse___old(
                         TSA & finalSA,                              // random access
                         TSA const & lexSA,                          // sequential scan
                         TText const &,                              // not needed
                         TSuffixModifier const & cyclic)
{
    typedef typename Iterator<TSA const, Standard>::Type    TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Size<TSA>::Type                        TSize;

    TSize   S = cyclic.span;
    TSize   L = length(lexSA);
    TSize   B = L / S;  // block size
    TSize   E = L % S;  // number of blocks of size B+1

    TLexSAIter sa       = begin(lexSA, Standard());
    TLexSAIter saEnd    = end(lexSA, Standard());
    TSAIter insert      = begin(finalSA, Standard());
    TSize blockNum, blockRank;

    // artificial border characters get sorted to the first S SA entries. Ignore them.
    for(sa+=S; sa < saEnd; ++sa, ++insert)
    {
        blockNum  = *sa > E*(B+1) ? E + (*sa - E*(B+1))/B : *sa/(B+1);
        blockRank = *sa > E*(B+1) ? (*sa - E*(B+1)) % B : *sa % (B+1);
        *insert = blockRank*S + blockNum;
    }
}




template < typename TLexText,
typename TSA,
typename TText,
typename TSuffixModifier>
inline void dislexTransform___old(
                                  TLexText & lexText,                         // random access (write)
                                  TSA const & partiallyOrderedSA,             // sequential scan
                                  TText const & origText,                     // random access
                                  TSuffixModifier const & cyclic)
{
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Value<TLexText>::Type                  TRank;
    typedef ModifiedString<typename Suffix<TText const
    >::Type, ModCyclicShape<TSuffixModifier> >          TModText;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;

    TSize   W = weight(cyclic);         // weight for comparison of consecutive entries
    TSize   S = cyclic.span;            // shape length
    TSize   N = length(origText);       // length of text = length of SA
    TSize   B = N/S + 1;                // block length
    TSize   E = N % S;                  // number of blocks of length B+1
    TSize   L = N+S; // = S*B + E       // length of lexText

    // allocate lexText
    resize(lexText, L, Exact());

    // scan along the SA
    TSAIter sa    = begin(partiallyOrderedSA, Standard());
    TSAIter saEnd = end(partiallyOrderedSA, Standard());

    TRank       rank    = static_cast<TRank>(S);        // rank 0..S-1 are reserved for delimiters
    TSAValue    txtPos  = *sa++;
    TSize       lexPos;
    for(; sa < saEnd; txtPos = *sa++)
    {
        // determine insert position in lexText
        lexPos = (txtPos%S)*B + txtPos/S + (txtPos%S < E ? txtPos%S : E);
        lexText[lexPos] = rank;

        // compare two consecutive values, this is probably slow
        TModText modStr1(suffix(origText, txtPos), cyclic);
        TModText modStr2(suffix(origText, *sa), cyclic);

        // if k-mer changes, increase the rank by 1
        if(prefix(modStr1, std::min(length(modStr1),W))
           != prefix(modStr2, std::min(length(modStr2),W)))
        ++rank;

        SEQAN_ASSERT_LEQ(   prefix(modStr1, std::min(length(modStr1),W)),
                         prefix(modStr2, std::min(length(modStr2),W)));
    }
    lexPos = (txtPos%S)*B + txtPos/S + (txtPos%S<E ? txtPos%S : E);
    lexText[lexPos] = rank;

    // manually insert delimiters between blocks:
    lexPos  = 0; // now 1 right of insert pos
    rank    = S-1;
    for(TSize p = 0; p < S; --rank, ++p)
    {
        lexPos += B + (p<E ? 1:0);
        lexText[lexPos-1] = S-p-1;
    }
    SEQAN_ASSERT_EQ(L, lexPos);
}


    /*
     
    // Old GappedTupler, but too complicated since tupleLen is user determined
     

    // --------------------------------------------------------------------------
    // Pipe < TInput, GappedTupler >                                   [Sequence]
    // --------------------------------------------------------------------------

    template <
    typename TInput,
    typename TShape,
    unsigned tupleLen,
    bool omitLast,
    typename TPack>
    struct Pipe< TInput, GappedTupler<TShape, tupleLen, omitLast, TPack> >
    {
        typedef typename Value<Pipe>::Type          TOutput;
        typedef typename Value<TOutput, 2 >::Type	TTuple;
        typedef typename Value<TTuple>::Type		TValue;
        typedef typename Size<TInput>::Type         TSize;

        // BuffSize gets rounded up to a multiple of span
        enum { BuffSize = (WEIGHT<TShape>::VALUE + tupleLen - 1)/WEIGHT<TShape>::VALUE * TShape::span };

        TInput      &in;
        TOutput     tmp;
        TSize       lastTuples;

        TValue      buffer[BuffSize];   // all elements will be shifted in ++
        // (ring buffer is complicated due to many if(p > TShape::span)
        // queries, maybe I will try that later)
        TSize       carePos[tupleLen];

        Pipe(TInput& _in): in(_in), buffer()
        {
            // TODO(meiers): These care positions of the shape are known at compile time
            //       They should be computed at compile time
            carePos[0] = TShape::loffset;
            for(TSize i=1, j=0; i<tupleLen; ++i, ++j)
            {
                if(j== WEIGHT<TShape>::VALUE) j=0;
                carePos[i] = carePos[i-1] + TShape::diffs[j];
            }
        }

        inline TOutput const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            if (eof(in)) --lastTuples;

            // it's just a jump to the left
            memmove(buffer, buffer+1, (BuffSize-1)*sizeof(TValue) );
            // memmove is probably faster than unrolled loop:
            // Loop<ShiftLeftWorker2_, BuffSize - 1>::run(this->buffer);

            if (lastTuples < TuplerNumberOfLastTuples_<tupleLen, omitLast>::VALUE)
                buffer[BuffSize - 1] = TValue();
            else
            {
                buffer[BuffSize - 1] = *in;
                ++in;
            }

            ++tmp.i1;
            _fillTmp2(); // a bit expensive, but easier to implement
            return *this;
        }


        inline void fill() {

            unsigned i;
            for(i = 0; i < BuffSize && !eof(in); ++i, ++in)
                buffer[i] = *in;

            // set lastTuples depending on the omitLast flag
            if (TuplerNumberOfLastTuples_<tupleLen, omitLast>::VALUE > BuffSize - i)
                lastTuples = TuplerNumberOfLastTuples_<tupleLen, omitLast>::VALUE - (BuffSize - i);
            else
                lastTuples = 0; // this will cause eof() of this pipe

            // fill remaining buffer, if it hasn't been filled yet
            for (; i < BuffSize; ++i)
                buffer[i] = TValue();

            // fill tmp
            tmp.i1 = 0;
            _fillTmp2();
        }

        inline void _fillTmp2()
        {
            // TODO: Use Loop struct?
            for(unsigned i = 0; i < tupleLen; ++i)
                tmp.i2[i] = buffer[carePos[i]];
        }
    };


    // --------------------------------------------------------------------------
    // Pipe < TInput, GappedTupler >                                      [Multi]
    // --------------------------------------------------------------------------


    template <
    typename TInput,
    typename TShape,
    unsigned tupleLen,
    bool omitLast,
    typename TPack,
    typename TPair,
    typename TLimitsString >
    struct Pipe< TInput, Multi<GappedTupler<TShape, tupleLen, omitLast, TPack>, TPair, TLimitsString> >
    {
        typedef typename Value<Pipe>::Type                              TOutput;
        typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
        typedef typename Value<TTuple>::Type							TValue;
        typedef typename Size<TInput>::Type                             TSize;
        typedef PairIncrementer_<TPair, TLimitsString>                  Incrementer;


        // BuffSize gets rounded up to a multiple of span
        enum { BuffSize = (WEIGHT<TShape>::VALUE + tupleLen - 1)/WEIGHT<TShape>::VALUE * TShape::span };

        TInput                      &in;
        Incrementer					localPos;
        TOutput                     tmp;
        TSize                       seqLength, lastTuples;
        TLimitsString const         &limits;

        TValue                      buffer[BuffSize];
        TSize                       carePos[tupleLen];

        template <typename TLimitsString_>
        // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
        Pipe(TInput& _in, TLimitsString_ &_limits):  in(_in), limits(_limits)
        {
            // TODO: These care positions of the shape are known at compile time...
            carePos[0] = TShape::loffset;
            for(TSize i=1, j=0; i<tupleLen; ++i, ++j)
            {
                if(j== WEIGHT<TShape>::VALUE) j=0;
                carePos[i] = carePos[i-1] + TShape::diffs[j];
            }
            std::cout << TuplerNumberOfLastTuples_<tupleLen, omitLast>::VALUE << std::endl;
        }

        inline TOutput const & operator*() const
        {
            //std::cout << "lastTuple:" << lastTuples << "  localPos:" << static_cast<TPair>(localPos) << "  tmp:" << tmp.i2 << "  eos:" << this->eos() << "  eof:" << std::endl;
            return tmp;
        }

        inline Pipe& operator++()
        {
            // process next sequence
            if (eos())
                if (--lastTuples == 0)
                {
                    assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
                    fill();
                    return *this;
                }

            // shift left 1 character
            memmove(buffer, buffer+1, (BuffSize-1)*sizeof(TValue) );
            assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);

            if (lastTuples < TuplerNumberOfLastTuples_<tupleLen, omitLast>::VALUE)
            {
                buffer[BuffSize - 1] = TValue();
            } else
            {
                buffer[BuffSize - 1] = *in;
                ++localPos;
                ++in;
            }

            _fillTmp2();
            return *this;
        }

        inline void fill()
        {
            do {
                unsigned i = 0;
                if (!eof(in))
                    do {
                        buffer[i] = *in;
                        ++in;
                        ++i;
                        ++localPos;
                    } while ((i < BuffSize) && !eos());
                lastTuples = TuplerNumberOfLastTuples_<tupleLen, omitLast>::VALUE;

                // eventually, reduce the number of half-filled tuples
                if (lastTuples <= BuffSize - i)
                    lastTuples = 0;
                else
                {
                    lastTuples -= BuffSize - i;
                    
                    // fill up with null chars
                    for(; i < BuffSize; ++i)
                        tmp.i2.i[i] = TValue();
                }
                
                if (lastTuples == 0)
                    assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
                
            } while ((lastTuples == 0) && !eof(in));
            
            assignValueI2(tmp.i1, 0);
            _fillTmp2();
        }
        
        inline bool eos() const
        {
            return (getValueI1(static_cast<TPair>(localPos)) > 0) && (getValueI2(static_cast<TPair>(localPos)) == 0);
        }
        
        inline void _fillTmp2()
        {
            // TODO: Use Loop struct?
            for(unsigned i = 0; i < tupleLen; ++i)
                tmp.i2[i] = buffer[carePos[i]];
        }
    };
    
    */
    

}

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_
