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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================
// The gapped suffix array as a subclass of IndexSa
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_GENREALIZED_SA_H
#define SEQAN_HEADER_INDEX_GENREALIZED_SA_H


namespace SEQAN_NAMESPACE_MAIN
{
    
// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct IndexSa;
    
// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(meiers): Remove this definition as soon as Ticket #1096 has been solved: http://trac.seqan.de/ticket/1096
template < typename TText, typename TSpec >
struct DefaultFinder< Index<TText, IndexSa<TSpec> > > {
    typedef EsaFindMlr Type;	// standard suffix array finder is mlr-heuristic
};

template <typename TSuffixMod, typename TSpec = void>
struct Generalized {};


template <typename TText, typename TSuffixMod, typename TSpec>
class Index<TText, IndexSa< Generalized<TSuffixMod, TSpec > > > :
    public Index<TText, IndexSa<TSpec> >
{
public:
    typedef Index<TText, IndexSa<TSpec> >            TSuper;
    typedef typename Suffix<Index>::Type             TSuffix;
    typedef typename Cargo<TSuffix>::Type            TModifierCargo;
    
    // derive           text, cargo, sa
    TModifierCargo      suffMod;
    
    Index() {}
    
    Index(Index & other) : suffMod(other.suffMod)
    {} //automatic call of TSuper's copy constructor

    Index(Index const & other) : suffMod(other.suffMod)
    {}

    template <typename TText_>
    Index(TText_ & _text) : TSuper(_text)
    {}
    
    template <typename TText_>
    Index(TText_ const & _text) : TSuper(_text)
    {}

    template <typename TText_>
    Index(TText_ & _text, TModifierCargo const & suffMod) : TSuper(_text), suffMod(suffMod)
    {}

    template <typename TText_>
    Index(TText_ const & _text, TModifierCargo const & suffMod) : TSuper(_text), suffMod(suffMod)
    {}
};



// ============================================================================
// Metafunctions
// ============================================================================
    
// ----------------------------------------------------------------------------
// DefaultIndexCreator
// ----------------------------------------------------------------------------

// TODO(meiers): Update to a better method
// TODO(meiers): Specialise this for CyclicShape
template < typename TText, typename TSuffixMod, typename TSpec>
struct DefaultIndexCreator<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > >, FibreSA> {
    typedef SAQSort Type;
};
    

// ----------------------------------------------------------------------------
// Metafunction Suffix                                                  general
// ----------------------------------------------------------------------------

// general modified suffix
template <typename TText, typename TSuffixMod, typename TSpec>
struct Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > >
{
    typedef ModifiedString<typename Suffix<TText>::Type, TSuffixMod>    Type;
};

// general modified suffix; const variant
template <typename TText, typename TSuffixMod, typename TSpec>
struct Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > const>
{
    typedef ModifiedString<typename Suffix<TText const>::Type, TSuffixMod >     Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function suffixModifier
// ----------------------------------------------------------------------------
    
template <typename TText, typename TSuffixMod, typename TSpec>
inline typename Reference<typename Cargo<typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > >::Type>::Type>::Type
suffixModifier(Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > &index)
{
    return index.suffMod;
}
// TODO(meiers): const version?

// ----------------------------------------------------------------------------
// Function suffix()                                                    general
// ----------------------------------------------------------------------------

template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin>
inline typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > >::Type
suffix(Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > & t, TPosBegin pos_begin)
{
    typedef typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > >::Type TModifiedSuffix;
    return TModifiedSuffix(suffix(indexText(t), pos_begin), suffixModifier(t));
}

// const variant
template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin>
inline typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > const>::Type
suffix(Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > const &t, TPosBegin pos_begin)
{
    typedef typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > const>::Type TModifiedSuffix;
    return TModifiedSuffix(suffix(indexText(t), pos_begin), suffixModifier(t));
}

// ----------------------------------------------------------------------------
// Function infix()                                                     general
// ----------------------------------------------------------------------------

// gapped Index SA
template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Prefix<typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > >::Type>::Type
infix(Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > &t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(t, pos_begin), pos_end - pos_begin);
}

// gapped Index SA; const variant
template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Prefix<typename Suffix<Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > const>::Type>::Type
infix(Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > const &t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(t, pos_begin), pos_end - pos_begin);
}

// ----------------------------------------------------------------------------
// Function indexCreate (FibreSA)                                       general
// ----------------------------------------------------------------------------

template <typename TText, typename TSuffixMod, typename TSpec, typename TSpecAlg>
inline bool indexCreate(Index<TText, IndexSa<Generalized<TSuffixMod, TSpec> > > &t, FibreSA, TSpecAlg const alg)
{
    resize(indexSA(t), length(indexRawText(t)), Exact());
    createGappedSuffixArray(indexSA(t), indexText(t), suffixModifier(t), alg);
    return true;
}

    
/*

    

// ----------------------------------------------------------------------------
// Function _findFirstIndex()                                    used in find()
// ----------------------------------------------------------------------------

template < typename TText, typename TCycShape, typename TSpec, typename TSpecFinder, typename TPattern >
inline void
_findFirstIndex(Finder< Index<TText, IndexSa<Generalized<TCycShape, TSpec> > >, TSpecFinder > &finder,
                TPattern const &pattern,
                EsaFindMlr const)
{
    typedef Index<TText, IndexSa<Generalized<TCycShape, TSpec> > > TIndex;
    TIndex &index = haystack(finder);
    indexRequire(index, EsaSA());
    finder.range = equalRangeSAIterator(index, indexSA(index), pattern); // handle the index, not the text!
}

// TODO(sascha): _findFirstIndex for Lcp and STree variant, see index_find_esa.h







///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             experimental stuff                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // TODO(sascha): Whole BottomUp iterator!
    // TODO(sascha): Missing functions for Iterator, e.g. repLength
    

// copy 'n paste from index_sa_stree.h
template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it,
                          TString const & pattern, TSize & lcp)
{
    typedef Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >   TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Size<TSA const>::Type                              TSASize;
    typedef typename Iterator<TSA const, Standard>::Type                TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>                   TSearchTreeIterator;
    
    if (_isLeaf(it, HideEmptyEdges()))
        return false;
    
    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    // TText const & text = indexText(index);                                        // CHANGE(sascha): handle index instead of text
        
    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);
    Pair<TSAIterator> range = _equalRangeSA(index, node, pattern, value(it).repLen); // CHANGE(sascha): handle index instead of text
    
    if (range.i1 >= range.i2)
        return false;
    
    // Save vertex descriptor.
    _historyPush(it);
    
    // Update range, lastChar and repLen.
    value(it).range.i1 = range.i1 - begin(sa, Standard());
    value(it).range.i2 = range.i2 - begin(sa, Standard());
    value(it).lastChar = back(pattern);
    value(it).repLen += length(pattern);
        
    lcp = length(pattern);
    
    return true;
}


template < typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec >
inline typename Prefix<typename Suffix<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > >::Type>::Type
representative(Iter<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >, VSTree<TSpec> > &it)
{
    return prefix(suffix(container(it), getOccurrence(it)), repLength(it));
}

template < typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec >
inline typename Prefix<typename Suffix<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > const>::Type>::Type
representative(Iter<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > const, VSTree<TSpec> > const &it)
{
    return prefix(suffix(container(it), getOccurrence(it)), repLength(it));
}


// seperate specializations for String and StringSet needed to avoid ambiguities
template <typename TPos, typename TText, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > const>::Type
suffixLength(TPos pos, Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > const &index)
{
    return length(suffix(index, pos));
}
    
template <typename TPos, typename TString, typename TSSetSpec, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<StringSet<TString, TSSetSpec>, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > const>::Type
suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > const &index)
{
    return length(suffix(index, pos));
}

template <typename TPos, typename TText, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > >::Type
suffixLength(TPos pos, Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > &index)
{
    return length(suffix(index, pos));
}
template <typename TPos, typename TString, typename TSSetSpec, typename TShapeSpec, typename TIndexSpec>
inline typename Size<Index<StringSet<TString, TSSetSpec>, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > >::Type
suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, IndexSa<Generalized<TShapeSpec, TIndexSpec> > > &index)
{
    return length(suffix(index, pos));
}


template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;
    
    // TODO(esiragusa): use HideEmptyEdges()
    if (_isLeaf(it, HideEmptyEdges()))
        return false;
    
    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    
    // TODO(esiragusa): check nodeHullPredicate
    
    Pair<typename Size<TIndex>::Type> saRange = range(it);
    
    // TODO(esiragusa): remove this check.
    if (saRange.i1 >= saRange.i2) return false;
    
    // Skip $-edges.
    while (suffixLength(saAt(saRange.i1, index), index) <= value(it).repLen)
    {
        // TODO(esiragusa): remove this check and ++saRange.i1 in loop.
        // Interval contains only $-edges.
        if (++saRange.i1 >= saRange.i2)
            return false;
    }
    
    // Get first and last characters in interval.
    // CHANGED(sascha): get char via the suffix instead of textAt
    TAlphabet cLeft  = value(suffix(index, saAt(saRange.i1,     index)), value(it).repLen);
    TAlphabet cRight = value(suffix(index, saAt(saRange.i2 - 1, index)), value(it).repLen);
    
    // Save vertex descriptor.
    _historyPush(it);
    
    // Update left range.
    value(it).range.i1 = saRange.i1;
    
    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);
        
        // CHANGED(sascha): handle index instead of text
        TSAIterator upperBound = _upperBoundSA(index, node, cLeft, value(it).repLen);
        
        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }
    
    // Update child repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;
        
    return true;
}
    
  
    
    
template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goRight(Iter<Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<Generalized<TShapeSpec, TIndexSpec> > >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;
            
    if (isRoot(it))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
            
    Pair<typename Size<TIndex>::Type> saRange;
    saRange.i1 = value(it).range.i2;
    saRange.i2 = (_isSizeInval(value(it).parentRight)) ? length(sa) : value(it).parentRight;

    if (saRange.i1 >= saRange.i2) return false;

    // Change repLen to parent repLen.
    value(it).repLen--;

    // TODO(esiragusa): don't check for empty edges (do it in goDown)
    // Skip $-edges.
    while (suffixLength(saAt(saRange.i1, index), index) <= value(it).repLen)
    {
        // Interval contains only $-edges.
        if (++saRange.i1 >= saRange.i2)
            return false;
    }

    // Get first and last characters in interval.
    // CHANGED(sascha): get char via the suffix instead of textAt
    TAlphabet cLeft  = value(suffix(index, saAt(saRange.i1,     index)), value(it).repLen);
    TAlphabet cRight = value(suffix(index, saAt(saRange.i2 - 1, index)), value(it).repLen);

    SEQAN_ASSERT_NEQ(ordValue(cLeft), ordValue(value(it).lastChar));

    // Update left range.
    value(it).range.i1 = saRange.i1;

    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);
        
        // CHANGED(sascha): handle index instead of text
        TSAIterator upperBound = _upperBoundSA(index, node, cLeft, value(it).repLen);
        
        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }

    // Update repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;

    return true;
}
    
    
    
    
    
    
    */
    
    
    
    
    
    
    
    
    

}

#endif // SEQAN_HEADER_INDEX_GENREALIZED_SA_H

