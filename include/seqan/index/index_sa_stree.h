// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INDEX_SA_STREE_H_
#define SEQAN_INDEX_SA_STREE_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TString, typename TSpec>
class SearchTreeIterator;

struct SortedList;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = void>
struct IndexSa {};

/**
.Spec.IndexSa:
..summary:An index based on a suffix array.
..cat:Index
..general:Class.Index
..signature:Index<TText, IndexSa<> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..include:seqan/index.h
*/

/*!
 * @class IndexSa
 * @extends Index
 * @headerfile <seqan/index.h>
 * @brief An index based on a suffix array.
 * @signature template <typename TText, typename TSpec>
 *            class Index<TText, IndexSa<TSpec> >;
 *
 * @tparam TText The type of the underlying @link TextConcept text @endlink.
 * @tparam TSpec A tag for specialization purposes.
 */

template <typename TText, typename TSpec>
class Index<TText, IndexSa<TSpec> >
{
public:
    typename Member<Index, FibreText>::Type         text;
    typename Fibre<Index, FibreSA>::Type            sa;
    typename Cargo<Index>::Type                     cargo;

    Index() {}

    Index(Index & other) :
        text(other.text),
        sa(other.sa)
    {}

    Index(Index const & other) :
        text(other.text),
        sa(other.sa)
    {}

    template <typename TText_>
    Index(TText_ & _text) :
        text(_text)
    {}

    template <typename TText_>
    Index(TText_ const & _text) :
        text(_text)
    {}
};

template <typename TSize, typename TAlphabet>
struct VertexSA : public VertexEsa<TSize>
{
    typedef VertexEsa<TSize>                        TBase;

    TSize       repLen;
    TAlphabet   lastChar;

    SEQAN_HOST_DEVICE
    VertexSA() :
        TBase(),
        repLen(0),
        lastChar(0)
    {}

    SEQAN_HOST_DEVICE
    VertexSA(MinimalCtor) :
        TBase(MinimalCtor()),
        repLen(0),
        lastChar(0)
    {}

    SEQAN_HOST_DEVICE
    VertexSA(VertexSA const & other) :
        TBase(other),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}
};

template <typename TSize, typename TAlphabet>
struct HistoryStackSA_
{
    Pair<TSize> range;
    TSize       repLen;
    TAlphabet   lastChar;
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TSpec>
struct VertexDescriptor<Index<TText, IndexSa<TSpec> > >
{
    typedef Index<TText, IndexSa<TSpec> >           TIndex;
    typedef typename Size<TIndex>::Type             TSize;
    typedef typename Value<TIndex>::Type            TAlphabet;

    typedef VertexSA<TSize, TAlphabet>              Type;
};

template <typename TText, typename TIndexSpec, typename TSpec>
struct HistoryStackEntry_<Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > >
{
private:
    typedef Index<TText, IndexSa<TIndexSpec> >      TIndex;
    typedef typename Size<TIndex>::Type             TSize;
    typedef typename Value<TIndex>::Type            TAlphabet;

public:
    typedef HistoryStackSA_<TSize, TAlphabet>       Type;
};

template <typename TText, typename TIndexSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > >
{
    typedef typename Value<Index<TText, IndexSa<TIndexSpec> > >::Type Type;
};

template < typename TText, typename TIndexSpec >
struct DefaultFinder< Index<TText, IndexSa<TIndexSpec> > >
{
    typedef FinderMlr Type;	// standard suffix array finder is mlr-heuristic
};


// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TIndexSpec>
SEQAN_HOST_DEVICE inline void _indexRequireTopDownIteration(Index<TText, IndexSa<TIndexSpec> > & index)
{
    indexRequire(index, FibreSA());
}

template <typename TText>
void _indexRequireTopDownIteration(Index<TText, IndexSa<InfixSegment> > &)
{
    // The SA fibre must be provided by calling setHost explicitely.
}

template <typename TText, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline typename SAValue<Index<TText, IndexSa<TIndexSpec> > >::Type
_lastOccurrence(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    if (_isSizeInval(value(it).range.i2))
        return back(indexSA(container(it)));
    else
        return saAt(value(it).range.i2 - 1, container(it));
}


// is this the first child of the parent node?
template <typename TIndex, typename TSpec>
inline bool
isFirstChild(Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > const & it)
{
    if (!empty(it.history))
    {
        typename Size<TIndex>::Type parentLeft = back(it.history).range.i1;
        if (value(it).range.i1 == parentLeft)
            return true;

        TIndex const & index = container(it);
        return suffixLength(saAt(value(it).range.i1 - 1, index), index) + 1 == value(it).repLen;
    }
    else
        return value(it).range.i1 == 0;

}

// is this the last child of the parent node?
template <typename TIndex, typename TSpec>
inline bool
isLastChild(Iter<TIndex, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).range.i2 == value(it).parentRight;
}

template <typename TIndex, typename TSpec>
inline bool
isLastChild(Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > const & it)
{
    if (!empty(it.history))
        return value(it).range.i2 == back(it.history).range.i2;
    else
        return _isSizeInval(value(it).range.i2);

}

// go down and skip non-branching nodes
// this function emulates suffix tree traversal
template <typename TIndex, typename TSpec>
inline bool
goDownSkipSingletons(Iter<TIndex, VSTree< TopDown<TSpec> > > &it)
{
    do {
        if (!goDown(it))
            return false;
    } while (isLastChild(it));
    return true;
}

// is this a leaf? (hide empty $-edges)
template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_isLeaf(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const & it,
        VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >                      TIndex;
    typedef typename SAValue<TIndex>::Type                          TOcc;

//    if (_isLeaf(value(it))) return true;

    TIndex const & index = container(it);

    typename Size<TIndex>::Type lcp = repLength(it);

    // if the last suffix in the interval is larger than the lcp,
    // not all outgoing edges are empty (uses lex. sorting)
    TOcc oc = _lastOccurrence(it);

    return suffixLength(oc, index) == lcp;
//  return getSeqOffset(oc, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(oc, stringSetLimits(index)), index);
}

template <typename TIndex, typename TSize, typename TAlphabet>
SEQAN_HOST_DEVICE inline typename Size<TIndex>::Type
repLength(TIndex const &, VertexSA<TSize, TAlphabet> const & vDesc)
{
    return vDesc.repLen;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Value<Index<TText, IndexSa<TIndexSpec> > >::Type
parentEdgeFirstChar(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

template <typename TText, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline void goRoot(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > & it)
{
    _historyClear(it);
    clear(it);
    if (!empty(indexSA(container(it))))
        _setSizeInval(value(it).range.i2);
    value(it).repLen = 0;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    // TODO(esiragusa): goDown including empty $-edges
    return _goDown(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;

#ifdef SEQAN_DEBUG
    std::cout << "goDown" << std::endl;
#endif

    if (_isLeaf(it, HideEmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);

    // TODO(esiragusa): check nodeHullPredicate

#ifdef SEQAN_DEBUG
    std::cout << "parent: " << value(it).range << std::endl;
#endif

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

    typedef typename SuffixFunctor<TIndex, typename Value<TSA>::Type>::result_type TSuffix;
    SuffixFunctor<TIndex, typename Value<TSA>::Type> dereferer(index);

    TAlphabet cLeft, cRight;

    // switch for gapped suffixes; note that both should result in the same characters,
    // the former one might just have less overhead.
    if (IsSameType<TIndexSpec, void>::VALUE)
    {
        cLeft = textAt(posAdd(saAt(saRange.i1, index), value(it).repLen), index);
        cRight= textAt(posAdd(saAt(saRange.i2 - 1, index), value(it).repLen), index);
    }
    else
    {
        // generate suffixes and the take the char at position 'repLen'
        cLeft  = value(dereferer(saAt(saRange.i1, index)),    value(it).repLen);
        cRight = value(dereferer(saAt(saRange.i2 -1, index)), value(it).repLen);
    }


#ifdef SEQAN_DEBUG
    std::cout << "char: " << Pair<TAlphabet>(cLeft, cRight) << std::endl;
#endif

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

        SuffixFunctor<TIndex, typename Value<TSA>::Type> dereferer(index);
        TSAIterator upperBound = _upperBoundSA(dereferer, node, cLeft, value(it).repLen);

        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }
    
    // Update child repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;

#ifdef SEQAN_DEBUG
    std::cout << "child: " <<  value(it).range << std::endl;
#endif

    return true;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges>
inline bool _goRight(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                     VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;

#ifdef SEQAN_DEBUG
    std::cout << "goRight" << std::endl;
#endif

    if (isRoot(it))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);

#ifdef SEQAN_DEBUG
    std::cout << "current: " << value(it).range << std::endl;
#endif

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

    typedef typename SuffixFunctor<TIndex, typename Value<TSA>::Type>::result_type TSuffix;
    SuffixFunctor<TIndex, typename Value<TSA>::Type> dereferer(index);

    TAlphabet cLeft, cRight;

    // switch for gapped suffixes; note that both should result in the same characters,
    // the former one might just have less overhead.
    if (IsSameType<TIndexSpec, void>::VALUE)
    {
        cLeft = textAt(posAdd(saAt(saRange.i1, index), value(it).repLen), index);
        cRight= textAt(posAdd(saAt(saRange.i2 - 1, index), value(it).repLen), index);
    }
    else
    {
        // generate suffixes and the take the char at position 'repLen'
        cLeft  = value(dereferer(saAt(saRange.i1, index)),    value(it).repLen);
        cRight = value(dereferer(saAt(saRange.i2 -1, index)), value(it).repLen);
    }



    SEQAN_ASSERT_NEQ(ordValue(cLeft), ordValue(value(it).lastChar));

#ifdef SEQAN_DEBUG
    std::cout << "char: " << Pair<TAlphabet>(cLeft, cRight) << std::endl;
#endif

    // Update left range.
    value(it).range.i1 = saRange.i1;

    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);

        TSAIterator upperBound = _upperBoundSA(dereferer, node, cLeft, value(it).repLen);

        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }

    // Update repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;

#ifdef SEQAN_DEBUG
    std::cout << "sibling: " <<  value(it).range << std::endl;
#endif

    return true;
}

// NOTE(meiers): Shortcut for going down edges character-wise. If the borders are labeled already
//               with the specific char no binary search has to be performed.
template <typename TText, typename TIndexSpec, typename TSpec, typename TValue>
inline bool _isEdge(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, TValue c)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;

    if(isRoot(it)) return false;
    if(value(it).range.i2 == value(it).range.i1) return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TSASize repLen = value(it).repLen;
    TSAIterator l = begin(sa, Standard()) + value(it).range.i1;
    TSAIterator r = begin(sa, Standard()) + value(it).range.i2 - 1;
    
    typedef SuffixFunctor<TIndex, typename Value<TSA>::Type>    TSufFunctor;
    typedef typename TSufFunctor::result_type                   TSuffix;
    typedef typename Iterator<TSuffix const>::Type              TSufIter;

    TSufFunctor dereferer(index);
    TSufIter lit = begin(dereferer(*l)) + repLen, len = end(dereferer(*l));
    TSufIter rit = begin(dereferer(*r)) + repLen, ren = end(dereferer(*r));

    if (lit < len && rit < ren && *lit == c && *rit == c)
        return true;
 
    return false;
}


template <typename TText, typename TIndexSpec, typename TSpec, typename TValue>
inline bool _goDownChar(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, TValue c)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;

    // Save vertex descriptor.
    _historyPush(it);

    TIndex const & index   = container(it);
    TSA const &    sa      = indexSA(index);
    TSAIterator    saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize        saLen   = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;

    TSearchTreeIterator node(saBegin, saLen);

    typedef SuffixFunctor<TIndex, typename Value<TSA>::Type> TSuf;
    TSuf  dereferer(index);

    // Check whether this is an edge, not a node of the suffix tree
    // NOTE(meiers): disable this branch because it is slower
    if ( false   &&  _isEdge(it,c))
    {
        value(it).lastChar = c;
        value(it).repLen++;
        return true;
    }
    else
    { 
        Pair<TSAIterator> range = _equalRangeSA(dereferer, node, c, value(it).repLen);

        if (!(range.i1 < range.i2))
            return false;

        // Save vertex descriptor.
        _historyPush(it);
    
        // Update range, lastChar and repLen.
        value(it).range.i1 = range.i1 - begin(sa, Standard());
        value(it).range.i2 = range.i2 - begin(sa, Standard());
        value(it).lastChar = c;
        value(it).repLen++;
    
        return true;
    }
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                          TString const & pattern, TSize & lcp)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TSA const>::Type                  TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;


    if (_isLeaf(it, HideEmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);

    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);
    SuffixFunctor<TIndex, typename Value<TSA>::Type> dereferer(index);
    Pair<TSAIterator> range = _equalRangeSA(dereferer, node, pattern, value(it).repLen);

    if (range.i1 >= range.i2)
        return false;
    
    // Save vertex descriptor.
    _historyPush(it);

    value(it).range.i1 = range.i1 - begin(sa, Standard());
    value(it).range.i2 = range.i2 - begin(sa, Standard());
    value(it).repLen += length(pattern);
    value(it).lastChar = back(pattern);

    lcp = length(pattern);

    return true;
}

// return vertex descriptor of parent's node
template <typename TText, typename TIndexSpec, typename TSpec>
inline typename VertexDescriptor<Index<TText, IndexSa<TIndexSpec> > >::Type
nodeUp(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > const & it)
{
    typedef Index<TText, IndexSa<TIndexSpec> > TIndex;
    
    if (!empty(it.history))
    {
        typename VertexDescriptor<TIndex>::Type desc;
        desc.range = back(it.history).range;
        desc.repLen = back(it.history).repLen;
        desc.lastChar = back(it.history).lastChar;
        if (length(it.history) >= 2)
            desc.parentRight = backPrev(it.history).range.i2;
        else
            desc.parentRight = value(it).parentRight;
        return desc;
    } else
        return value(it);
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool _goUp(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
#ifdef SEQAN_DEBUG
    std::cout << "goUp" << std::endl;
#endif

    if (!empty(it.history))
    {
        value(it).range = back(it.history).range;
        value(it).repLen = back(it.history).repLen;
        value(it).lastChar = back(it.history).lastChar;
        eraseBack(it.history);
        if (!empty(it.history))
            value(it).parentRight = back(it.history).range.i2;
        return true;
    }
    return false;
}

// ============================================================================

template <typename TText, typename TIndexSpec, typename TSpec>
inline void _historyPush(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    it._parentDesc = value(it);
    value(it).parentRight = value(it).range.i2;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline void _historyPush(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    typedef Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;
    typename HistoryStackEntry_<TIter>::Type h;
    h.range = value(it).range;
    h.repLen = value(it).repLen;
    h.lastChar = value(it).lastChar;
    appendValue(it.history, h);

    value(it).parentRight = value(it).range.i2;
}

// ============================================================================

template <typename TText, typename TSpec>
inline void clear(Index<TText, IndexSa<TSpec> > & index)
{
    clear(getFibre(index, FibreSA()));
}

template <typename TObject, typename TSpec>
inline bool open(Index<TObject, IndexSa<TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".txt");
    if ((!open(getFibre(index, FibreText()), toCString(name), openMode)) &&
        (!open(getFibre(index, FibreText()), fileName, openMode))) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    return true;
}

template <typename TObject, typename TSpec>
inline bool open(Index<TObject, IndexSa<TSpec> > & index, const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TObject, IndexSa<TSpec> > >::VALUE);
}

template <typename TObject, typename TSpec>
inline bool save(Index<TObject, IndexSa<TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;
    
    name = fileName;    append(name, ".txt");
    if ((!save(getFibre(index, FibreText()), toCString(name), openMode)) &&
        (!save(getFibre(index, FibreText()), fileName, openMode))) return false;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    return true;
}

template <typename TObject, typename TSpec>
inline bool save(Index<TObject, IndexSa<TSpec> > & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TObject, IndexSa<TSpec> > >::VALUE);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INDEX_SA_STREE_H_
