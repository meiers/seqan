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

#ifndef SEQAN_HEADER_INDEX_GAPPED_SA_H
#define SEQAN_HEADER_INDEX_GAPPED_SA_H

// TODO(meiers): Fibre access to modifierCargo

namespace SEQAN_NAMESPACE_MAIN
{
    
// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct IndexSa;
template <typename T, typename TSAValue> struct SuffixFunctor;
    
// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(meiers): Remove this definition as soon as Ticket #1096 has been solved: http://trac.seqan.de/ticket/1096
template < typename TText, typename TSpec >
struct DefaultFinder< Index<TText, IndexSa<TSpec> > > {
    typedef EsaFindMlr Type;	// standard suffix array finder is mlr-heuristic
};

template <typename TSuffixMod, typename TSpec = void>
struct Gapped {};

// ----------------------------------------------------------------------------
// SuffixFunctor
// ----------------------------------------------------------------------------

template <typename TText, typename TSAValue, typename TModifier, typename TSpec>
struct SuffixFunctor < Index<TText, IndexSa<Gapped<TModifier, TSpec> > >, TSAValue> :
    std::unary_function<TSAValue, typename Suffix<Index<TText, IndexSa<Gapped<TModifier, TSpec> > > const>::Type>
{
    typedef Index<TText, IndexSa<Gapped<TModifier, TSpec> > >   TIndex;
    typedef typename Fibre<TIndex, FibreText>::Type             TText2;
    typedef typename Suffix<TIndex const>::Type                 result_type;
    typedef typename Cargo<result_type>::Type                   TModCargo;

    TText2 const    &text;
    TModCargo const &modCargo;

    SuffixFunctor(TIndex const &index) :
        text(indexText(index)),
        modCargo(index.modifierCargo)
    {}

    result_type
    operator() (TSAValue const &pos) const
    {
        return result_type(suffix(text, pos), modCargo);
    }
};


// Carries a modifierCargo
// of type TModifierCargo a.k.a. Cargo<Suffix<Index> >
template <typename TText, typename TSuffixMod, typename TSpec>
class Index<TText, IndexSa< Gapped<TSuffixMod, TSpec> > > :
    public Index<TText, IndexSa<TSpec> >
{
public:
    typedef Index<TText, IndexSa<TSpec> >           TBase;
    typedef typename Suffix<Index>::Type            TSuffix;
    typedef typename Cargo<TSuffix>::Type           TModifierCargo;
    
    // derive           text, cargo, sa
    TModifierCargo      modifierCargo;
    
    Index() {}
    
    Index(Index & other) :
        TBase(static_cast<TBase &>(other)),
        modifierCargo(other.modifierCargo)
    {}

    Index(Index const & other) :
        TBase(static_cast<TBase const &>(other)),
        modifierCargo(other.modifierCargo)
    {}

    template <typename TText_>
    Index(TText_ & _text) :
        TBase(_text)
    {}
    
    template <typename TText_>
    Index(TText_ const & _text) :
        TBase(_text)
    {}

    template <typename TText_>
    Index(TText_ & _text, TModifierCargo const & modifierCargo) :
        TBase(_text),
        modifierCargo(modifierCargo)
    {}

    template <typename TText_>
    Index(TText_ const & _text, TModifierCargo const & modifierCargo) :
        TBase(_text),
        modifierCargo(modifierCargo)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// DefaultIndexCreator
// ----------------------------------------------------------------------------

// TODO(meiers): Update to a better method like RadixSort
// TODO(meiers): Specialise this for CyclicShape to be Dislex
template < typename TText, typename TShape, typename TSpec>
struct DefaultIndexCreator<Index<TText, IndexSa<Gapped<ModCyclicShape<TShape>,TSpec> > >, FibreSA>
{
    typedef Index<TText, IndexSa<Gapped<ModCyclicShape<TShape>,TSpec> > >   TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                           TSA;

    typedef typename AllowsFastRandomAccess<TSA>::Type                      TRandomSA;
    typedef typename AllowsFastRandomAccess<TText>::Type                    TRandomText;
    typedef typename If<And<TRandomText,TRandomSA>,
                        SAQSort,
                        DislexExternal<TShape> >::Type                      Type;
};

// ----------------------------------------------------------------------------
// Metafunction Suffix                                                  general
// ----------------------------------------------------------------------------

// general modified suffix
template <typename TText, typename TSuffixMod, typename TSpec>
struct Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >
{
    typedef ModifiedString<typename Suffix<TText>::Type, TSuffixMod>    Type;
};

// general modified suffix; const variant
template <typename TText, typename TSuffixMod, typename TSpec>
struct Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>
{
    typedef ModifiedString<typename Suffix<TText const>::Type, TSuffixMod >     Type;
};

// ----------------------------------------------------------------------------
// Metafunction Infix                                                   general
// ----------------------------------------------------------------------------

// general modified suffix
template <typename TText, typename TSuffixMod, typename TSpec>
struct Infix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >
{
    typedef typename Prefix<typename Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type>::Type Type;
};

// general modified suffix; const variant
template <typename TText, typename TSuffixMod, typename TSpec>
struct Infix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>
{
    typedef typename Prefix<typename Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function suffix()
// ----------------------------------------------------------------------------

template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin>
inline typename Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type
suffix(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > & t, TPosBegin pos_begin)
{
    typedef typename Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type TModifiedSuffix;
    return TModifiedSuffix(suffix(indexText(t), pos_begin), t.modifierCargo);
}

// const variant
template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin>
inline typename Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type
suffix(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const &t, TPosBegin pos_begin)
{
    typedef typename Suffix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type TModifiedSuffix;
    return TModifiedSuffix(suffix(indexText(t), pos_begin), t.modifierCargo);
}

// ----------------------------------------------------------------------------
// Function suffixLength()
// ----------------------------------------------------------------------------

// TODO(meiers): There is a problem with ambigous calls to suffixLength, see definitions in index_base, too.
/*
template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin>
SEQAN_HOST_DEVICE inline typename Size<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type
suffixLength(TPosBegin pos, Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > &index)
{
    return length(suffix(index, pos) );
}

template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin>
SEQAN_HOST_DEVICE inline typename Size<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type
suffixLength(TPosBegin pos, Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const &index)
{
    return length(suffix(index, pos) );
}
*/
// ----------------------------------------------------------------------------
// Function infixWithLength()
// ----------------------------------------------------------------------------

template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type
infixWithLength(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > &index, TPosBegin pos_begin, TSize length)
{
    return prefix(suffix(index, pos_begin), length);
}

template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type
infixWithLength(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const &index, TPosBegin pos_begin, TSize length)
{
    return prefix(suffix(index, pos_begin), length);
}

// ----------------------------------------------------------------------------
// Function infix()                                                     general
// ----------------------------------------------------------------------------

template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type
infix(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > &t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(t, pos_begin), pos_end - pos_begin);
}

// const variant
template <typename TText, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type
infix(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const &t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(t, pos_begin), pos_end - pos_begin);
}

// ----------------------------------------------------------------------------
// Function infix()                                      for Index on StringSet
// ----------------------------------------------------------------------------

template <typename TText, typename TTextSpec, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Index<StringSet<TText, TTextSpec>, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type
infix(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > &t, TPosBegin pos_begin, TPosEnd pos_end)
{

    return prefix(suffix(t, pos_begin), pos_end.i2 - pos_begin.i2);
}

// const variant
template <typename TText, typename TTextSpec, typename TSuffixMod, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Index<StringSet<TText, TTextSpec>, IndexSa<Gapped<TSuffixMod, TSpec> > > const>::Type
infix(Index<StringSet<TText, TTextSpec>, IndexSa<Gapped<TSuffixMod, TSpec> > > const &t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return prefix(suffix(t, pos_begin), pos_end.i2 - pos_begin.i2);
}

// ----------------------------------------------------------------------------
// Function indexCreate (FibreSA)                                       general
// ----------------------------------------------------------------------------

// TODO(meiers): more parameters, see createSuffixArray(... K, ...)
template <typename TText, typename TSuffixMod, typename TSpec, typename TSpecAlg>
inline bool indexCreate(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > &t, FibreSA, TSpecAlg const alg)
{
    resize(indexSA(t), length(indexRawText(t)), Exact());
    createGappedSuffixArray(indexSA(t), indexText(t), t.modifierCargo, TSuffixMod(), alg);
    return true;
}

/*
 // ----------------------------------------------------------------------------
 // Function suffixModifier
 // ----------------------------------------------------------------------------

 template <typename TText, typename TSuffixMod, typename TSpec>
 inline typename Cargo<Suffix<>::Type>::Type &
 suffixModifier(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > & t)
 {
 return t.suffMod;
 }

 // TODO(meiers): Replace & by seqan Reference type

 template <typename TText, typename TSuffixMod, typename TSpec>
 inline typename Cargo<Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > >::Type const &
 suffixModifier(Index<TText, IndexSa<Gapped<TSuffixMod, TSpec> > > const & t)
 {
 return t.suffMod;
 }
 */
}

#endif // SEQAN_HEADER_INDEX_GAPPED_SA_H

