// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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
// This file contains type overloadings.
// ==========================================================================

#ifndef APP_YARA_MISC_TYPES_H_
#define APP_YARA_MISC_TYPES_H_

using namespace seqan;

// ============================================================================
// Yara Limits
// ============================================================================

// ----------------------------------------------------------------------------
// Class YaraBits
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct YaraBits
{
    static const unsigned CONTIG_ID   = 8;
    static const unsigned CONTIG_SIZE = 30;
    static const unsigned READ_ID     = 21;
    static const unsigned READ_SIZE   = 14;
    static const unsigned ERRORS      = 6;
};

// ----------------------------------------------------------------------------
// Class YaraLimits
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct YaraLimits
{
    static const unsigned CONTIG_ID   = Power<2, YaraBits<TSpec>::CONTIG_ID>::VALUE - 1;
    static const unsigned CONTIG_SIZE = Power<2, YaraBits<TSpec>::CONTIG_SIZE>::VALUE - 1;
    static const unsigned READ_ID     = Power<2, YaraBits<TSpec>::READ_ID>::VALUE - 1;
    static const unsigned READ_SIZE   = Power<2, YaraBits<TSpec>::READ_SIZE>::VALUE - 1;
    static const unsigned ERRORS      = Power<2, YaraBits<TSpec>::ERRORS>::VALUE - 1;
};

// ============================================================================
// String Types
// ============================================================================

// ----------------------------------------------------------------------------
// String Spec
// ----------------------------------------------------------------------------

//#ifndef YARA_INDEXER
//typedef MMap<>  YaraStringSpec;
//#else
typedef Alloc<> YaraStringSpec;
//#endif

// ============================================================================
// Index Types
// ============================================================================

// ----------------------------------------------------------------------------
// Index Text Type
// ----------------------------------------------------------------------------

typedef StringSet<String<Dna5, Packed<YaraStringSpec> >, Owner<ConcatDirect<> > >   YaraContigs;
typedef StringSet<String<Dna>, Owner<ConcatDirect<> > >                             YaraContigsFM;

// ----------------------------------------------------------------------------
// FM Index Fibres
// ----------------------------------------------------------------------------

struct YaraFMIndexConfig
{
    typedef TwoLevels<void>    TValuesSpec;
    typedef Naive<void>        TSentinelsSpec;

    static const unsigned SAMPLING = 10;
};

typedef FMIndex<void, YaraFMIndexConfig>        YaraIndexSpec;
typedef Index<YaraContigsFM, YaraIndexSpec>     YaraIndex;

// ----------------------------------------------------------------------------
// FM Index Size
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct Size<YaraIndex>
{
    typedef __uint32 Type;
};
}

// ----------------------------------------------------------------------------
// Default Index Fibre Specs
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct DefaultIndexStringSpec<YaraContigsFM>
{
    typedef YaraStringSpec Type;
};

template <>
struct DefaultIndexStringSpec<CompressedSA<YaraContigsFM, void, YaraFMIndexConfig> >
{
    typedef YaraStringSpec Type;
};
}

// ----------------------------------------------------------------------------
// Suffix Array Value Type
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct StringSetPosition<YaraContigs>
{
    typedef Pair<__uint8, __uint32, Pack> Type;
};

template <>
struct StringSetPosition<YaraContigsFM>
{
    typedef Pair<__uint8, __uint32, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// FibreLF Size
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec, typename TConfig>
struct Size<LF<YaraContigsFM, TSpec, TConfig> >
{
    typedef __uint32 Type;
};
}

// ----------------------------------------------------------------------------
// Rank Dictionary Size
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct Size<RankDictionary<Dna, TwoLevels<TSpec> > >
{
    typedef __uint32 Type;
};

template <typename TSpec>
struct Size<RankDictionary<bool, TwoLevels<TSpec> > >
{
    typedef __uint32 Type;
};

template <typename TSpec>
struct Size<RankDictionary<bool, Naive<TSpec> > >
{
    typedef __uint32 Type;
};
}

// ----------------------------------------------------------------------------
// Rank Dictionary Fibre Specs
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct RankDictionaryFibreSpec<RankDictionary<Dna, TwoLevels<TSpec> > >
{
    typedef YaraStringSpec Type;
};

template <typename TSpec>
struct RankDictionaryFibreSpec<RankDictionary<bool, TwoLevels<TSpec> > >
{
    typedef YaraStringSpec Type;
};

template <typename TSpec>
struct RankDictionaryFibreSpec<RankDictionary<bool, Naive<TSpec> > >
{
    typedef YaraStringSpec Type;
};
}

// ----------------------------------------------------------------------------
// CSA Size
// ----------------------------------------------------------------------------

//namespace seqan {
//// TODO(esiragusa): Overload Size<CSA> instead of Size<SparseString>
//template <typename TValueString>
//struct Size<SparseString<TValueString, void> >
//{
//    typedef __uint32    Type;
//};
//}

#endif  // #ifndef APP_YARA_MISC_TYPES_H_
