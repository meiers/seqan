// ==========================================================================
//                                gappedIndex
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

#ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_FIND_H_
#define CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_FIND_H_

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

using namespace seqan;


/*!
 * Test the Finder Interface of gapped Suffixes
 * Idea: 
 * 1. Given a text and a set of queries, search the queries in the 
 *    text using an existing (correct) method. 
 * 2. Modify the queries such that the don't care positions are removed.
 * 3. Search these modified queries in the gapped Index. The results we 
 *    obtain here should contain all the hits from the first search, but
 *    might contain more as less characters have to match. Further,
 *    the order might differ.
 */
template <typename TText, typename TQuerySet, typename TCyclicShape>
void testFinder(TText const & text, TQuerySet const & querySet, TCyclicShape const & shape, bool verbose)
{
    typedef typename Value<TQuerySet const>::Type               TQueryString;
    typedef typename Iterator<TQuerySet const, Standard>::Type  TQuerySetIter;
    typedef ModifiedString<TQueryString, ModCyclicShape<TCyclicShape> > TMod;

    TQuerySet modSet;
    for (TQuerySetIter quIt = begin(querySet); quIt != end(querySet); ++quIt)
    {
        TMod m(*quIt, shape);
        appendValue(modSet, m);
    }

    // Find using an existing method
    Index<TText const> correctIndex(text);
    Finder<Index<TText const> > correctFinder(correctIndex);

    typedef Index<TText const, IndexSa<Gapped<ModCyclicShape<TCyclicShape> > > > TIndex;
    TIndex index(text, shape);
    indexCreate(index, FibreSA(), Dislex<>());
    Finder<TIndex> finder(index);

    // for all queries
    TQuerySetIter qu = begin(querySet);
    TQuerySetIter modQu = begin(modSet);
    for(; qu != end(querySet); ++qu, ++modQu)
    {
        typedef typename Position<Finder<Index<TText const> > >::Type TPos;
        String<TPos> res, modRes;

        clear(correctFinder);
        goBegin(correctFinder);
        while (find(correctFinder, *qu))
            appendValue(res, beginPosition(correctFinder));

        clear(finder);
        goBegin(finder);

        while (find(finder, *modQu))
            appendValue(modRes, beginPosition(finder));

        // compare results
        if (res != modRes)
        {
            std::sort (begin(res), end(res));
            std::sort (begin(modRes), end(modRes));

            typedef typename Iterator<String<TPos> >::Type TIt;
            TIt it = begin(res), modIt = begin(modRes);
            while (it != end(res))
            {
                if (modIt == end(modRes)) break;

                // skip mod results that are not in the correct results
                while (modIt != end(modRes) && *modIt != *it)
                    ++modIt;

                ++it;
            }

            if (verbose) std::cout << "(" << length(res) << ":" << length(modRes) << "),";
            SEQAN_ASSERT_EQ_MSG(it, end(res),"Not all results found by the correct index have been found by the gapped index");
        } else {
            if (verbose) std::cout << "(=),";
        }
    }
    if (verbose) std::cout << std::endl;
}


template <typename TText, typename TShape>
void call_testFinder(TText &text, TShape const & shape, bool verbose = false)
{
    // get randomText
    typedef typename Value<TText>::Type TAlph;
    String<TAlph> s;
    generateText(s, 10000);
    clear(text);
    assign(text, s);

    // get randomQueries
    StringSet<String<TAlph> > queries;
    generatePattern(queries, s, 100);

    // call test
    testFinder(text, queries, shape, verbose);
}

template <typename TText, typename TSpec, typename TShape>
void call_testFinder(StringSet<TText, TSpec> &text, TShape const & shape, bool verbose = false)
{
    // get randomText
    typedef typename Value<TText>::Type TAlph;
    StringSet<String<TAlph> > s;
    generateText(s, 30, 500);
    clear(text);
    assign(text, s);

    // get randomQueries
    StringSet<String<TAlph> > queries;
    generatePattern(queries, s, 80);

    // call test
    testFinder(text, queries, shape, verbose);
}


struct _ShapeDefs
{
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<> >, 1> >       S_10;
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,2> >, 1> >    S_11010;
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,1,1> >, 2> >  S_111100;
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<4> >, 0> >      S_10001;
    CyclicShape<FixedShape<1,GappedShape<HardwiredShape<> >, 0> >       S_01;
    CyclicShape<FixedShape<2,GappedShape<HardwiredShape<1> >, 0> >      S_0011;

    CyclicShape<GenericShape> s_10;
    CyclicShape<GenericShape> s_11010;
    CyclicShape<GenericShape> s_111100;
    CyclicShape<GenericShape> s_10001;
    CyclicShape<GenericShape> s_01;
    CyclicShape<GenericShape> s_0011;

    _ShapeDefs()
    {
        stringToCyclicShape(s_10,       "10");
        stringToCyclicShape(s_11010,    "11010");
        stringToCyclicShape(s_111100,   "111100");
        stringToCyclicShape(s_10001,    "10001");
        stringToCyclicShape(s_01,       "01");
        stringToCyclicShape(s_0011,     "0011");
    }
};



// Test Find for Shape 10

SEQAN_DEFINE_TEST(test_gappedIndex_find_10_DnaString)
{
    _ShapeDefs SD;
    DnaString str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_Dna5String)
{
    _ShapeDefs SD;
    Dna5String str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_Peptide)
{
    _ShapeDefs SD;
    Peptide str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_CharString)
{
    _ShapeDefs SD;
    CharString str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_DnaString_Set)
{
    _ShapeDefs SD;
    StringSet<DnaString> str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_Dna5String_Set)
{
    _ShapeDefs SD;
    StringSet<Dna5String> str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_Peptide_Set)
{
    _ShapeDefs SD;
    StringSet<Peptide> str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_CharString_Set)
{
    _ShapeDefs SD;
    StringSet<CharString> str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}


// special
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_Finite)
{
    _ShapeDefs SD;
    String<SimpleType<char, Finite<256> > > str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10_Finite_Set)
{
    _ShapeDefs SD;
    StringSet<String<SimpleType<unsigned, Finite<256> > > > str;
    call_testFinder(str, SD.S_10);
    call_testFinder(str, SD.s_10);
}



// Test Find for Shape 11010

SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_DnaString)
{
    _ShapeDefs SD;
    DnaString str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_Dna5String)
{
    _ShapeDefs SD;
    Dna5String str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_Peptide)
{
    _ShapeDefs SD;
    Peptide str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_CharString)
{
    _ShapeDefs SD;
    CharString str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_DnaString_Set)
{
    _ShapeDefs SD;
    StringSet<DnaString> str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_Dna5String_Set)
{
    _ShapeDefs SD;
    StringSet<Dna5String> str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_Peptide_Set)
{
    _ShapeDefs SD;
    StringSet<Peptide> str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_11010_CharString_Set)
{
    _ShapeDefs SD;
    StringSet<CharString> str;
    call_testFinder(str, SD.S_11010);
    call_testFinder(str, SD.s_11010);
}



// Test Find for Shape 111100

SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_DnaString)
{
    _ShapeDefs SD;
    DnaString str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_Dna5String)
{
    _ShapeDefs SD;
    Dna5String str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_Peptide)
{
    _ShapeDefs SD;
    Peptide str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_CharString)
{
    _ShapeDefs SD;
    CharString str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_DnaString_Set)
{
    _ShapeDefs SD;
    StringSet<DnaString> str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_Dna5String_Set)
{
    _ShapeDefs SD;
    StringSet<Dna5String> str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_Peptide_Set)
{
    _ShapeDefs SD;
    StringSet<Peptide> str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_111100_CharString_Set)
{
    _ShapeDefs SD;
    StringSet<CharString> str;
    call_testFinder(str, SD.S_111100);
    call_testFinder(str, SD.s_111100);
}



// Test Find for Shape 10001

SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_DnaString)
{
    _ShapeDefs SD;
    DnaString str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_Dna5String)
{
    _ShapeDefs SD;
    Dna5String str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_Peptide)
{
    _ShapeDefs SD;
    Peptide str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_CharString)
{
    _ShapeDefs SD;
    CharString str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_DnaString_Set)
{
    _ShapeDefs SD;
    StringSet<DnaString> str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_Dna5String_Set)
{
    _ShapeDefs SD;
    StringSet<Dna5String> str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_Peptide_Set)
{
    _ShapeDefs SD;
    StringSet<Peptide> str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_10001_CharString_Set)
{
    _ShapeDefs SD;
    StringSet<CharString> str;
    call_testFinder(str, SD.S_10001);
    call_testFinder(str, SD.s_10001);
}


// Test Find for Shape 01

SEQAN_DEFINE_TEST(test_gappedIndex_find_01_DnaString)
{
    _ShapeDefs SD;
    DnaString str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_Dna5String)
{
    _ShapeDefs SD;
    Dna5String str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_Peptide)
{
    _ShapeDefs SD;
    Peptide str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_CharString)
{
    _ShapeDefs SD;
    CharString str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_DnaString_Set)
{
    _ShapeDefs SD;
    StringSet<DnaString> str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_Dna5String_Set)
{
    _ShapeDefs SD;
    StringSet<Dna5String> str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_Peptide_Set)
{
    _ShapeDefs SD;
    StringSet<Peptide> str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_01_CharString_Set)
{
    _ShapeDefs SD;
    StringSet<CharString> str;
    call_testFinder(str, SD.S_01);
    call_testFinder(str, SD.s_01);
}



// Test Find for Shape 0011

SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_DnaString)
{
    _ShapeDefs SD;
    DnaString str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_Dna5String)
{
    _ShapeDefs SD;
    Dna5String str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_Peptide)
{
    _ShapeDefs SD;
    Peptide str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_CharString)
{
    _ShapeDefs SD;
    CharString str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_DnaString_Set)
{
    _ShapeDefs SD;
    StringSet<DnaString> str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_Dna5String_Set)
{
    _ShapeDefs SD;
    StringSet<Dna5String> str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_Peptide_Set)
{
    _ShapeDefs SD;
    StringSet<Peptide> str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}
SEQAN_DEFINE_TEST(test_gappedIndex_find_0011_CharString_Set)
{
    _ShapeDefs SD;
    StringSet<CharString> str;
    call_testFinder(str, SD.S_0011);
    call_testFinder(str, SD.s_0011);
}


#endif  // CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_FIND_H_
