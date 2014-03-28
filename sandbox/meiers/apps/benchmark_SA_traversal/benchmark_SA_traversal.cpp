// ==========================================================================
//                           benchmark_SA_traversal
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

#include <iostream>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class SATraversalOptions
// --------------------------------------------------------------------------

struct SATraversalOptions
{
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<> >, 1> > TShape_1;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1> >, 1> > TShape_2;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,7> >, 0> > TShape_3;
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,3,1,2> >, 1> > TShape_4;
    typedef CyclicShape<FixedShape<0,ShapePatternHunter, 0> > TShape_5;
    typedef CyclicShape<FixedShape<0,ShapePatternHunter, 2> > TShape_6;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<2,1,1,3,2,2,1,1,2,1,2,2,1,1,2,2,2,1> >, 0> > TShape_7;

    CharString file;
    void print()
    {
        std::cout << "RunSACA Options" << std::endl;
        std::cout << "    file:      " << file << std::endl;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(SATraversalOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("Benchmark SA traversal");
    setShortDescription(parser, "Load a file, build the suffix array and then traverse it randomly or search for random patterns");
    setVersion(parser, "0.1");
    setDate(parser, "April 2014");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIFASTA FILE\\fP\"");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    // Extract option values.
    getArgumentValue(options.file, parser, 0);

    return ArgumentParser::PARSE_OK;
}


// --------------------------------------------------------------------------
// Function benchmarkSearch funcitons
// --------------------------------------------------------------------------

template <typename TIndex, typename TPatternSet, typename TSize>
void benchmarkSearchUngapped(TIndex & index, TPatternSet const & set, TSize len)
{
    Finder<TIndex, EsaFindMlr> finder(index);

    typedef typename Iterator<TPatternSet const, Standard>::Type TSetIter;

    TSetIter setIt = begin(set, Standard());
    TSetIter setEnd = end(set, Standard());
    for(; setIt != setEnd; ++setIt)
    {
        clear(finder);
        
        find(finder, prefix(*setIt, len));
    }
}

template <typename TIndex, typename TPatternSet, typename TSize, typename TShape>
void benchmarkSearch(TIndex & index, TPatternSet const & set, TSize len, TShape const & shape)
{
    Finder<TIndex, EsaFindMlr> finder(index);

    typedef typename Iterator<TPatternSet const, Standard>::Type TSetIter;
    typedef ModifiedString<typename Value<TPatternSet const>::Type, ModCyclicShape<TShape> > TModStr;

    TSetIter setIt = begin(set, Standard());
    TSetIter setEnd = end(set, Standard());
    for(; setIt != setEnd; ++setIt)
    {
        clear(finder);
        find(finder, prefix(TModStr(*setIt, shape), len));
    }
}


template <typename TIndex, typename TPatternSet, typename TSize>
void benchmarkSearchUngapped_stree(TIndex & index, TPatternSet const & set, TSize len)
{
    typedef typename Iterator<TPatternSet const, Standard>::Type TSetIter;
    typedef typename Iterator<typename Value<TPatternSet const>::Type, Standard>::Type TPatternIter;
    typedef typename Iterator<TIndex, TopDown<> >::Type TTreeIter;

    TTreeIter treeIt(index);

    TSetIter setIt = begin(set, Standard());
    TSetIter setEnd = end(set, Standard());
    for(; setIt != setEnd; ++setIt)
    {
        TPatternIter patIt = begin(*setIt);
        TPatternIter patEnd = patIt + len;
        for(; patIt !=patEnd; ++patIt)
            goDown(treeIt, *patIt);
    }
}

template <typename TIndex, typename TPatternSet, typename TSize, typename TShape>
void benchmarkSearch_stree(TIndex & index, TPatternSet const & set, TSize len, TShape const & shape)
{
    typedef typename Iterator<TPatternSet const, Standard>::Type    TSetIter;
    typedef ModifiedString<typename Value<TPatternSet const>::Type, ModCyclicShape<TShape> > TModStr;
    typedef typename Iterator<TModStr, Standard>::Type             TPatternIter;
    typedef typename Iterator<TIndex, TopDown<> >::Type             TTreeIter;

    TTreeIter treeIt(index);

    TSetIter setIt = begin(set, Standard());
    TSetIter setEnd = end(set, Standard());
    for(; setIt != setEnd; ++setIt)
    {
        TModStr mod(*setIt);
        TPatternIter patIt = begin(mod);
        TPatternIter patEnd = patIt + len;
        for(; patIt !=patEnd; ++patIt)
            goDown(treeIt, *patIt);
    }
}



// --------------------------------------------------------------------------
// Function callAllTests()
// --------------------------------------------------------------------------


template <typename TText>
void callAllTests(TText const & str, SATraversalOptions const &)
{
    // generate random patterns into a string set

    const unsigned NUMPATTERNS = 100000; // set to 100000
    const unsigned PATTERNLENGTH[] = {1,2,3,4,5,10,20,50,100,200,500,1000,2000,5000}; // size=14
    unsigned plenSize = sizeof(PATTERNLENGTH) / sizeof(unsigned);

    Rng<MersenneTwister> rng(94792847);
    unsigned maxPos = length(str) - 15000;

    // Generate a lot of random strings!
    StringSet<TText> patterns;
    resize(patterns, NUMPATTERNS);
    for (unsigned x = 0; x < NUMPATTERNS; ++x)
    {
        unsigned p = pickRandomNumber(rng) % maxPos;
        patterns[x] = prefix(suffix(str, p), 10000);
    }

    // Ungapped Index
    {
        Index<TText, IndexSa<> > index(str);
        indexRequire(index, FibreSA());

        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearchUngapped(index, patterns, PATTERNLENGTH[x]);
            std::cout << "ungapped\tmlr\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }

        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearchUngapped_stree(index, patterns, PATTERNLENGTH[x]);
            std::cout << "ungapped\tstree\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }
    }

    // Fixed Shape 110
    {
        typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1> >, 1> > TShape;
        Index<TText, IndexSa<Gapped<ModCyclicShape<TShape> > > > index(str);
        indexRequire(index, FibreSA());

        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearch(index, patterns, PATTERNLENGTH[x], TShape());
            std::cout << "Fixed 110\tmlr\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }

        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearch_stree(index, patterns, PATTERNLENGTH[x], TShape());
            std::cout << "Fixed 110\tstree\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }
    }

    // Fixed Shape 101110010101110110101110101011 (19/30)
    {
        typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<2,1,1,3,2,2,1,1,2,1,2,2,1,1,2,2,2,1> >, 0> > TShape;
        Index<TText, IndexSa<Gapped<ModCyclicShape<TShape> > > > index(str);
        indexRequire(index, FibreSA());

        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearch(index, patterns, PATTERNLENGTH[x], TShape());
            std::cout << "Fixed19/30\tmlr\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }

        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearch_stree(index, patterns, PATTERNLENGTH[x], TShape());
            std::cout << "Fixed19/30\tstree\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }
    }


// generic probably way too slow... won't finish before my deadline xD

/*
    // Generic Shape 110
    {
        typedef CyclicShape<GenericShape> TShape;
        TShape shape;
        stringToCyclicShape(shape, "110");
        Index<TText, IndexSa<Gapped<ModCyclicShape<TShape> > > > index(str, shape);
        indexRequire(index, FibreSA());
        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearch(index, patterns, PATTERNLENGTH[x], shape);
            std::cout << "Generic 110\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }
    }

    // Generic Shape 101110010101110110101110101011
    {
        typedef CyclicShape<GenericShape> TShape;
        TShape shape;
        stringToCyclicShape(shape, "101110010101110110101110101011");
        Index<TText, IndexSa<Gapped<ModCyclicShape<TShape> > > > index(str, shape);
        indexRequire(index, FibreSA());
        for (unsigned x = 0; x < plenSize; ++x)
        {
            double teim = sysTime();
            benchmarkSearch(index, patterns, PATTERNLENGTH[x], shape);
            std::cout << "Generic19/30\t" << PATTERNLENGTH[x] << "\t" << sysTime() - teim << std::endl;
        }
    }
 */

}

// --------------------------------------------------------------------------
// Function castAndRun()
// --------------------------------------------------------------------------

template <typename T>
void castAndRun(T & text, SATraversalOptions const & options)
{
    static bool arr[256]; // 0-initialized
    typedef typename Iterator<T, Standard>::Type TIter;
    for(TIter it=begin(text); it != end(text); ++it)
        arr[ordValue(*it)] = true;
    unsigned sigma = 0;
    for (unsigned i = 0; i< 256; ++i)
        if (arr[i]) ++sigma;

    std::cout << "# length " << length(text) << std::endl;
    std::cout << "# alphabet size: " << sigma << " (assuming 4: DNA, 5: DNA5, <=24: AminoAcid, >24: char)" << std::endl;

    if (sigma ==4) { DnaString str; move(str,text); callAllTests(str, options); }
    if (sigma ==5) { Dna5String str; move(str,text); callAllTests(str, options); }
    if (sigma >5 && sigma < 25) { Peptide str; move(str,text); callAllTests(str, options); }
    if (sigma > 24) { CharString str; move(str,text); callAllTests(str, options); }
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    SATraversalOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    options.print();


    // Read file
    std::ifstream f(toCString(options.file));
    if (!f.is_open())
    {
        std::cout << "Cannot open " << f << std::endl;
        return 1;
    }
    CharString str;
    std::string line;
    while ( getline (f,line) )
    {
        if (line[0] == '>')
            continue;
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);
        str += static_cast<CharString>(line);
        if (length(str) > 100000000) break;
    }
    if (length(str) > 100000000)
        resize(str, 100000000); // take only 100000

    // cast file to the smallest alphabet and run tests on it
    castAndRun(str, options);

    f.close();
    return 0;
}

