// ==========================================================================
//                         benchmark_sa_construction
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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>

#include <time.h>

using namespace seqan;


int globalWrongMethods = 0;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

struct AppOptions
{
    CharString infile;
    CharString mode;
    bool qsort;
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqanBuildSuffixArray");
    // Set short description, version, and date.
    setShortDescription(parser, "Builds a (gapped) suffix array with different methods");
    setVersion(parser, "Sascha.0.1");
    setDate(parser, "July 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "Builds an index (ESA) of the specified genome and writes it to disk.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "IN"));

    addOption(parser, seqan::ArgParseOption("n", "noqsort", "Disable QuickSort (for large input files)"));
    addOption(parser, seqan::ArgParseOption("m", "mode", "What kind of benchmark?",
                                            seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "mode", "correct external runtime linear medium");
    setDefaultValue(parser, "m", "correct");

    addTextSection(parser, "Modes");
    addText(parser, "Correct mode (-m correct): Checks whether different Methods get the "
                    "same result. Runs all exisiting methods.");
    addText(parser, "External mode (-m external): Only runs the external versions of Skew7, "
                    "Skew3 and Dislex. You can monitor the external memory consumption "
                    "using lsof by giving the PID to lsof and grepping for the tmp directory");
    addText(parser, "Runtime mode (-m runtime): Run all methods and measure their runtime. "
                    "Does not specify whether the internal or external algorithms are "
                    "called, but this happens in all modi except for external");
    addText(parser, "Linear mode (-m linear): Like runtime mode, but only for linear "
                    "time methods, so Skew and dislex. ");
    addText(parser, "Medium mode (-m medium): Somewhere between runtime mode and linear mode. "
                    "It runs the linear methods and inplace radix sort. This can be used on "
                    "larger files than runtime, but those should not be too repeat- (or N)-rich "
                    "as this usually kills radix Sort");



    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // disable qicksort
    options.qsort = ! isSet(parser, "noqsort");

    // set outfile to infile.db if not outfile was specified
    seqan::String<char> outfile;
    seqan::getArgumentValue(options.infile, parser, 0);
    seqan::getOptionValue(options.mode, parser, "mode");

    outfile = options.	infile;
    outfile += ".index";

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function createAndCheckSACA()
// --------------------------------------------------------------------------

template <typename TIndex, typename TAlgo, typename TSuffAr, typename TLabel>
void createAndCheckSACA(
                   TIndex & index,
                   TAlgo const & algo,
                   TSuffAr const & correctSA,
                   TLabel const & labelPattern,
                   TLabel const & labelAlgorithm,
                   TLabel const & labelText)
{
    double start = sysTime();
    indexCreate(index, EsaSA(), algo);

    std::cout << labelAlgorithm << "\t" << labelPattern << "\t" << labelText << "\t" << sysTime() - start << "\t";

    if(length(correctSA) == length(indexSA(index)) && length(indexText(index)) > 0)
    {
        typename Size<TIndex>::Type errs =0;
        typename Iterator<typename Fibre<TIndex, EsaSA>::Type>::Type it2 = begin(indexSA(index));
        for(typename Iterator<TSuffAr const>::Type it=begin(correctSA); it < end(correctSA); ++it, ++it2)
        //if(*it != *it2) std::cout << std::endl << "   |-> err " << ++errs << " at pos " << it - begin(correctSA) << ": correct " << *it << " , seen " << *it2;
        if(*it != *it2) ++errs;
        std::cout << "errors:" << errs << std::endl;
        if (errs > 0)
            ++globalWrongMethods;
    } else {
        std::cout << "/" << std::endl;
    }
}

// --------------------------------------------------------------------------
// Function externalBenchmark()
// --------------------------------------------------------------------------

template <typename TIndex, typename TAlgo, typename TLabel>
void externalBenchmark( TIndex & index,
                        TAlgo const & algo,
                        TLabel const & labelPattern,
                        TLabel const & labelAlgorithm,
                        TLabel const & labelText)
{
    double start = sysTime();
    _createSuffixArrayPipelining(indexSA(index), indexText(index), algo);
    std::cout << labelAlgorithm << "\t" << labelPattern << "\t" << labelText << "\t" << sysTime() - start << std::endl;
}

template <typename TIndex, typename TAlgo, typename TLabel>
void runtime(  TIndex & index,
               TAlgo const & algo,
               TLabel const & labelPattern,
               TLabel const & labelAlgorithm,
               TLabel const & labelText)
{
    double start = sysTime();
    indexCreate(index, EsaSA(), algo);
    std::cout << labelAlgorithm << "\t" << labelPattern << "\t" << labelText << "\t" << sysTime() - start << std::endl;
}


// --------------------------------------------------------------------------
// Function callBenchmarksExternal()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void callBenchmarksExternal(StringSet<TString, TSpec> const & set) {

    // header
    std::cout << "# Mode: External" << std::endl << std::endl;
    std::cout <<  "algorithm         pattern          text    runtime" << std::endl;

    // labels for output
    CharString ungapped  = "ungapped______";
    CharString shape     = "11000100101110";
    CharString shape111  = "10x1          ";
    CharString skew7     = "Skew7         ";
    CharString skew3     = "Skew3         ";
    CharString dislex7   = "Dislex + Skew7";
    CharString dislex3   = "Dislex + Skew3";
    CharString skewLex   = "Skew7 lexText ";
    CharString string    = "String";
    CharString stringset = "StrSet";

    {   //- Ungapped Indices ----------------------------------------------------------------
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            externalBenchmark(index, Skew7(), ungapped, skew7, stringset);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            externalBenchmark(index, Skew7(), ungapped, skew7, string);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            externalBenchmark(index, Skew3(), ungapped, skew3, string);
        }
    }

    { //----- Gapped Indices: 11000100101110 ----------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            externalBenchmark(index, DislexExternal<TShape, Skew7>(), shape, dislex7, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            externalBenchmark(index, DislexExternal<TShape, Skew7>(), shape, dislex7, string);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            externalBenchmark(index, DislexExternal<TShape, Skew3>(), shape, dislex3, string);
        }
        {
            typedef String<typename Size<TString>::Type> TLexText;
            TLexText SA, lexText;

            resize(SA, length(concat(set)));
            _initializeSA(SA, concat(set));

            inplaceRadixSort(SA, concat(set), WEIGHT<TShape>::VALUE +1, TShape(), ModCyclicShape<TShape>());

            // disLexTransformation
            // TODO(meiers): Without the static_cast (copy!!) there is an problem in instantiation of the modified iterator
            _dislex(lexText, SA, static_cast<TString>(concat(set)), TShape());
            clear(SA);

            typedef Index<TLexText const, IndexSa<> > TIndex;
            TIndex index(lexText);
            externalBenchmark(index, Skew7(), shape, skewLex, string);
        }
    }
    { //----- Gapped Indices: 1111111111 ----------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,1,1,1,1,1,1,1,1> >, 0> > TShape;
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            externalBenchmark(index, DislexExternal<TShape, Skew7>(), shape111, dislex7, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            externalBenchmark(index, DislexExternal<TShape, Skew7>(), shape111, dislex7, string);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            externalBenchmark(index, DislexExternal<TShape, Skew3>(), shape111, dislex3, string);
        }
        {
            typedef String<typename Size<TString>::Type> TLexText;
            TLexText SA, lexText;

            resize(SA, length(concat(set)));
            _initializeSA(SA, concat(set));

            inplaceRadixSort(SA, concat(set), WEIGHT<TShape>::VALUE +1, TShape(), ModCyclicShape<TShape>());

            // disLexTransformation
            // TODO(meiers): Without the static_cast (copy!!) there is an problem in instantiation of the modified iterator
            _dislex(lexText, SA, static_cast<TString>(concat(set)), TShape());
            clear(SA);

            typedef Index<TLexText const, IndexSa<> > TIndex;
            TIndex index(lexText);
            externalBenchmark(index, Skew7(), shape111, skewLex, string);
        }
    }
}






// --------------------------------------------------------------------------
// Function callBenchmarksForCorrectness()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void callBenchmarksForCorrectness(StringSet<TString, TSpec> const & set) {

    // header
    std::cout << "# Mode: Correctness check" << std::endl << std::endl;
    std::cout <<  "algor.\tpattern_______\ttext   \truntime\tcheck" << std::endl;



    // labels for output
    CharString ungapped  = "ungapped______";
    CharString shape101  = "101___________";
    CharString shape2    = "11000100101110";
    CharString shape3    = "0001010_______";
    CharString skew7     = "Skew7";
    CharString skew3     = "Skew3";
    CharString saqsort   = "QSort";
    CharString dislex    = "Dislex";
    CharString dislexExt = "Extern";
    CharString radix     = "Radix";
    CharString string    = "String";
    CharString stringset = "StrSet";


    std:: cout << "Whole String Set" << std::endl;
    {   //- Ungapped Indices ----------------------------------------------------------------

        String<Pair<typename Size<TString>::Type> > correctSA1;
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, Skew7(), correctSA1, ungapped, skew7, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, Skew3(), correctSA1, ungapped, skew3, stringset);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, SAQSort(), correctSA1, ungapped, saqsort, stringset);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, InplaceRadixSort(), correctSA1, ungapped, radix, stringset);
        }
    }
    {   //- Gapped Indices: 101 ----------------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<2> >, 0> > TShape;
        String<Pair<typename Size<TString>::Type> > correctSA1;
        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, InplaceRadixSort(), correctSA1, shape101, radix, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, SAQSort(), correctSA1, shape101, saqsort, stringset);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, Dislex<Skew7>(), correctSA1, shape101, dislex, stringset);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, DislexExternal<TShape>(), correctSA1, shape101, dislexExt, stringset);
        }
    }
    { //----- Gapped Indices: 11000100101110 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
        String<Pair<typename Size<TString>::Type> > correctSA1;
        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, InplaceRadixSort(), correctSA1, shape2, radix, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, SAQSort(), correctSA1, shape2, saqsort, stringset);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, Dislex<Skew7>(), correctSA1, shape2, dislex, stringset);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, DislexExternal<TShape>(), correctSA1, shape2, dislexExt, stringset);
        }
    }
    { //----- Gapped Indices: 0001010 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<3,GappedShape<HardwiredShape<2> >, 1> > TShape;
        String<Pair<typename Size<TString>::Type> > correctSA1;
        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, InplaceRadixSort(), correctSA1, shape3, radix, stringset);

            correctSA1 = indexSA(index);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, SAQSort(), correctSA1, shape3, saqsort, stringset);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, Dislex<Skew7>(), correctSA1, shape3, dislex, stringset);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            createAndCheckSACA(index, DislexExternal<TShape>(), correctSA1, shape3, dislexExt, stringset);
        }
    }


    std:: cout << "Now for single strings String Set" << std::endl;
    
    {   //- Ungapped Indices ----------------------------------------------------------------

        StringSet<String<typename Size<TString>::Type> > correctSAs;
        resize(correctSAs, length(set));

        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (set[i]);
            createAndCheckSACA(index, Skew7(), correctSAs[i], ungapped, skew7, string);

            correctSAs[i] = indexSA(index);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (set[i]);
            createAndCheckSACA(index, Skew3(), correctSAs[i], ungapped, skew3, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (set[i]);
            createAndCheckSACA(index, SAQSort(), correctSAs[i], ungapped, saqsort, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, InplaceRadixSort(), correctSAs[i], ungapped, radix, string);
        }
    }

    {   //- Gapped Indices: 101 ----------------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<2> >, 0> > TShape;
        StringSet<String<typename Size<TString>::Type> > correctSAs;
        resize(correctSAs, length(set));

        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, InplaceRadixSort(), correctSAs[i], shape101, radix, string);

            correctSAs[i] = indexSA(index);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, SAQSort(), correctSAs[i], shape101, saqsort, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, Dislex<Skew7>(), correctSAs[i], shape101, dislex, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, DislexExternal<TShape>(), correctSAs[i], shape101, dislexExt, string);
        }
    }
    { //----- Gapped Indices: 11000100101110 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
        StringSet<String<typename Size<TString>::Type> > correctSAs;
        resize(correctSAs, length(set));

        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, InplaceRadixSort(), correctSAs[i], shape2, radix, string);

            correctSAs[i] = indexSA(index);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, SAQSort(), correctSAs[i], shape2, saqsort, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, Dislex<Skew7>(), correctSAs[i], shape2, dislex, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, DislexExternal<TShape>(), correctSAs[i], shape2, dislexExt, string);
        }
    }
    { //----- Gapped Indices: 0001010 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<3,GappedShape<HardwiredShape<2> >, 1> > TShape;
        StringSet<String<typename Size<TString>::Type> > correctSAs;
        resize(correctSAs, length(set));

        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, InplaceRadixSort(), correctSAs[i], shape3, radix, string);

            correctSAs[i] = indexSA(index);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, SAQSort(), correctSAs[i], shape3, saqsort, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, Dislex<Skew7>(), correctSAs[i], shape3, dislex, string);
        }
        for (unsigned i=0; i<length(set); ++i)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set[i]);
            createAndCheckSACA(index, DislexExternal<TShape>(), correctSAs[i], shape3, dislexExt, string);
        }
    }


    std::cout << std::endl << std::endl;
    std::cout << "Number of test cases that are incorrect: " << globalWrongMethods << std::endl;

}


// --------------------------------------------------------------------------
// Function callBenchmarks()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void callBenchmarks(StringSet<TString, TSpec> const & set, int level=3) {

    // header
    std::cout << "# Mode: Runtime measurment" << std::endl << std::endl;
    std::cout <<  "algor.\tpattern_______\ttext  \truntime" << std::endl;

    // labels for output
    CharString ungapped  = "ungapped______";
    CharString shape101  = "101___________";
    CharString shape2    = "11000100101110";
    CharString shape3    = "0001010_______";
    CharString skew7     = "Skew7";
    CharString skew3     = "Skew3";
    CharString saqsort   = "QSort";
    CharString dislex    = "Dislex";
    CharString dislexExt = "Extern";
    CharString radix     = "Radix";
    CharString string    = "String";
    CharString stringset = "StrSet";


    {   //- Ungapped Indices ----------------------------------------------------------------

        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            runtime(index, Skew7(), ungapped, skew7, stringset);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            runtime(index, Skew7(), ungapped, skew7, string);
        }
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            runtime(index, Skew3(), ungapped, skew3, stringset);
        }
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            runtime(index, Skew3(), ungapped, skew3, string);
        }
        if (level>2)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            runtime(index, SAQSort(), ungapped, saqsort, stringset);
        }
        if (level>2)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index (concat(set));
            runtime(index, SAQSort(), ungapped, saqsort, string);
        }
        if (level>1)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<> > TIndex;
            TIndex index(set);
            runtime(index, InplaceRadixSort(), ungapped, radix, stringset);
        }
        if (level>1)
        {
            typedef Index<TString, IndexSa<> > TIndex;
            TIndex index(concat(set));
            runtime(index, InplaceRadixSort(), ungapped, radix, string);
        }
    }
    {   //- Gapped Indices: 101 ----------------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<2> >, 0> > TShape;

        if (level>1)
        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, InplaceRadixSort(), shape101, radix, stringset);
        }
        if (level>1)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, InplaceRadixSort(), shape101, radix, string);
        }
        if (level>2)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, SAQSort(), shape101, saqsort, stringset);
        }
        if (level>2)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, SAQSort(), shape101, saqsort, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, Dislex<Skew7>(), shape101, dislex, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, Dislex<Skew7>(), shape101, dislex, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, DislexExternal<TShape>(), shape101, dislexExt, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, DislexExternal<TShape>(), shape101, dislexExt, string);
        }
    }


    { //----- Gapped Indices: 11000100101110 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;

        if (level>1)
        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, InplaceRadixSort(), shape2, radix, stringset);
        }
        if (level>1)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, InplaceRadixSort(), shape2, radix, string);
        }
        if (level>2)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, SAQSort(), shape2, saqsort, stringset);
        }
        if (level>2)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, SAQSort(), shape2, saqsort, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, Dislex<Skew7>(), shape2, dislex, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, Dislex<Skew7>(), shape2, dislex, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, DislexExternal<TShape>(), shape2, dislexExt, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, DislexExternal<TShape>(), shape2, dislexExt, string);
        }

    }

    {   //----- Gapped Indices: 0001010 ---------------------------------------------------------

        typedef CyclicShape<FixedShape<3,GappedShape<HardwiredShape<2> >, 1> > TShape;

        if (level>1)
        {
            typedef Index<StringSet<TString, TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, InplaceRadixSort(), shape3, radix, stringset);
        }
        if (level>1)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, InplaceRadixSort(), shape3, radix, string);
        }
        if (level>2)
        {
            typedef Index<StringSet<TString, TSpec> const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, SAQSort(), shape3, saqsort, stringset);
        }
        if (level>2)
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, SAQSort(), shape3, saqsort, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, Dislex<Skew7>(), shape3, dislex, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, Dislex<Skew7>(), shape3, dislex, string);
        }
        {
            typedef Index<StringSet<TString,TSpec>, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(set);
            runtime(index, DislexExternal<TShape>(), shape3, dislexExt, stringset);
        }
        {
            typedef Index<TString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
            TIndex index(concat(set));
            runtime(index, DislexExternal<TShape>(), shape3, dislexExt, string);
        }
    }

}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;


    // start main program ///////////////////////////////////////////////////

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;



    // New reading:
    // use the sequenceStream
    seqan::SequenceStream seqStream(seqan::toCString(options.infile));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file." << std::endl;
        return 1;
    }
    if (readAll(ids, seqs, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from " << options.infile << std::endl;
        return 1;
    }
    

    time_t t; time(&t);
    std::cout << "############# Call Benchmarks for Seqan's SACAs ###############" << std::endl;
    std::cout << "# input file : " << options.infile << std::endl;
    std::cout << "# timestamp  : " << ctime(&t);
    std::cout << "# alphabet size: " << (int)ValueSize<Dna5>::VALUE << std::endl;
    std::cout << "# dataset: " << length(seqs) << " strings, total length " << length(concat(seqs)) << std::endl;


    if (options.mode == "correct")
        callBenchmarksForCorrectness(seqs);
    if (options.mode == "external")
        callBenchmarksExternal(seqs);
    if (options.mode == "runtime")
        callBenchmarks(seqs, 3);
    if (options.mode == "medium")
        callBenchmarks(seqs, 2);
    if (options.mode == "linear")
        callBenchmarks(seqs, 1);

    
    return 0;
}