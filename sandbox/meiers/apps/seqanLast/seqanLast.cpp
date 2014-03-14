// ==========================================================================
//                                 seqanLast
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
// Parts of this code are copied (and or modified) from the Stellar source code
// by Birte Kehr, others might be taken from Martin Friths LAST code
// (last.cbrc.jp).

#include <cstdlib>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/align_extend.h>
#include "seqanLast_IO.h"
#include "seqanLast_core.h"

//#include <seqan/align_extend.h>

using namespace seqan;


// =============================================================================
// Classes
// =============================================================================

// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
// Function importAndRun()
// -----------------------------------------------------------------------------

template <unsigned K, typename TShape>
int importAndRun(SeqanLastOptions &options,
                 TShape const &)
{
    // typedefs happen here (generic functions follow)
    typedef Align<typename Infix<TMMapStringSet const>::Type>       TAlignment;
    typedef SeqanLastMatch<Size<TMMapStringSet>::Type, TAlignment>  TMatch;


    // Import Suffix Array
    Index<TMMapStringSet, IndexSa<Gapped<ModCyclicShape<TShape> > > > suffixArray;          // Index Type
    if (!open(suffixArray, toCString(options.databaseName)))
        return 1;
    if (options.verbosity>1)
        std::cout << "Loaded suffix array with " << length(indexSA(suffixArray)) <<
        " entries" << std::endl;


    // Import Database sequences
    CharString f = options.databaseName;    append(f, ".txt");
    if (!open(indexText(suffixArray), toCString(f)))
        return 1;
    if (options.verbosity>1)
        std::cout << "Loaded database with " << length(indexText(suffixArray)) <<
        " sequences and a total length of " << lengthSum(indexText(suffixArray)) << std::endl;


    // Import Database ids
    StringSet<CharString> dbIds;
    if (!_readIdFile(dbIds, options) || length(dbIds) != length(indexText(suffixArray)))
        return 1;
    if (options.verbosity>1)
        std::cout << "Loaded database ids." << std::endl;


    // Load Q-Gram table
    Index<TMMapStringSet, IndexQGram<UngappedShape<K> > > hashTab;
    f = options.databaseName;    append(f, ".dir");
    if (!open(indexDir(hashTab), toCString(f)))
        return 1;
    if (options.verbosity>1)
        std::cout << "Loaded hashTable with " << length(indexDir(hashTab)) << " entries" << std::endl;


    // Import Query
    TMMapStringSet querySet;
    StringSet<CharString> queryIds;
    if (!_importSequences(querySet, queryIds, options.queryFile, options.verbosity))
        return 1;


    // Prepare Scores
    Score<int, Simple> scoreMatrix(options.matchScore, options.mismatchScore, options.gapExtendScore,
                                   options.gapExtendScore + options.gapOpenScore);

    
    

    // Do the main work: alignments
    if (options.verbosity) std::cout << "Start searching..." << std::endl;
    String<TMatch> matchContainer;
    linearLastal(matchContainer,
                 suffixArray,
                 hashTab,
                 querySet,
                 options.frequency,
                 scoreMatrix,
                 options.gaplessXDrop,
                 options.gappedXDrop,
                 options.gaplessThreshold,
                 options.gappedThreshold,
                 options.verbosity);


    // Output
    std::ofstream file;
	file.open(toCString(options.outputFile), ::std::ios_base::out);
	if (options.outputFile == "stdout" || !file.is_open()) {
        if (options.outputFile != "stdout")
            std::cout << "Could not open \"" << options.outputFile << "\" to write output. Using stdout instead." << std::endl;
        std::cout << "================================================================================" << std::endl;
        _outputMatches(matchContainer, dbIds, queryIds, std::cout, options.verbosity);
	} else {
        std::cout << "Writing output to " << options.outputFile << "..." << std::endl;
        _outputMatches(matchContainer, dbIds, queryIds, file, options.verbosity);
    }
    file.close();

    std::cout << "ENDE" << std::endl;
    return 0;
}

// -----------------------------------------------------------------------------
// Function _lastChoice()
// -----------------------------------------------------------------------------

// 2: Choose k
template <typename TShape>
int _lastChoice2(SeqanLastOptions &options,
                 TShape const &)
{
    switch (options.k)
    {
        case 2:  return importAndRun <2> (options, TShape());
        case 3:  return importAndRun <3> (options, TShape());
        case 4:  return importAndRun <4> (options, TShape());
        case 5:  return importAndRun <5> (options, TShape());
        case 6:  return importAndRun <6> (options, TShape());
        case 7:  return importAndRun <7> (options, TShape());
        case 8:  return importAndRun <8> (options, TShape());
        case 9:  return importAndRun <9> (options, TShape());
        case 10: return importAndRun <10>(options, TShape());
        case 11: return importAndRun <11>(options, TShape());
        case 12: return importAndRun <12>(options, TShape());
        default:
            std::cerr << "No valid k-mer size chosen. Exit" << std::endl;
            return 2;
    }
}

// 1: Choose Shape
int _lastChoice1(SeqanLastOptions &options)
{
    switch (options.shapeChoice)
    {
        case 1: return _lastChoice2(options, Shape1() );
        case 2: return _lastChoice2(options, Shape2() );
        default:
            std::cerr << "No valid shape chosen. Exit" << std::endl;
            return 1;
    }
}


// -----------------------------------------------------------------------------
// Function main()
// -----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // get options
    SeqanLastOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Output Options:
    if(options.verbosity > 1)
        options.print();

    // Load property file
    if (!_readPropertyFile(options))
        return 1;

    return _lastChoice1(options);
}
