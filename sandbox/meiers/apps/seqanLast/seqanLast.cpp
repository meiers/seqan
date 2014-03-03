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
// Function main()
// -----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    typedef Dna5 TAlphabet;
    typedef String<TAlphabet>                            TSeq;
    typedef StringSet<TSeq, Owner<> >                    TSeqSet; // TODO: MAke this concat direct ??
    typedef Align<typename Infix<TSeq const>::Type>      TAlignment;
    typedef SeqanLastMatch<Size<TSeq>::Type, TAlignment> TMatch;

    // get options
    SeqanLastOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    if (!_readPropertyFile(options))
        return 1;


    // Import Index
    typedef StringSet<String<Dna5, External<> >, Owner<ConcatDirect<> > > TStringSet;
    Index<TStringSet, IndexSa<> > suffixArray;
    CharString f = options.databaseName;    append(f, ".txt");
    if (!open(indexText(suffixArray), toCString(f)))
        return 2;
    if (options.verbosity>1)
        std::cout << "Loaded database with " << length(indexText(suffixArray)) << " sequences and a total length of "
        << lengthSum(indexText(suffixArray)) << std::endl;
    if (!open(suffixArray, toCString(options.databaseName)))
        return 2;
    if (options.verbosity>1)
        std::cout << "Loaded suffix array with " << length(indexSA(suffixArray)) << " entries" << std::endl;






    return 0;
/*
    // import query sequences
    TSeqSet queries;
    StringSet<CharString> queryIds;
    if (!_importSequences(queries, queryIds, options.queryFile, options.verbosity))
        return 1;


    // Output Options:
    if(options.verbosity > 1)
        options.print();

    // Create Indices:
    if (options.verbosity) std::cout << "Building libraries..." << std::endl;
    typedef Index<TSeqSet, IndexSa<> > TIndex;
    TIndex index(databases);

    typedef Index<TSeqSet, IndexQGram<Shape<Dna5, UngappedShape<8> > > > TTable;
    TTable table(databases);

    // Prepare Scores
    Score<int, Simple> scoreMatrix(options.matchScore, options.mismatchScore, options.gapExtendScore,
                                   options.gapExtendScore + options.gapOpenScore);


    // Do the main work: alignments
    if (options.verbosity) std::cout << "Start searching..." << std::endl;
    String<TMatch> results;
    linearLastal(results, index, table, queries, options.frequency, scoreMatrix, options.gaplessXDrop,
                 options.gappedXDrop, options.gaplessThreshold, options.gappedThreshold, options.verbosity);


    // Output
    std::ofstream file;
	file.open(toCString(options.outputFile), ::std::ios_base::out | ::std::ios_base::app);
	if (options.outputFile == "" || !file.is_open()) {
        if (options.outputFile != "")
            std::cout << "Could not open \"" << options.outputFile << "\" to write output. Using stdout instead." << std::endl;
        std::cout << "================================================================================" << std::endl;
        _outputMatches(results, databaseIds, queryIds, std::cout, options.verbosity);
	} else {
        std::cout << "Writing output to " << options.outputFile << "..." << std::endl;
        _outputMatches(results, databaseIds, queryIds, file, options.verbosity);
    }
    file.close();
    
    
    return 0;
 
 */
}
