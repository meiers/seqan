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
    typedef StringSet<String<TAlphabet>, Owner<> > TSeqSet; // TODO: MAke this concat direct ??
    typedef SeqanLastMatch<Size<TSeqSet>::Type, Gaps<Infix<String<TAlphabet> >::Type, ArrayGaps> > Match;

    // get options
    SeqanLastOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // import database sequence
    TSeqSet databases;
    StringSet<CharString> databaseIDs;
    if (!_importSequences(databases, databaseIDs, options.databaseFile, options.verbosity))
        return 1;

    // import query sequences
    TSeqSet queries;
    StringSet<CharString> queryIDs;
    if (!_importSequences(queries, queryIDs, options.queryFile, options.verbosity))
        return 1;


    // Get Indices
    typedef Index<TSeqSet, IndexSa<> > TIndex;
    TIndex index(databases);
    typedef Index<TSeqSet, IndexQGram<Shape<AminoAcid, UngappedShape<2> > > > TTable;
    TTable table(databases);
    // TODO: Make hash table suit the SA


    // Prepare Scores
    Score<int, Simple> scoreMatrix(options.matchScore, options.mismatchScore, options.gapExtendScore,
                                   options.gapExtendScore + options.gapOpenScore);

    // Output Options:
    if(options.verbosity > 1)
        options.print();


    // Do the main work: alignments
    String<Match> results;
    linearLastal(results, index, table, queries, options.frequency, scoreMatrix, options.gaplessXDrop,
                 options.gappedXDrop, options.gaplessThreshold, options.gappedThreshold, options.verbosity);
    
    
    
    return 0;
}
