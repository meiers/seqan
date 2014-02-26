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

#include "seqanLast_IO.h"
#include "seqanLast_core.h"

using namespace seqan;

// =============================================================================
// Classes
// =============================================================================

// -----------------------------------------------------------------------------
// Class SeqanLastDbOptions
// -----------------------------------------------------------------------------

template <typename Shape, unsigned Q, typename Size1, typename Size2>
struct SeqanLastDbOptions
{
    int verbosity; // 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    CharString databaseFile;

    SeqanLastOptions() : verbosity(1)
    {}

    void print() // TODO(meiers): could make an operator<< out of this
    {
        std::cout << "Files:" << std::endl;
        std::cout << "   database: " << databaseFile << std::endl;
    }
};

// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
// Function adaptedCreateQGramIndexDirOnly()
// -----------------------------------------------------------------------------

template <
    typename TDir,
    typename TBucketMap,
    typename TText,
    typename TShape>
void _insertMissingQGrams(TDir &dir, TBucketMap &bucketMap, TText const &text, TShape &shape)
{
    typedef typename Size<TText>::Type TPos;
    typedef typename Iterator<TText const,Standard>::Type TIter;

    if (!length(text)) return;

    std::cout << shape.span << std::endl;
    TPos start = length(text) < shape.span ? 0 : length(text)-shape.span+1;
    TIter it = begin(text, Standard()) + start;
    ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];

    for(++start, ++it; start < length(text); ++start, ++it)
    {
        ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];
    }
}

template <
    typename TDir,
    typename TBucketMap,
    typename TText,
    typename TSetSpec,
    typename TShape>
void _insertMissingQGrams(TDir &dir, TBucketMap &bucketMap, StringSet<TText, TSetSpec> const &textSet, TShape &shape)
{
    typedef typename Size<TText>::Type TPos;
    typedef typename Iterator<TText const,Standard>::Type TIter;
    typedef typename Iterator<StringSet<TText, TSetSpec> const,Standard>::Type TSetIter;

    if (!length(textSet)) return;

    for (TSetIter set=begin(textSet,Standard()); set != end(textSet,Standard()); ++set)
    {
        TText const & text = *set;
        if (!length(text)) continue;

        TPos start = length(text) < shape.span ? 0 : length(text)-shape.span+1;
        TIter it = begin(text, Standard()) + start;
        ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];

        for(++start, ++it; start < length(text); ++start, ++it)
        {
            ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];
        }
    }
}


template <
    typename TDir,
    typename TBucketMap,
    typename TText,
    typename TShape,
    typename TStepSize >
void adaptedCreateQGramIndexDirOnly(
    TDir &dir,
    TBucketMap &bucketMap,
    TText const &text,
    TShape &shape,
    TStepSize stepSize)
{
	SEQAN_CHECKPOINT

    // 1. clear counters
    _qgramClearDir(dir, bucketMap);

    // 2. count q-grams
    _qgramCountQGrams(dir, bucketMap, text, shape, 1);

    // New part: Insert missing q-grams.
    _insertMissingQGrams(dir, bucketMap, text, shape);
    
    // 3. cumulative sum (Step 4 is ommited)
    _qgramCummulativeSumAlt(dir, False());
}

// -----------------------------------------------------------------------------
// Function lastdb()
// -----------------------------------------------------------------------------

template <typename TShape, unsigned Q, typename TStringSet>
void lastdb(TStringSet const & str, CharString const & outputName)
{
    std::cout << "Building Suffix array... " << std::endl;

}


// -----------------------------------------------------------------------------
// Function main()
// -----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{

    // only Dna5 supported
    typedef Dna5 TAlphabet;

    typedef String<TAlphabet>                            TSeq;
    typedef StringSet<TSeq, Owner<ConcatDirect<> > >     TSeqSet; // TODO: MAke this concat direct ??

    // get options
    SeqanLastDbOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // import database sequence
    TSeqSet databases;
    StringSet<CharString> databaseIds;
    if (!_importSequences(databases, databaseIds, options.databaseFile, options.verbosity))
        return 1;

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
    indexRequire(index, FibreSA());
    typedef Index<TSeqSet, IndexQGram<Shape<Dna5, UngappedShape<8> > > > TTable;
    TTable table(databases);
    resize(indexDir(table), _fullDirLength(table), Exact());
    adaptedCreateQGramIndexDirOnly(indexDir(table), indexBucketMap(table), indexText(table), indexShape(table), getStepSize(table));


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
}
