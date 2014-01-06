// ==========================================================================
//                               seqanLast_IO.h
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
// =============================================================================
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// =============================================================================

#ifndef CORE_APPS_SEQANLAST_SEQANLAST_IO_H_
#define CORE_APPS_SEQANLAST_SEQANLAST_IO_H_
using namespace seqan;

// =============================================================================
// Forwards
// =============================================================================

// =============================================================================
// Tags, Classes, Enums
// =============================================================================

// -----------------------------------------------------------------------------
// Class SeqanLastOptions
// -----------------------------------------------------------------------------

struct SeqanLastOptions
{
    int verbosity; // 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.

    unsigned frequency;
    int matchScore;
    int mismatchScore;
    int gapOpenScore;
    int gapExtendScore;
    int gappedXDrop;
    int gaplessXDrop;
    int gaplessThreshold;
    int gappedThreshold;

    CharString databaseFile;
    CharString queryFile;

    SeqanLastOptions() : verbosity(1)
    {}
};

// =============================================================================
// Metafunctions
// =============================================================================

// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
// Function _setParser()
// -----------------------------------------------------------------------------

void _setParser(ArgumentParser & parser)
{
    setShortDescription(parser, "Seqan version of the LAST aligner");
    setDate(parser, "February 2014");
    setVersion(parser, "0.1");
    setCategory(parser, "Local Alignment");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIFASTA FILE 1\\fP> <\\fIFASTA FILE 2\\fP>");

    addDescription(parser,
                   "SeqanLast is a reimplemtation of the core functionality of the LAST aligner "
                   "developed by Martin Frith (last.cbrc.jp, Kielbasa et al 2008: \"Adaptive "
                   "seeds tame genomic sequence comparison\"). It uses adaptive gapped seeds "
                   "to find local similarities, which are then verified using certain heuristics.");
    addDescription(parser,
                   "Right now: Input two FASTA files. Later: Input a database and one (query) FASTA."
                   "TODO: (1) supports ONLY DNA right now. (2) Default values are not chosen with any "
                   "thought. (3) Reading Database instead of 2 fasta files. (4) Output options");
    addDescription(parser, "(c) 2013-2014 by Sascha Meiers");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE 1"));
    setValidValues(parser, 0, "fa fasta");  // allow only fasta files as input
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE 2"));
    setValidValues(parser, 1, "fa fasta");  // allow only fasta files as input


    addOption(parser, ArgParseOption("v", "verbose", "Set verbosity mode."));
    addOption(parser, ArgParseOption("V", "very-verbose", "Set stronger verbosity mode."));
    addOption(parser, ArgParseOption("Q", "quiet", "No output, please."));

    addSection(parser, "Main Options");

    addOption(parser, ArgParseOption("F", "frequency", "Initial match frequency for adaptive seeds",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "F", "10");
                      setMinValue(parser, "F", "0");
    addOption(parser, ArgParseOption("r", "matchScore", "Match score for verification phase",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "r", "1");
    addOption(parser, ArgParseOption("q", "mismatchScore", "Mismatch score for verification phase",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "q", "-1");
    addOption(parser, ArgParseOption("a", "gapOpenScore", "Gap opening score for verification phase",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "a", "-7");
    addOption(parser, ArgParseOption("b", "gapExtendScore", "Gap extension score for verification phase. "
                      "Gaps of length k are scored with a+k*b", ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "b", "-1");
    addOption(parser, ArgParseOption("x", "gappedXDrop", "Maximal x-drop score during the gapped alignment",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "x", "10");
    addOption(parser, ArgParseOption("y", "gaplessXDrop", "Maximal x-drop score during the gapless extension",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "y", "14");
    addOption(parser, ArgParseOption("d", "gaplessThreshold", "Minimum score of the gapless alignment to continue",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "d", "12");
    addOption(parser, ArgParseOption("e", "gappedThreshold", "Minimum score of the gapped alignment to continue",
                      ArgParseArgument::INTEGER));
                      setDefaultValue(parser, "e", "25");

    addTextSection(parser, "References");
    addText(parser, "Kielbasa et al.: Adaptive seeds tame genomic sequence comparison. 2008");
}

// -----------------------------------------------------------------------------
// Function parseCommandLine()
// -----------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(SeqanLastOptions & options, int argc, char const ** argv)
{
    // initialise a parser for seqanLast and read in the command line
    ArgumentParser parser;
    _setParser(parser);
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract mandatory values.
    getArgumentValue(options.databaseFile, parser, 0);
    getArgumentValue(options.queryFile, parser, 1);

    // Extract optional values.
    getOptionValue(options.frequency, parser, "frequency");
    getOptionValue(options.matchScore, parser, "matchScore");
    getOptionValue(options.mismatchScore, parser, "mismatchScore");
    getOptionValue(options.gapOpenScore, parser, "gapOpenScore");
    getOptionValue(options.gapExtendScore, parser, "gapExtendScore");
    getOptionValue(options.gappedXDrop, parser, "gappedXDrop");
    getOptionValue(options.gaplessXDrop, parser, "gaplessXDrop");
    getOptionValue(options.gaplessThreshold, parser, "gaplessThreshold");
    getOptionValue(options.gappedThreshold, parser, "gappedThreshold");

    // Extract verbosity options.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return seqan::ArgumentParser::PARSE_OK;
}

// -----------------------------------------------------------------------------
// Function _checkUniqueId()
// -----------------------------------------------------------------------------

template <typename TId>
inline bool
_checkUniqueId(std::set<TId> & uniqueIds, TId & id)
{
    TId shortId;
    typedef typename Iterator<TId>::Type TIterator;

    TIterator it = begin(id);
    TIterator itEnd = end(id);
    while (it != itEnd && *it > 32)
    {
        appendValue(shortId, *it);
        ++it;
    }

    if (uniqueIds.count(shortId) == 0)
    {
        uniqueIds.insert(shortId);
        return 1;
    }

    return 0;
}

// -----------------------------------------------------------------------------
// Function _importSequences()
// -----------------------------------------------------------------------------

template <typename TSeqSet, typename TIdSet>
inline bool
_importSequences(TSeqSet & seqs, TIdSet & ids, CharString const & fileName, int verbosity)
{
    typedef typename Value<TSeqSet>::Type TSequence;
    typedef typename Value<TIdSet>::Type  TId;

    MultiSeqFile multiSeqFile;
    if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
    {
        if(verbosity)
            std::cerr << "Failed to open " << fileName << "." << std::endl;
        return false;
    }

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    unsigned seqCount = length(multiSeqFile);
    reserve(seqs, seqCount, Exact());
    reserve(ids, seqCount, Exact());

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    TSequence seq;
    TId idstr;
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seq, multiSeqFile[i], format);
        assignSeqId(idstr, multiSeqFile[i], format);

        idsUnique &= _checkUniqueId(uniqueIds, idstr);

        appendValue(seqs, seq, Generous());
        appendValue(ids, idstr, Generous());
    }

    if(verbosity) std::cout << "Loaded " << seqCount << " sequence" <<
        ((seqCount > 1) ? "s" : "") << " from " << fileName << std::endl;
    if (!idsUnique)
        if(verbosity) std::cerr << "WARNING: Non-unique IDs in " <<
            fileName << " ids. Output can be ambigous.\n";
    return true;
}

#endif  // #ifndef CORE_APPS_SEQANLAST_SEQANLAST_IO_H_