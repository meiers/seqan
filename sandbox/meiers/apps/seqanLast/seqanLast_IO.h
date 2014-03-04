// ============================================================================
//                               seqanLast_IO.h
// ============================================================================
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

#ifndef SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_IO_H_
#define SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_IO_H_
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

    CharString databaseName;
    CharString queryFile;
    CharString outputFile;

    // options from lastdb
    int shapeChoice;
    int k;

    SeqanLastOptions() : verbosity(1)
    {}

    void print() // TODO(meiers): could make an operator<< out of this
    {
        std::cout << "Files:" << std::endl;
        std::cout << "   database: " << databaseName << std::endl;
        std::cout << "   query:    " << queryFile << std::endl;
        std::cout << "   output:   " << outputFile << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "   frequency:        " << frequency << std::endl;
        std::cout << "   matchScore:       " << matchScore << std::endl;
        std::cout << "   mismatchScore:    " << mismatchScore << std::endl;
        std::cout << "   gapOpenScore:     " << gapOpenScore << std::endl;
        std::cout << "   gapExtendScore:   " << gapExtendScore << std::endl;
        std::cout << "   gappedXDrop:      " << gappedXDrop << std::endl;
        std::cout << "   gaplessThreshold: " << gaplessThreshold << std::endl;
        std::cout << "   gappedThreshold:  " << gappedThreshold << std::endl;
        std::cout << "LastDB options:" << std::endl;
        std::cout << "   k:                " << k << std::endl;
        std::cout << "   shapeChoice:      " << shapeChoice << std::endl;
    }
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
    setDate(parser, "March 2014");
    setVersion(parser, "0.1");
    setCategory(parser, "Local Alignment");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIDATABASE\\fP> <\\fIFASTA FILE 2\\fP>");

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

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "DATABASE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "QUERY FASTA FILE"));
    setValidValues(parser, 1, "fa fasta");  // allow only fasta files as input

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output", "File to write results into (gff format)",
                                     ArgParseArgument::OUTPUTFILE));
    setDefaultValue(parser, "o", "");
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
    getArgumentValue(options.databaseName, parser, 0);
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
    getOptionValue(options.outputFile, parser, "output");

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
    typedef typename Value<typename Value<TSeqSet>::Type>::Type  TAlph;
    typedef String<TAlph> TSequence;
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

// -----------------------------------------------------------------------------
// Function _writeMatchGff() and many more from Birte Kehrs STELLAR
// -----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// Determines the length and the number of matches of two alignment rows
template<typename TRow, typename TSize>
inline void
_analyzeAlignment(TRow const & row0, TRow const & row1, TSize & aliLen, TSize & matches) {
	TSize pos = 0;
    SEQAN_ASSERT_EQ(length(row0), length(row1));
    TSize end0 = length(row0);
    TSize end1 = length(row1);

	matches = 0;
    while ((pos < end0) && (pos < end1)) {
        if (!isGap(row0, pos) && !isGap(row1, pos)) {
            if (value(row0, pos) == value(row1, pos)) {
                ++matches;
            }
        }
        ++pos;
    }

	aliLen = _max(end0, end1);
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the identity of two alignment rows (percentage of matching positions).
template<typename TRow>
double
_computeIdentity(TRow const & row0, TRow const & row1) {
    typedef typename Size<TRow>::Type TSize;
    TSize matches, aliLen;
	_analyzeAlignment(row0, row1, aliLen, matches);

    return floor(1000000.0 * matches / aliLen) / 10000.0;
}

///////////////////////////////////////////////////////////////////////////////
// Computes a CIGAR string and mutations from rows of StellarMatch.
template<typename TRow, typename TString>
void
_getCigarLine(TRow const & row0, TRow const & row1, TString & cigar, TString & mutations) {
    SEQAN_CHECKPOINT
    typedef typename Size<TRow>::Type TSize;

    TSize pos = 0;

    SEQAN_ASSERT_EQ(length(row0), length(row1));
    TSize dbEndPos = length(row0);
    TSize queryEndPos = length(row1);

    bool first = true;
    TSize readBasePos = pos + clippedBeginPosition(row1);
    TSize readPos = 0;
	while (pos < dbEndPos || pos < queryEndPos) {
		int matched = 0;
		int inserted = 0;
		int deleted = 0;
		while (pos != dbEndPos && pos != queryEndPos &&
               !isGap(row0, pos) && !isGap(row1, pos)) {
            ++readPos;
			if (value(row0, pos) != value(row1, pos)) {
				if (first) first = false;
				else mutations << ",";
				mutations << readPos << value(source(row1), readBasePos);
			}
			++readBasePos;
			++pos;
			++matched;
		}
		if (matched > 0) cigar << matched << "M" ;
		while (pos < dbEndPos && isGap(row1, pos)) {
			++pos;
			++deleted;
		}
		if (deleted > 0) cigar << deleted << "D";
		while (pos < queryEndPos && isGap(row0, pos)) {
			++pos;
			++readPos;
			if (first) first = false;
			else mutations << ",";
			mutations << readPos << value(source(row1), readBasePos);
			++readBasePos;
			++inserted;
		}
		if (inserted > 0) cigar << inserted << "I";
	}
}

template<typename TId, typename TScore, typename TRow, typename TFile>
void
_writeMatchGff(TId const & databaseID,
               TId const & patternID,
               bool databaseStrand,
               TScore score, //lenAdjustment
               TRow const & row0,
               TRow const & row1,
               TFile & file)
{
    typedef typename Value<typename Source<TRow>::Type>::Type TAlphabet;

    for (typename Position<TId>::Type i = 0; i < length(databaseID) && value(databaseID, i) > 32; ++i) {
        file << value(databaseID, i);
    }

    file << "\tseqanLast";
    file << "\tseqanLast_match";

    if (databaseStrand) {
        file << "\t" << beginPosition(row0) + beginPosition(source(row0)) + 1;
        file << "\t" << endPosition(row0) + beginPosition(source(row0));
    } else {
        file << "\t" << length(source(row0)) - (endPosition(row0) + beginPosition(source(row0))) + 1;
        file << "\t" << length(source(row0)) - (beginPosition(row0) + beginPosition(source(row0)));
    }

    file << "\t" << _computeIdentity(row0, row1);

    file << "\t" << (databaseStrand ? '+' : '-');

    file << "\t.\t";
    for (typename Position<TId>::Type i = 0; i < length(patternID) && value(patternID, i) > 32; ++i) {
        file << value(patternID, i);
    }

    file << ";score=" << score;

    file << ";seq2Range=" << beginPosition(row1) + beginPosition(source(row1)) + 1;
    file << "," << endPosition(row1) + beginPosition(source(row1));

    //    if (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
    //        file << ";eValue=" << _computeEValue(row0, row1, lengthAdjustment);

    std::stringstream cigar, mutations;
    _getCigarLine(row0, row1, cigar, mutations);
    file << ";cigar=" << cigar.str();
    file << ";mutations=" << mutations.str();
    file << "\n";
}



///////////////////////////////////////////////////////////////////////////////
// Calls _writeMatchGff for each match in StringSet of String of matches.
//   = Writes matches in gff format to a file.
template<typename TMatches, typename TIds, typename TOutStream>
bool _outputMatches(TMatches const & matches,
    TIds const & dbIds,
    TIds const & quIds,
    TOutStream & stream,
    int verbosity)
{
    typedef typename Value<TMatches>::Type TMatch;
    typedef typename TMatch::Type TAlign;
    typedef typename Value<typename Rows<TAlign>::Type>::Type TRow;
    typedef typename Source<TRow>::Type TSeq;
    typedef typename Size<TSeq>::Type TSize;
    for(typename Iterator<TMatches const, Standard>::Type it = begin(matches); it!= end(matches); ++it)
    {
        _writeMatchGff(dbIds[it->dbId], quIds[it->quId], true, it->score, row(it->align,0), row(it->align,1), stream);
    }

    if (!verbosity)
        std::cout << "# matches = " << length(matches) << std::endl;
}




bool _readPropertyFile(SeqanLastOptions & options)
{
    CharString fileName = options.databaseName; append(fileName, ".prt");
    std::ifstream file(toCString(fileName));

    bool b_shape=false,b_k=false;

    std::string line;
    if (file.is_open())
    {
        while ( getline(file,line) )
        {
            std::stringstream l(line);
            std::string key;
            getline(l, key, '=');
            std::string value;
            getline(l, value, '=');

            if (key == "shape")
            {
                if (b_shape) {
                    if(options.verbosity>0) std::cerr << "Double definition of shape in property file. Use the first one." << std::endl;
                    continue;
                }
                b_shape = true;
                int v = std::atoi(value.c_str());
                if (v <=0 || v >2)
                {

                } else {
                    options.shapeChoice = v;
                }
            }
            else if (key == "k")
            {
                if (b_k) {
                    if(options.verbosity>0) std::cerr << "Double definition of k-mer size in property file. Use the first one." << std::endl;
                    continue;
                }
                b_k = true;
                int v = std::atoi(value.c_str());
                if (v <=0 || v >12)
                {
                    if(options.verbosity>0) std::cerr << "Strange k-mer size in property file: " << line << std::endl;
                    return false;
                } else {
                    options.k = v;
                }
            }
            else {
                if(options.verbosity>0) std::cerr << "Unknown key \"" << key << "\" in property file. Ignore it." << std::endl;
                continue;
            }
        }
        file.close();
        return true;
    }
    else
    {
        if(options.verbosity>0) std::cerr << "Cannot read property file \"" << fileName << "\" of the database." << std::endl;
        return false;
    }
}

template <typename TIdSet>
bool _readIdFile(TIdSet & set, SeqanLastOptions & options)
{
    CharString f = options.databaseName;    append(f, ".ids");
    std::ifstream file(toCString(f));
    std::string line;
    if (file.is_open())
    {
        while ( getline(file,line) )
        {
            appendValue(set, line, Generous());
        }
        return true;
    }
    else
    {
        if(options.verbosity>0) std::cerr << "Cannot read ID file \"" << f << "\" of the database." << std::endl;
        return false;
    }
}

#endif  // #ifndef SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_IO_H_
