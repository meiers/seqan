// ==========================================================================
//                              compare_matches
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: bkehr, edited by Sascha Meiers
// ==========================================================================

#include <iostream>
#include <fstream>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>

using namespace seqan;


//////////////////////////////////////////////////////////////////////////////
// Tags

struct TagLastZ_;
typedef Tag<TagLastZ_> const LastZ;

struct TagLast_;
typedef Tag<TagLast_> const Last;

struct TagBlast_;
typedef Tag<TagBlast_> const Blast;

//struct TagBlat_;
//typedef Tag<TagBlat_> const Blat;

struct TagGff_;
typedef Tag<TagGff_> const Gff;


template <typename TPosition>
struct _Match 
{
	typedef TPosition TPos;

	CharString seq1;
	CharString seq2;
    int score;                  // new

	TPos begin1;
	TPos begin2;
	TPos end1;
	TPos end2;

	bool strand;

	_Match() {
		seq1 = "";
		seq2 = "";
		strand = true;
	}
};

template <typename TPos>
int scoreDiff(_Match<TPos> const & a, _Match<TPos> const & b)
{
    return abs(a.score - b.score);
}

template <typename TPos>
int distance(_Match<TPos> const & a, _Match<TPos> const & b)
{
    return abs(a.begin1 - b.begin1) + abs(a.begin2 - b.begin2) + abs(a.end1 - b.end1) + abs(a.end2 - b.end2);
}

template<typename TPos>
void
write(_Match<TPos> & match) {
	std::cout << match.seq1 << " < " << match.begin1 << " , " << match.end1 << " > ";
	std::cout << match.seq2 << " < " << match.begin2 << " , " << match.end2 << " >";
    if (match.score >0)
        std::cout << "  score: " << match.score;
    std::cout << std::endl;
}

template<typename TReader, typename TMatch>
bool
readRecord(TReader & reader, TMatch & match, Gff) {
    CharString buffer;
    int res;

    skipWhitespaces(reader);

    // read column 1: name of sequence 1
    res = readUntilWhitespace(match.seq1, reader);
    if (res) return res;

    skipWhitespaces(reader);

	// skip columns 2 and 3
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

	// read column 4: begin position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin1 = lexicalCast<unsigned>(buffer) - 1;

    skipWhitespaces(reader);

	// read column 5: end position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end1 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // skip column 6
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

	// read column 7: strand
	if (value(reader) == '+') match.strand = true;
	else if (value(reader) == '-') match.strand = false;
	else {
		skipLine(reader);
		return 1;
	}
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

	// skip column 8
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

	// read name of sequence 2
    res = readUntilChar(match.seq2, reader, ';');
    if (res) return res;


    bool found_seq2Range = false;
    bool found_score = false;
    match.score = -999;

    // skip until key "seq2Range" or "length"
    while (!atEnd(reader) )
    {
        if (buffer == "seq2Range")
        {
            clear(buffer);
            res = readDigits(buffer, reader);
            if (res) return res;
            match.begin2 = lexicalCast<unsigned>(buffer) -1;

            // skip ','
            skipChar(reader, ',');

            // read end position 2
            clear(buffer);
            res = readDigits(buffer, reader);
            if (res) return res;
            match.end2 = lexicalCast<unsigned>(buffer);

            found_seq2Range = true;
        }
        else if (buffer == "length" && !found_seq2Range)
        {
            match.begin2 = 1;

            // length
            clear(buffer);
            res = readDigits(buffer, reader);
            if (res) return res;
            match.end2 = lexicalCast<unsigned>(buffer);
        }
        else if (buffer == "score" && !found_score)
        {
            clear(buffer);
            res = readDigits(buffer, reader);
            if (res) return res;
            match.score = lexicalCast<unsigned>(buffer);

            found_score = true;
        }
        // if both was found, leave
        if (found_seq2Range && found_score)
        {
            skipLine(reader);
            return 0;
        }

        skipUntilChar(reader, ';');
        skipChar(reader, ';');

        clear(buffer);
        res = readUntilChar(buffer, reader, '=');
        if (res) return res;

        res = skipChar(reader, '=');
        if (res) return res;
    }
	return 1;
}

template<typename TReader, typename TMatch>
bool
readRecord(TReader & reader, TMatch & match, Blast) {
    CharString buffer;
    int res;

    match.score = -999;
    skipWhitespaces(reader);

    // skip comment lines and read query sequence identifier
	while (value(reader) == '#') {
        res = skipChar(reader, '#');
        if (res) return res;

        skipWhitespaces(reader);

        clear(buffer);
        res = readUntilWhitespace(buffer, reader);
        if (res) return res;

		if (buffer == "Subject:") {
            skipWhitespaces(reader);

            // read query sequence identifyer
            res = readUntilWhitespace(match.seq1, reader);
            if (res) return res;

            skipLine(reader);
            return 0;
		}

        skipLine(reader);
	}

	// read column: subject sequence identifier
    res = readUntilWhitespace(match.seq2, reader);
    if (res) return res;

    skipWhitespaces(reader);

	// skip column 2 - 6
	for (int i = 0; i < 5; ++i) {
        skipUntilWhitespace(reader);
        skipWhitespaces(reader);
	}

	// read column: begin position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

	// read column: end position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

	// read column: begin position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin1 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

	// read column: end position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end1 = lexicalCast<unsigned>(buffer);

    // determine orientation
	match.strand = true;
	if (match.begin1 > match.end1) {
		typename TMatch::TPos help = match.begin1;
		match.begin1 = match.end1;
		match.end1 = help;
		match.strand = !match.strand;
	}
	if (match.begin2 > match.end2) {
		typename TMatch::TPos help = match.begin2;
		match.begin2 = match.end2;
		match.end2 = help;
		match.strand = !match.strand;
	}

    skipLine(reader);
	return 0;
}

template<typename TReader, typename TMatch>
bool
readRecord(TReader & reader, TMatch & match, Blat) {
    CharString buffer;
    int res;

    match.score = -999;
    skipWhitespaces(reader);

	// read column: subject sequence identifier
    res = readUntilWhitespace(match.seq2, reader);
    if (res) return res;

    skipWhitespaces(reader);

	// read column: subject sequence identifier
    res = readUntilWhitespace(match.seq1, reader);
    if (res) return res;

    skipWhitespaces(reader);

    // skip column 2 - 5
	for (int i = 0; i < 4; ++i) {
        skipUntilWhitespace(reader);
        skipWhitespaces(reader);
	}

	// read column: begin position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

	// read column: end position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

	// read column: begin position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin1 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

	// read column: end position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end1 = lexicalCast<unsigned>(buffer);

	match.strand = true;
	if (match.begin1 > match.end1) {
		typename TMatch::TPos help = match.begin1;
		match.begin1 = match.end1;
		match.end1 = help;
		match.strand = !match.strand;
	}
	if (match.begin2 > match.end2) {
		typename TMatch::TPos help = match.begin2;
		match.begin2 = match.end2;
		match.end2 = help;
		match.strand = !match.strand;
	}

	skipLine(reader);
	return 0;
}



template<typename TReader, typename TMatch>
bool
readRecord(TReader & reader, TMatch & match, LastZ) {
    CharString buffer;
    int res;

    match.score = -999;
    skipWhitespaces(reader);

    // read column: subject sequence identifier
    res = readUntilWhitespace(match.seq1, reader);
    if (res) return res;

    skipWhitespaces(reader);

    // read column: begin position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin1 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: end position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end1 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: strand 1
    if (value(reader) == '+') match.strand = true;
    else if (value(reader) == '-') match.strand = false;
    else {
        skipLine(reader);
        return 1;
    }
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

    // read column: query sequence identifier
    res = readUntilWhitespace(match.seq2, reader);
    if (res) return res;

    skipWhitespaces(reader);

    // read column: begin position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: end position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: strand 2
    if (value(reader) == '-') match.strand = !match.strand;
    else if (value(reader) == '+') {}
    else {
        skipLine(reader);
        return 1;
    }

    skipLine(reader);
    return 0;
}

template<typename TReader, typename TMatch>
bool
readRecord(TReader & reader, TMatch & match, Last) {
    CharString buffer;
    int res;

    skipWhitespaces(reader);

    // read first column: score
    res = readDigits(buffer, reader);
    if (res) return res;
    match.score = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: subject sequence identifier
    res = readUntilWhitespace(match.seq1, reader);
    if (res) return res;

    skipWhitespaces(reader);

    // read column: begin position 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin1 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: length 1
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end1 = match.begin1 + lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: strand 1
    if (value(reader) == '+') match.strand = true;
    else if (value(reader) == '-') match.strand = false;
    else {
        skipLine(reader);
        return 1;
    }
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

    // skip one column
    skipUntilWhitespace(reader);
    skipWhitespaces(reader);

    // read column: query sequence identifier
    res = readUntilWhitespace(match.seq2, reader);
    if (res) return res;

    skipWhitespaces(reader);

    // read column: begin position 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.begin2 = lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);

    // read column: length 2
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res) return res;
    match.end2 = match.begin2 + lexicalCast<unsigned>(buffer);

    skipWhitespaces(reader);
    
    // read column: strand 2
    if (value(reader) == '-') match.strand = !match.strand;
    else if (value(reader) == '+') {}
    else {
        skipLine(reader);
        return 1;
    }

    skipLine(reader);
    return 0;
}

template<typename TMatch>
bool
readMatches(CharString & fileName,
			String<TMatch> & matches,
			Blast)
{
    std::fstream file(toCString(fileName), std::ios::in | std::ios::binary);
    if (!file.good()) {
        std::cerr << "ERROR: Could not open " << fileName << "!\n";
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > recordReader(file);

    // Read seq1 identifier from # lines
	TMatch dummyMatch;
    int res = readRecord(recordReader, dummyMatch, Blast());
    if (res) return res;

    // Read remaining matches line by line
    while (!atEnd(recordReader)) {
        if (value(recordReader) == '#') {
            skipLine(recordReader);
        }
        else {

            TMatch match;
            match.strand = true;
            match.seq1 = dummyMatch.seq1;

            res = readRecord(recordReader, match, Blast());
            if (res) return res;
        
            appendValue(matches, match);
        }
	}

    if (!atEnd(recordReader)) return 1;

    return 0;
}

template<typename TMatch, typename TTag>
bool
readMatches(CharString & fileName,
			String<TMatch> & matches,
			TTag tag)
{
    std::fstream file(toCString(fileName), std::ios::in | std::ios::binary);
    if (!file.good()) {
        std::cerr << "ERROR: Could not open " << fileName << "!\n";
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > recordReader(file);

    // Read matches line by line
    while (!atEnd(recordReader)) {
        if (value(recordReader) == '#') {
            skipLine(recordReader);
        }
        else {
            TMatch match;
            match.strand = true;

            int res = readRecord(recordReader, match, tag);
            if (res) return res;
        
            appendValue(matches, match);
        }
	}

    if (!atEnd(recordReader)) return 1;

    return 0;
}

template<typename TMatch, typename TMapping>
void
mapMatches(String<TMatch> & subjectMatches, String<TMatch> & otherMatches, TMapping & mapping)
{
	// interval tree for epsilon matches
	typedef IntervalAndCargo<unsigned, unsigned> TInterval;
	String<TInterval> intervals;
	for (unsigned i = 0; i < length(subjectMatches); ++i) {
		TMatch m = value(subjectMatches, i);
		appendValue(intervals, TInterval(m.begin1, m.end1, i));
	}
	IntervalTree<unsigned, unsigned> subjectMatchTree(intervals);

	typedef typename Iterator<String<TMatch> >::Type TIterator;
	TIterator otherIt = begin(otherMatches);
	TIterator otherEnd = end(otherMatches);

	// fill the list of overlapping matches
	while (otherIt != otherEnd) {
		// find overlaps in seq1
		String<unsigned> seq1Overlaps;
		findIntervals(subjectMatchTree, (*otherIt).begin1, (*otherIt).end1, seq1Overlaps);
		
		// check strand and overlap in seq2
		Iterator<String<unsigned> >::Type overlapIt = begin(seq1Overlaps);
		while (overlapIt != end(seq1Overlaps)) {
			if ((*otherIt).seq1 == subjectMatches[*overlapIt].seq1 && (*otherIt).seq2 == subjectMatches[*overlapIt].seq2 && 
				(*otherIt).strand == subjectMatches[*overlapIt].strand) {
				if ((*otherIt).begin2 < subjectMatches[(*overlapIt)].end2 && 
					(*otherIt).end2 > subjectMatches[(*overlapIt)].begin2) {
					// insert match into mapping
					appendValue(mapping[(*overlapIt)], position(otherIt, otherMatches));
				}
			}
			++overlapIt;
		}
		++otherIt;
	}
}

template<typename TPos, typename TMatch>
double
computeCoverage(TMatch & subjectMatch, String<TPos> & map, String<TMatch> & otherMatches) {
	// sort overlapping matches by begin position in seq1
	String<Pair<TPos> > others;
	typename Iterator<String<TPos> >::Type oIt = begin(map);
	while (oIt != end(map)) {
		TPos i = 0;
		while (i < length(others)) {
			if (otherMatches[(*oIt)].begin1 < others[i].i1) {
				break;
			}
			++i;
		}
		replace(others, i, i, Pair<TPos>(otherMatches[(*oIt)].begin1, otherMatches[(*oIt)].end1));
		++oIt;
	}

	// compute number of epsMatch positions that are not covered by other matches
	TPos uncoveredPosSeq1 = 0;
	TPos subjectPos = subjectMatch.begin1;
	typename Iterator<String<Pair<TPos> > >::Type otherIt = begin(others);
	while (otherIt != end(others)) {
		if ((*otherIt).i1 <= subjectPos) {
			subjectPos = _max(subjectPos, (*otherIt).i2);
		} else {
			uncoveredPosSeq1 += (*otherIt).i1 - subjectPos;
			subjectPos = (*otherIt).i2;
		}
		++otherIt;
	}
	if (subjectPos < subjectMatch.end1) {
		uncoveredPosSeq1 += subjectMatch.end1 - subjectPos;
	}

	// compute coverage of epsMatch
	TPos subjectLength = subjectMatch.end1 - subjectMatch.begin1;
	return double(subjectLength - uncoveredPosSeq1) * 100.0 / (double)subjectLength;
}


template<typename TMatch, typename TMapping, typename TBins>
void
binMatches(String<TMatch> & epsMatches, String<TMatch> & otherMatches, TMapping & mapping, TBins & bins) {
	// iterate over eps-matches
	typedef typename Iterator<TMapping>::Type TIterator;
	TIterator mapIt = begin(mapping);
	TIterator mapEnd = end(mapping);

	while (mapIt != mapEnd) {
		// compute coverage of epsilon-match
		TMatch epsMatch = epsMatches[position(mapIt, mapping)];
		double coverageSeq1 = computeCoverage(epsMatch, *mapIt, otherMatches);


		// increment corresponding bin counter
		if (coverageSeq1 == 0.0) {
			write(epsMatch);
			++bins[0];
		}
		else if (coverageSeq1 < 10.0) {
			//write(epsMatch);
			++bins[1];
		}
		else if (coverageSeq1 < 50.0) {
			//write(epsMatch);
			++bins[2];
		}
		else if (coverageSeq1 < 80.0) {
			//write(epsMatch);
			++bins[3];
		}
		else if (coverageSeq1 < 90.0) {
			//write(epsMatch);
			++bins[4];
		}
		else ++bins[5];
		++mapIt;
	}
}



template<typename TPos, typename TMatch>
Tuple<int, 5>
analyze(TMatch & subjectMatch, String<TPos> & map, String<TMatch> & otherMatches, bool verbose) \
{
    Tuple<int, 5> result;
    int & exactHits = result.i[0];
    int & scoreDist = result.i[1];
    int & sameStart = result.i[2];
    int & sameEnd   = result.i[3];
    int & overlap   = result.i[4];

    exactHits = 0;
    scoreDist = 0;
    sameStart = 0;
    sameEnd   = 0;
    overlap   = 0;


	typename Iterator<String<TPos> >::Type oIt = begin(map);
	while (oIt != end(map))
    {
		if (distance(subjectMatch, otherMatches[*oIt]) == 0)
        {
            if (scoreDiff(subjectMatch, otherMatches[*oIt]) == 0)
                ++exactHits;
            else {
                if (verbose) std::cout << "Score difference: " << std::endl;
                if (verbose) std::cout << "ref:  "; write(subjectMatch);
                if (verbose) std::cout << "comp: "; write(otherMatches[*oIt]);
                ++scoreDist;
            }
        }
        else
        {
            if ( subjectMatch.begin1 == otherMatches[*oIt].begin1 && subjectMatch.begin2 == otherMatches[*oIt].begin2)
            {
                if (verbose) std::cout << "Same start: " << std::endl;
                if (verbose) std::cout << "ref:  "; write(subjectMatch);
                if (verbose) std::cout << "comp: "; write(otherMatches[*oIt]);
                ++sameStart;
            }
            else if ( subjectMatch.end1 == otherMatches[*oIt].end1 && subjectMatch.end2 == otherMatches[*oIt].end2)
            {
                if (verbose) std::cout << "Same end: " << std::endl;
                if (verbose) std::cout << "ref:  "; write(subjectMatch);
                if (verbose) std::cout << "comp: "; write(otherMatches[*oIt]);
                ++sameEnd;
            }
        }
        ++oIt;
	}
    return result;
}


template<typename TMatch, typename TMapping>
void
exactAnalysis(String<TMatch> & epsMatches, String<TMatch> & otherMatches, TMapping & mapping, bool verbose = false) {
	// iterate over eps-matches
	typedef typename Iterator<TMapping>::Type TIterator;
	TIterator mapIt = begin(mapping);
	TIterator mapEnd = end(mapping);

    Tuple<int, 5> r;
    std::fill(&(r.i[0]), &(r.i[0])+5, 0);
	while (mapIt != mapEnd)
    {
		// compute coverage of epsilon-match
		TMatch epsMatch = epsMatches[position(mapIt, mapping)];
		Tuple<int, 5> p = analyze(epsMatch, *mapIt, otherMatches, verbose);
        r.i[0] += p.i[0];
        r.i[1] += p.i[1];
        r.i[2] += p.i[2];
        r.i[3] += p.i[3];
        r.i[4] += p.i[4];
        ++mapIt;
	}
    std::cout << "exactHits: " << r.i[0] << std::endl;
    std::cout << "scoreDiff: " << r.i[1] << std::endl;
    std::cout << "sameStart: " << r.i[2] << std::endl;
    std::cout << "sameEnd:   " << r.i[3] << std::endl;
    std::cout << "overlap:   " << r.i[4] << std::endl;

}



void
_setupParser(ArgumentParser & parser) {

    setShortDescription(parser, "comparison of local alignments");
    setVersion(parser, "0.2");
    setDate(parser, "Jul 2013");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIREFMATCHFILE\\fP \\fICOMPMATCHFILE\\fP");

    addDescription(parser,
                   "Computes the coverage of matches in reference file by matches from comparison file and "
                   "counts the number of matches with a coverage <10%, 10-50%, 50-80%, 80-90%, and >90%.");
    addDescription(parser, "(c) 2013 by Birte Kehr");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "REFMATCHFILE")); // file containing reference matches , e.g. all eps-matches
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "COMPMATCHFILE")); // file containing comparison matches

    addSection(parser, "Main Options");
    addOption(parser, ArgParseOption("fr", "file_format1", "format of file containing reference matches", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("fc", "file_format2", "format of file containing comparison matches", ArgParseArgument::STRING));
    setValidValues(parser, "fr", "gff last lastz blat blast");
    setValidValues(parser, "fc", "gff last lastz blat blast");

    addOption(parser, ArgParseOption("e", "extra", "Count how many matches are exactly equal"));
}

int main(int argc, const char *argv[]) {
	typedef _Match<unsigned> TMatch;

    CharString refFile, compFile;
    CharString refFormat, compFormat;
    String<TMatch> refMatches, compMatches;
    bool       extraAnalysis = false;

    ArgumentParser parser("compare_matches");
    _setupParser(parser);
    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
    if (res == ArgumentParser::PARSE_OK)
    {
        getArgumentValue(refFile, parser, 0);
        getArgumentValue(compFile, parser, 1);
        getOptionValue(refFormat, parser, "fr");
        getOptionValue(compFormat, parser, "fc");
        if (isSet(parser, "extra"))
            extraAnalysis = true;
    }

    if (res != ArgumentParser::PARSE_OK) {
        return res == ArgumentParser::PARSE_ERROR; // 1 - error, 0 - otherwise
    }

    //std::cout << "Reference file: " << refFile << std::endl;
    //std::cout << "Compared file : " << compFile << std::endl;

    // read subject-matches file
    if (refFormat == "gff")
    {
        if (readMatches(refFile, refMatches, Gff())) return false;
    }
    else if (refFormat == "lastz")
    {
        if (readMatches(refFile, refMatches, LastZ())) return false;
    }
    else if (refFormat == "last")
    {
        if (readMatches(refFile, refMatches, Last())) return false;
    }
    else if (refFormat == "blat")
    {
        if (readMatches(refFile, refMatches, Blat())) return false;
    }
    else if (refFormat == "blast")
    {
        if (readMatches(refFile, refMatches, Blast())) return false;
    }
    else
    {
        std::cerr << "Unknown file format: " << refFormat << std::endl;
        return 1;
    }

    // read comparison match file
    if (compFormat == "gff")
    {
        if(readMatches(compFile, compMatches, Gff())) return false;
    }
    else if (compFormat == "lastz")
    {
        if (readMatches(compFile, compMatches, LastZ())) return false;
    }
    else if (compFormat == "last")
    {
        if (readMatches(compFile, compMatches, Last())) return false;
    }
    else if (compFormat == "blat")
    {
        if (readMatches(compFile, compMatches, Blat())) return false;
    }
    else if (compFormat == "blast")
    {
        if (readMatches(compFile, compMatches, Blast())) return false;
    }
    else
    {
        std::cerr << "Unknown file format: " << compFormat << std::endl;
        return 1;
    }

    std::cout << "# subject matches: " << length(refMatches) << std::endl;
    std::cout << "# other matches: " << length(compMatches) << std::endl;

    // for all epsMatches: list of overlapping matches (indices from otherMatches)
    typedef String<String<unsigned> > TMapping;
    TMapping mapping;
    resize(mapping, length(refMatches));

    mapMatches(refMatches, compMatches, mapping);

    // bins that count the number of matches with a certain coverage
    typedef Size<String<TMatch> >::Type TSize;
    String<TSize> bins;
    resize(bins, 6, 0);

    // compute coverage and bin matches
    binMatches(refMatches, compMatches, mapping, bins);

    // output bins
    for (TSize i = 0; i < length(bins); ++i) {
        std::cout << "bin " << i << ": " << value(bins, i) << " matches" << std::endl;
    }

    // new analysis
    if (extraAnalysis)
        exactAnalysis(refMatches, compMatches, mapping);

    return 0;
}
