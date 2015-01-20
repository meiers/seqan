// ==========================================================================
//                                  SAK
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Swiss Army Knife tool... "It slices, it dices and it makes the laundry!"
//
// Rewrite of the original sak tool.
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// Class SakOptions
// --------------------------------------------------------------------------

struct SakOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to FASTA/FASTQ file.
    seqan::CharString inFastxPath;

    // Path to output file.
    seqan::CharString outPath;

    // Set if one sequence is to be retrieved.
    seqan::String<__uint64> seqIndices;

    // Set if multiple sequences are to be retrieved.
    seqan::String<seqan::Pair<__uint64> > seqIndexRanges;

    // Set if output is to be limited to an infix.
    __uint64 seqInfixBegin;
    __uint64 seqInfixEnd;

    // Whether or not to reverse-complement the result.
    bool reverseComplement;

    // Maximal length of sequence characters to print.
    __uint64 maxLength;

    // Prefix of read names to output if not empty.
    seqan::CharString readPattern;

    // Line length configuration etc.
    seqan::SequenceOutputOptions seqOutOptions;

    SakOptions() :
        verbosity(1),
        seqInfixBegin(seqan::maxValue<__uint64>()),
        seqInfixEnd(seqan::maxValue<__uint64>()),
        reverseComplement(false),
        maxLength(seqan::maxValue<__uint64>())
    {
    }
};

// --------------------------------------------------------------------------
// Function parseRange()
// --------------------------------------------------------------------------

template <typename TNum>
bool parseRange(TNum & beginPos, TNum & endPos, seqan::CharString const & rangeStr)
{
    seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type reader = directionIterator(rangeStr, seqan::Input());

    // Parse out begin position.
    seqan::CharString buffer;
    readUntil(buffer, reader, seqan::EqualsChar<'-'>(), seqan::EqualsChar<','>());
    if (!lexicalCast(beginPos, buffer))
        return false;

    if (atEnd(reader))
        return false;

    skipOne(reader);    // Skip '-'.

    // Parse out end position.
    clear(buffer);
    readUntil(buffer, reader, seqan::False(), seqan::EqualsChar<','>());

    if (empty(buffer))
    {
        endPos = seqan::maxValue<TNum>();
        return true;
    }

    if (!lexicalCast(endPos, buffer))
        return false;

    return (beginPos <= endPos);
}

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(SakOptions & options,
          int argc,
          char ** argv)
{
    seqan::ArgumentParser parser("sak");
    setShortDescription(parser, "Slicing and dicing of FASTA/FASTQ files..");
    setVersion(parser, "0.2");
    setDate(parser, "November 2012");
    setCategory(parser, "Utilities");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fB-o\\fP \\fIOUT.{fa,fq}\\fP] \\fIIN.{fa,fq}\\fP");
    addDescription(parser, "\"It slices, it dices and it makes the laundry!\"");
    addDescription(parser, "Original SAK tool by David Weese. Rewrite by Manuel Holtgrewe.");

    // The only argument is the input file.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "IN"));

    // Only FASTA and FASTQ files are allowed as input.
    setValidValues(parser, 0, seqan::SeqFileIn::getFileExtensions());

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-path",
                                            "Path to the resulting file.  If omitted, result is printed to stdout in FastQ format.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "FASTX"));
    setValidValues(parser, "out-path", seqan::SeqFileOut::getFileExtensions());
    addOption(parser, seqan::ArgParseOption("rc", "revcomp", "Reverse-complement output."));
    addOption(parser, seqan::ArgParseOption("l", "max-length", "Maximal number of sequence characters to write out.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));

    addSection(parser, "Filter Options");
    addOption(parser, seqan::ArgParseOption("s", "sequence", "Select the given sequence for extraction by 0-based index.",
                                            seqan::ArgParseOption::INTEGER, "NUM", true));
    addOption(parser, seqan::ArgParseOption("sn", "sequence-name", "Select sequence with name prefix being \\fINAME\\fP.",
                                            seqan::ArgParseOption::STRING, "NAME", true));
    addOption(parser, seqan::ArgParseOption("ss", "sequences",
                                            "Select sequences \\fIfrom\\fP-\\fIto\\fP where \\fIfrom\\fP and \\fIto\\fP "
                                            "are 0-based indices.",
                                            seqan::ArgParseArgument::STRING, "RANGE", true));
    addOption(parser, seqan::ArgParseOption("i", "infix",
                                            "Select characters \\fIfrom\\fP-\\fIto\\fP where \\fIfrom\\fP and \\fIto\\fP "
                                            "are 0-based indices.",
                                            seqan::ArgParseArgument::STRING, "RANGE", true));

    addOption(parser, seqan::ArgParseOption("ll", "line-length",
                                            "Set line length in output file.  See section \\fILine Length\\fP for details.",
                                            seqan::ArgParseArgument::INTEGER, "LEN", false));
    setMinValue(parser, "line-length", "-1");

    addTextSection(parser, "Line Length");
    addText(parser,
            "You can use the setting \\fB--line-length\\fP for setting the resulting line length.  By default, "
            "sequences in FASTA files are written with at most 70 characters per line and sequences in FASTQ files are "
            "written without any line breaks.  The quality sequence in FASTQ file is written in the same way as the "
            "residue sequence.");
    addText(parser,
            "The default is selected with a \\fB--line-length\\fP value of \\fI-1\\fP and line breaks can be disabled "
            "with a value of \\fI0\\fP.");

    addTextSection(parser, "Usage Examples");
    addListItem(parser, "\\fBsak\\fP \\fB-s\\fP \\fI10\\fP \\fIIN.fa\\fP",
                "Cut out 11th sequence from \\fIIN.fa\\fP and write to stdout as FASTA.");
    addListItem(parser, "\\fBsak\\fP \\fB-ss\\fP \\fI10-12\\fP \\fB-ss\\fP \\fI100-200\\fP \\fIIN.fq\\fP",
                "Cut out 11th up to and including 12th and 101th up to and including 199th sequence from \\fIIN.fq\\fP "
                "and write to stdout as FASTA.");

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.inFastxPath, parser, 0);

    seqan::CharString tmp;
    getOptionValue(tmp, parser, "out-path");

    if (isSet(parser, "out-path"))
        getOptionValue(options.outPath, parser, "out-path");

    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    if (isSet(parser, "sequence"))
    {
        std::vector<std::string> sequenceIds = getOptionValues(parser, "sequence");
        for (unsigned i = 0; i < seqan::length(sequenceIds); ++i)
        {
            unsigned idx = 0;
            if (!seqan::lexicalCast(idx, sequenceIds[i]))
            {
                std::cerr << "ERROR: Invalid sequence index " << sequenceIds[i] << "\n";
                return seqan::ArgumentParser::PARSE_ERROR;
            }
            appendValue(options.seqIndices, idx);
        }
    }

    if (isSet(parser, "sequences"))
    {
        std::vector<std::string> sequenceRanges = getOptionValues(parser, "sequences");
        seqan::CharString buffer;
        for (unsigned i = 0; i < seqan::length(sequenceRanges); ++i)
        {
            seqan::Pair<__uint64> range;
            if (!parseRange(range.i1, range.i2, sequenceRanges[i]))
            {
                std::cerr << "ERROR: Invalid range " << sequenceRanges[i] << "\n";
                return seqan::ArgumentParser::PARSE_ERROR;
            }
            appendValue(options.seqIndexRanges, range);
        }
    }

    if (isSet(parser, "infix"))
    {
        seqan::CharString buffer;
        getOptionValue(buffer, parser, "infix");
        if (!parseRange(options.seqInfixBegin, options.seqInfixEnd, buffer))
        {
            std::cerr << "ERROR: Invalid range " << buffer << "\n";
            return seqan::ArgumentParser::PARSE_ERROR;
        }
    }

    options.reverseComplement = isSet(parser, "revcomp");

    if (isSet(parser, "max-length"))
        getOptionValue(options.maxLength, parser, "max-length");

    if (isSet(parser, "sequence-name"))
        getOptionValue(options.readPattern, parser, "sequence-name");

    getOptionValue(options.seqOutOptions.lineLength, parser, "line-length");

    return res;
}

// ---------------------------------------------------------------------------
// Function yesNo()
// ---------------------------------------------------------------------------

char const * yesNo(bool b)
{
    if (b)
        return "YES";
    else
        return "NO";
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    double startTime = 0;

    // Parse command line.
    SakOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 2)
    {
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    " << options.verbosity << "\n"
                  << "IN           " << options.inFastxPath << "\n"
                  << "OUT          " << options.outPath << "\n"
                  << "INFIX BEGIN  " << options.seqInfixBegin << "\n"
                  << "INFIX END    " << options.seqInfixEnd << "\n"
                  << "MAX LEN      " << options.maxLength << "\n"
                  << "READ PATTERN " << options.readPattern << "\n"
                  << "REVCOMP      " << yesNo(options.reverseComplement) << "\n"
                  << "SEQUENCES\n";
        for (unsigned i = 0; i < length(options.seqIndices); ++i)
            std::cerr << "  SEQ  " << options.seqIndices[i] << "\n";
        for (unsigned i = 0; i < length(options.seqIndexRanges); ++i)
            std::cerr << "  SEQS " << options.seqIndexRanges[i].i1 << "-" << options.seqIndexRanges[i].i2 << "\n";
    }

    // -----------------------------------------------------------------------
    // Open Files.
    // -----------------------------------------------------------------------
    seqan::SeqFileIn inFile;
    seqan::SeqFileOut outFile;

    bool openRes = false;
    if (!empty(options.inFastxPath))
        openRes = open(inFile, toCString(options.inFastxPath));
    else
        openRes = open(inFile, std::cin);
    if (!openRes)
    {
        std::cerr << "ERROR: Problem opening input file.\n";
        return 1;
    }

    if (!empty(options.outPath))
        openRes = open(outFile, toCString(options.outPath));
    else
        openRes = open(outFile, std::cout, seqan::Fastq());
    if (!openRes)
    {
        std::cerr << "ERROR: Problem opening output file.\n";
        return 1;
    }

    // Compute index of last sequence to write if any.
    __uint64 endIdx = seqan::maxValue<__uint64>();
    for (unsigned i = 0; i < length(options.seqIndices); ++i)
        if (endIdx == seqan::maxValue<__uint64>() || endIdx > options.seqIndices[i] + 1)
            endIdx = options.seqIndices[i] + 1;
    for (unsigned i = 0; i < length(options.seqIndexRanges); ++i)
        if (endIdx == seqan::maxValue<__uint64>() || endIdx > options.seqIndexRanges[i].i2)
            endIdx = options.seqIndexRanges[i].i2;
    if (options.verbosity >= 2)
        std::cerr << "Sequence end idx: " << endIdx << "\n";

    // -----------------------------------------------------------------------
    // Read and Write Filtered.
    // -----------------------------------------------------------------------
    startTime = seqan::sysTime();

    unsigned idx = 0;
    __uint64 charsWritten = 0;
    seqan::CharString id;
    seqan::CharString seq;
    seqan::CharString quals;
    while (!atEnd(inFile) && charsWritten < options.maxLength && idx < endIdx)
    {
        try
        {
            readRecord(id, seq, quals, inFile);
        }
        catch (seqan::ParseError const & e)
        {
            std::cerr << "ERROR: Problem reading file: " << e.what() << "\n";
            return 1;
        }

        // Check whether to write out sequence.
        bool writeOut = false;
        if (empty(options.seqIndices) && empty(options.seqIndexRanges))
            writeOut = true;
        // One of options.seqIndices.
        if (!writeOut)
        {
            for (unsigned i = 0; i < length(options.seqIndices); ++i)
            {
                if (options.seqIndices[i] == idx)
                {
                    writeOut = true;
                    break;
                }
            }
        }
        // One of options.seqIndexRanges.
        if (!writeOut)
        {
            for (unsigned i = 0; i < length(options.seqIndexRanges); ++i)
            {
                if (idx >= options.seqIndexRanges[i].i1 && idx < options.seqIndexRanges[i].i2)
                {
                    writeOut = true;
                    break;
                }
            }
        }
        // Name pattern matches.
        if (!writeOut && !empty(options.readPattern))
        {
            unsigned l = length(options.readPattern);
            if (l > length(id))
                l = length(id);
            if (prefix(id, l) == prefix(options.readPattern, l))
                writeOut = true;
        }

        // Write out if we want this.
        if (writeOut)
        {
            // Get begin and end index of infix to write out.
            __uint64 infixBegin = 0;
            if (options.seqInfixBegin != seqan::maxValue<__uint64>())
                infixBegin = options.seqInfixBegin;
            if (infixBegin > length(seq))
                infixBegin = length(seq);
            __uint64 infixEnd = length(seq);
            if (options.seqInfixEnd < length(seq))
                infixEnd = options.seqInfixEnd;
            if (infixEnd < infixBegin)
                infixEnd = infixBegin;
            if (options.verbosity >= 3)
                std::cerr << "INFIX\tbegin:" << infixBegin << "\tend:" << infixEnd << "\n";

            if (options.reverseComplement)
            {
                seqan::Dna5String seqCopy = seq;
                reverseComplement(seqCopy);
                reverse(quals);
                infixEnd = length(seq) - infixEnd;
                infixBegin = length(seq) - infixBegin;
                std::swap(infixEnd, infixBegin);

                writeRecord(outFile, id, infix(seqCopy, infixBegin, infixEnd), infix(quals, infixBegin, infixEnd));
            }
            else
                writeRecord(outFile, id, infix(seq, infixBegin, infixEnd), infix(quals, infixBegin, infixEnd));
        }

        // Advance counter idx.
        idx += 1;
    }

    if (options.verbosity >= 2)
        std::cerr << "Took " << (seqan::sysTime() - startTime) << " s\n";

    return 0;
}
