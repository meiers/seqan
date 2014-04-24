// ==========================================================================
//                                 uselessApp
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


#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace seqan;

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

struct AppOptions
{
    seqan::String<char> infile;
    seqan::String<char> outfile;
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
    setShortDescription(parser, "Builds an ESA-index.");
    setVersion(parser, "Sascha.0.1");
    setDate(parser, "July 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "Builds an index (ESA) of the specified genome and writes it to disk.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "IN"));

    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "Output file (<IN>.index if not specified)",
                                            seqan::ArgParseArgument::OUTPUTFILE, "OUT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // set outfile to infile.db if not outfile was specified
    seqan::String<char> outfile;
    seqan::getArgumentValue(options.infile, parser, 0);

    outfile = options.	infile;
    outfile += ".index";
    getOptionValue(outfile, parser, "output");
    options.outfile = outfile;

    return seqan::ArgumentParser::PARSE_OK;
}






// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
using namespace seqan;
int main(int argc, char const ** argv)
{
    Dna5String ref;
    Dna5String comp;
    CharString id;

    SequenceStream seqStream("/Users/mauli87/Development/seqan-build/chr1.fa");
    if (isGood(seqStream) && !readRecord(id, ref, seqStream))
        std::cout << "read reference file: ok" << std::endl;

    clear(id);
    SequenceStream seqStream2("/Users/mauli87/Development/seqan-build/chr20_head10k.fa");
    if (isGood(seqStream2) && !readRecord(id, comp, seqStream2))
        std::cout << "read comparison file: ok" << std::endl;

    clear(id);

    if(argc != 5) {
        std::cout << "no positions given. Use R_from R_to C_from C_to" << std::endl;
        return 1;
    }

    unsigned rF = atoi(argv[1]);
    unsigned rT = atoi(argv[2]);
    unsigned cF = atoi(argv[3]);
    unsigned cT = atoi(argv[4]);

    std::stringstream ss;
    ss << rF;
    CharString s = "*";
    append(s,ss.str());

    std::cout << s;
    for(unsigned i=rF+length(s); i<rT; ++i)
    {
        if(i%100 == 0) std::cout << "|";
        else if(i%50 == 0) std::cout << ":";
        else if(i%10 ==0) std::cout << ".";
        else std::cout << " ";
    }
    std::cout << std::endl;

    //ref sequence
    std::cout << infix(ref,  rF, rT) << std::endl;


    // alignment symbols
    for(unsigned i=0; i < std::min(rT-rF, cT-cF); ++i)
        if(ref[rF+i] == comp[cF+i]) std::cout << "+";
        else std::cout << " ";
    std::cout << std::endl;


    // compare seq
    std::cout << infix(comp, cF, cT) << std::endl;


    return 0;
}

