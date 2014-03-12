// ==========================================================================
//                      benchmarkModCyclicShapeIterator
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
#include <seqan/modifier.h>

#include <time.h>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

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
    ArgumentParser parser("seqanBuildSuffixArray");
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
    //addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    //addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

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
// Function benchmarkIterator()
// --------------------------------------------------------------------------

template <typename TMessage, typename TModString>
void benchmarkIterator(TModString const & txt, TMessage const & label)
{
    typedef typename Iterator<TModString const, Standard>::Type TModIter;
    typename Size<TModString>::Type count = 0;
    TModIter begIt = begin(txt, Standard());
    TModIter endIt = end(txt, Standard());
    TModIter it;

    double start = sysTime();
    for(unsigned I = 0; I < 5; ++I) {
        // forward:
        for(it = begIt; it != endIt; ++it) if(*it == 'A') ++count;
        // backward:
        for(it = endIt; it != begIt; --it) if(*(it-1) == 'A') ++count;
    }
    std::cout << length(txt) << "\t" << count << "\t" << sysTime() - start << "s" << "\t" << label << std::endl;
}


// --------------------------------------------------------------------------
// Function callBenchmarks()
// --------------------------------------------------------------------------

template <typename TString>
void callBenchmarks(TString const & txt)
{
    typedef CyclicShape<GenericShape> ShapeGen;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<> >, 0> > Shape1;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,1> >, 0> > Shape3;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,1,1> >, 0> > Shape4;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,1,1,1,1,1,1,1,1,1,1,1,1,1,1> >, 0> > Shape16;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1> >, 0> > Shape17;
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,2,3,1,1,3,3> >, 0> > ShapeGap16; // 0110100111001001
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,2,3,1,1,3,3> >, 1> > ShapeGap17; // 01101001110010010

    std::cout << "|T|    \t#As    \ttime  \ttype" << std::endl;


    benchmarkIterator(txt, "noShape");


    ShapeGen s1; stringToCyclicShape(s1, "1");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s1),
                      "Generic:1");

    Shape1 s2;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<Shape1> >(txt, s2),
                      "Fixed:1");


    ShapeGen s3; stringToCyclicShape(s3, "111");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s3),
                      "Generic:111");

    Shape3 s4;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<Shape3> >(txt, s4),
                      "Fixed:111");


    ShapeGen s5; stringToCyclicShape(s5, "1111");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s5),
                      "Generic:1111");

    Shape4 s6;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<Shape4> >(txt, s6),
                      "Fixed:1111");


    ShapeGen s7; stringToCyclicShape(s7, "1111111111111111");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s7),
                      "Generic:16/16");

    Shape16 s8;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<Shape16> >(txt, s8),
                      "Fixed:16/16");


    ShapeGen s9; stringToCyclicShape(s7, "11111111111111111");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s9),
                      "Generic:17/17");

    Shape17 s10;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<Shape17> >(txt, s10),
                      "Fixed:17/17");



    ShapeGen s11; stringToCyclicShape(s11, "0110100111001001");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s11),
                      "Generic:8/16");

    ShapeGap16 s12;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGap16> >(txt, s12),
                      "Fixed:8/16");


    ShapeGen s13; stringToCyclicShape(s13, "01101001110010010");
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGen> >(txt, s13),
                      "Generic:8/17");

    ShapeGap17 s14;
    benchmarkIterator(ModifiedString<TString const, ModCyclicShape<ShapeGap17> >(txt, s14),
                      "Fixed:8/17");


}





// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // options.infile = "/Users/mauli87/Downloads/hs_alt_CHM1_1.0_chr22.fa";


    std::fstream in( seqan::toCString(options.infile), std::ios::binary | std::ios::in);
    //seqan::RecordReader<std::fstream, seqan::DoublePass<seqan::StringReader> > reader(in);
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

    seqan::StringSet<seqan::CharString> ids;
    typedef seqan::StringSet<seqan::Dna5String, Owner<ConcatDirect<> > > TSet;
    TSet seqs;

    if (!in.good())
    {
        std::cerr << "Couldn't open " << options.infile << std::endl;
        return 1;
    }
    if (read2(ids, seqs, reader, seqan::Fasta()) != 0)
        return 2;  // Could not record from file.


    time_t t; time(&t);
    std::cout << "# input file : " << options.infile << std::endl;
    std::cout << "# timestamp  : " << ctime(&t) << std::endl;

    callBenchmarks(concat(seqs));
    
    return 0;
}
