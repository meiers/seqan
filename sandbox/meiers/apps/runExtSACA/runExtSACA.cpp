// ==========================================================================
//                                  runSACA
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <iostream>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

namespace seqan {
template<> struct SAValue<String<Dna5, MMap<> > >
{
    typedef unsigned Type;
};

template<> struct SAValue<StringSet<String<Dna5, MMap<> >, Owner<ConcatDirect<> > > >
{
    typedef Pair<unsigned, unsigned, Pack> Type;
};

template <typename TIndexSpec> 
struct Fibre< Index< StringSet<String<Dna5, MMap<> >, Owner<ConcatDirect<> > > const, TIndexSpec>, FibreSA> 
{
    typedef StringSet<String<Dna5, MMap<> >, Owner<ConcatDirect<> > > TSet;
    typedef String<typename SAValue<TSet>::Type, MMap<> > Type;
};

    template <typename TIndexSpec>
    struct Fibre< Index< String<Dna5, MMap<> > const, TIndexSpec>, FibreSA>
    {
        typedef String<typename SAValue<String<Dna5, MMap<> > >::Type, MMap<> > Type;
    };
}

// --------------------------------------------------------------------------
// Class RunExtSACAOptions
// --------------------------------------------------------------------------

struct RunExtSACAOptions
{
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<> >, 1> > TShape_1;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1> >, 1> > TShape_2;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,7> >, 0> > TShape_3;
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,3,1,2> >, 1> > TShape_4;
    typedef CyclicShape<FixedShape<0, ShapePatternHunter, 0> > TShape_5;
    typedef CyclicShape<FixedShape<0, ShapePatternHunter, 2> > TShape_6;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<2,1,1,3,2,2,1,1,2,1,2,2,1,1,2,2,2,1> >, 0> > TShape_7;

    CharString file;
    CharString outfile;
    int shape;
    bool concatenation;
    
    void print()
    {
        std::cout << "RunSACA Options" << std::endl;
        std::cout << "    file:      " << file << std::endl;
        std::cout << "    shape:     ";
        CharString sh;
        switch (shape) {
            case 0: std::cout << "no shape" << std::endl; break;
            case 1: cyclicShapeToString(sh, TShape_1()); std::cout << sh << std::endl; break;
            case 2: cyclicShapeToString(sh, TShape_2()); std::cout << sh << std::endl; break;
            case 3: cyclicShapeToString(sh, TShape_3()); std::cout << sh << std::endl; break;
            case 4: cyclicShapeToString(sh, TShape_4()); std::cout << sh << std::endl; break;
            case 5: cyclicShapeToString(sh, TShape_5()); std::cout << sh << std::endl; break;
            case 6: cyclicShapeToString(sh, TShape_6()); std::cout << sh << std::endl; break;
            case 7: cyclicShapeToString(sh, TShape_7()); std::cout << sh << std::endl; break;
        }
    }
};


// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(RunExtSACAOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("runSACA");
    setShortDescription(parser, "Load a file and run an external SACA on it");
    setVersion(parser, "0.1");
    setDate(parser, "April 2014");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIFASTA FILE\\fP\" \"\\fIOUTFILE\\fP\"   ");
    addDescription(parser, "0:        no shape");
    addDescription(parser, "1:        10 (1/2)");
    addDescription(parser, "2:        110 (2/3)");
    addDescription(parser, "3:        110000001 (3/9)");
    addDescription(parser, "4:        011001101 (5/9)");
    addDescription(parser, "5:        111010010100110111 (11/18, Pattern Hunter)");
    addDescription(parser, "6:        11101001010011011100 (11/20, Pattern Hunter 2)");
    addDescription(parser, "7:        101110010101110110101110101011 (19/30)");
    
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE, "OUTPUTFILE"));
    addOption(parser, seqan::ArgParseOption("c", "concat", "Force text to be read as a single string by concatenating entries. If not set a string set will be used"));
    addOption(parser, seqan::ArgParseOption("s", "shape", "Specify the cyclic shape.", ArgParseArgument::INTEGER, "[0..7]"));
    setMinValue(parser, "shape", "0");
    setMaxValue(parser, "shape", "7");
    setDefaultValue(parser, "shape", "0");
    
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    
    // Extract option values.
    getOptionValue(options.shape, parser, "shape");
    getArgumentValue(options.file, parser, 0);
    getArgumentValue(options.outfile, parser, 1);
    options.concatenation = isSet(parser, "concat");

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function build_Index() => translate runtime to compile time arguments
// --------------------------------------------------------------------------

template <typename TText, typename TShape>
void build_Index(TText const & text, TShape const &, RunExtSACAOptions const & options)
{
    Index<TText const, IndexSa<Gapped<ModCyclicShape<TShape> > > > index(text);
    std::cout << "Open " << options.outfile << " for the SA" << std::endl;
    
    if (!open(indexSA(index), toCString(options.outfile), OPEN_CREATE | OPEN_RDWR)) std::cout << "cannot open " << options.outfile << " to store SA" << std::endl;
    
    indexCreate(index, FibreSA(), DislexExternal<TShape>());
    std::cout << "done creating the gapped Index" << std::endl;
}

template <typename TText>
void build_Index_ungapped(TText const & text, RunExtSACAOptions const & options)
{
    Index<TText const, IndexSa<> > index(text);
    if (!open(indexSA(index), toCString(options.outfile), OPEN_CREATE | OPEN_RDWR)) std::cout << "cannot open " << options.outfile << " to store SA" << std::endl;

    resize(indexSA(index), length(indexRawText(index)), Exact());
    _createSuffixArrayPipelining(indexSA(index), indexText(index), Skew7()); // NOTE: ONLY SKEW !!!! 
    std::cout << "done creating an ungapped Index" << std::endl;
    

}

// --------------------------------------------------------------------------
// Function call() => translate runtime to compile time arguments
// --------------------------------------------------------------------------



template <typename T>
void call(T const & text, RunExtSACAOptions const & options)
{
    switch (options.shape) {
        case 0: build_Index_ungapped(text, options); break;
        case 1: build_Index(text, RunExtSACAOptions::TShape_1(), options); break;
        case 2: build_Index(text, RunExtSACAOptions::TShape_2(), options); break;
        case 3: build_Index(text, RunExtSACAOptions::TShape_3(), options); break;
        case 4: build_Index(text, RunExtSACAOptions::TShape_4(), options); break;
        case 5: build_Index(text, RunExtSACAOptions::TShape_5(), options); break;
        case 6: build_Index(text, RunExtSACAOptions::TShape_6(), options); break;
        case 7: build_Index(text, RunExtSACAOptions::TShape_7(), options); break;  
   }
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    RunExtSACAOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    options.print();
    
    
    // Open memory mapped string.
    seqan::String<char, seqan::MMap<> > mmapString;
    if (!open(mmapString, toCString(options.file), seqan::OPEN_RDONLY))
        return 1;

    // Create RecordReader.
    seqan::RecordReader<seqan::String<char, seqan::MMap<> >,
                        seqan::DoublePass<seqan::StringReader> > reader(mmapString);

    // Read file in one pass.
    StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > > ids;
    StringSet<String<Dna5, MMap<> >, Owner<ConcatDirect<> > > seqs;

    double teim=sysTime();

    if (read2(ids, seqs, reader, Fasta()) != 0) {
        std::cout << "cannot read file" << std::endl;
        return 1;
    }
    clear(ids);
    std::cout << "Time for reading MMap String: " << sysTime() - teim << " (" << length(seqs) << "x, " << lengthSum(seqs) << ")" << std::endl;


    if (!options.concatenation)
    {
        std::cout <<"String Set" << std::endl;
        call(seqs, options);
    }
    else {
        std::cout <<"String" << std::endl;
        call(seqs.concat, options);
    }


    return 0;
}
