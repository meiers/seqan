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
#include <seqan/arg_parse.h>
#include <seqan/index.h>
using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class RunSACAOptions
// --------------------------------------------------------------------------

struct RunSACAOptions
{
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<> >, 1> > TShape_1;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1> >, 1> > TShape_2;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,7> >, 0> > TShape_3;
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,3,1,2> >, 1> > TShape_4;
    typedef CyclicShape<FixedShape<0,ShapePatternHunter, 0> > TShape_5;
    typedef CyclicShape<FixedShape<0,ShapePatternHunter, 2> > TShape_6;
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<2,1,1,3,2,2,1,1,2,1,2,2,1,1,2,2,2,1> >, 0> > TShape_7;

    CharString file;
    CharString algorithm;
    int shape;
    bool set;
    void print()
    {
        std::cout << "RunSACA Options" << std::endl;
        std::cout << "    file:      " << file << std::endl;
        std::cout << "    algorithm: " << algorithm << std::endl;
        std::cout << "    set?       " << set << std::endl;
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

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(RunSACAOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("runSACA");
    setShortDescription(parser, "Load a file and run a SACA on it");
    setVersion(parser, "0.1");
    setDate(parser, "April 2014");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIFASTA FILE\\fP\"");
    addDescription(parser, "Pass a fasta file that should be indexed and specify "
                "a SACA and a cyclic shape using the options -a and -s. The "
                "following shapes are available:");
    addDescription(parser, "0:        no shape");
    addDescription(parser, "1:        10 (1/2)");
    addDescription(parser, "2:        110 (2/3)");
    addDescription(parser, "3:        110000001 (3/9)");
    addDescription(parser, "4:        011001101 (5/9)");
    addDescription(parser, "5:        111010010100110111 (11/18, Pattern Hunter)");
    addDescription(parser, "6:        11101001010011011100 (11/20, Pattern Hunter 2)");
    addDescription(parser, "7:        101110010101110110101110101011 (19/30)");
    
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE"));

    addOption(parser, seqan::ArgParseOption("a", "algorithm", "Specify the SACA algorithm to use.", ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "algorithm", "DislexSkew7 DislexExternal InplaceRadixSort QSort Skew7 SkewExternal None");
    setDefaultValue(parser, "algorithm", "None");
    addOption(parser, seqan::ArgParseOption("s", "shape", "Specify the cyclic shape.", ArgParseArgument::INTEGER, "[0..7]"));
    setMinValue(parser, "shape", "0");
    setMaxValue(parser, "shape", "7");
    setDefaultValue(parser, "shape", "0");
    addOption(parser, seqan::ArgParseOption("S", "set", "Force text to be read as a string set. If this is not set, it will be interpret as a single string!"));
    
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    // Extract option values.
    if (isSet(parser, "algorithm"))
        {}
    getOptionValue(options.shape, parser, "shape");
    getOptionValue(options.algorithm, parser, "algorithm");
    options.set = isSet(parser, "set");
        
    options.set = false;
    if (isSet(parser, "set"))
        { options.set = true; }
    getArgumentValue(options.file, parser, 0);

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function build_Index() => translate runtime to compile time arguments
// --------------------------------------------------------------------------

template <typename TText, typename TShape, typename TAlg>
void build_Index(TText const & text, TAlg const &, TShape const &)
{
    Index<TText const, IndexSa<Gapped<ModCyclicShape<TShape> > > > index(text);
    indexCreate(index, FibreSA(), TAlg());
    std::cout << "done creating the gapped Index" << std::endl;
}

template <typename TText, typename TAlg>
void build_Index_ungapped(TText const & text, TAlg const &, bool external=false)
{
    Index<TText const, IndexSa<> > index(text);
    if (external)
    {
        resize(indexSA(index), length(indexRawText(index)), Exact());
        _createSuffixArrayPipelining(indexSA(index), indexText(index), Skew7()); // NOTE: ONLY SKEW !!!!
    }
    else
        indexCreate(index, FibreSA(), TAlg());
    std::cout << "done creating an ungapped Index" << std::endl;
}

// --------------------------------------------------------------------------
// Function call() => translate runtime to compile time arguments
// --------------------------------------------------------------------------


template <typename T, typename TShape>
void call2(T const & text, TShape const &, RunSACAOptions const & options)
{

    if (options.algorithm == "DislexSkew7")         build_Index(text, Dislex<Skew7>(), TShape());
    if (options.algorithm == "DislexExternal")      build_Index(text, DislexExternal<TShape>(), TShape());
    if (options.algorithm == "InplaceRadixSort")
    {
        if (TShape::span == 1)                      build_Index_ungapped(text, InplaceRadixSort());
        else                                        build_Index(text, InplaceRadixSort(), TShape());
    }
    if (options.algorithm == "QSort") {
        if (TShape::span == 1)                      build_Index_ungapped(text, SAQSort());
        else                                        build_Index(text, SAQSort(), TShape());
    }
    if (options.algorithm == "Skew7")               build_Index_ungapped(text, Skew7());
    if (options.algorithm == "SkewExternal")        build_Index_ungapped(text, Skew7(), true);
    if (options.algorithm == "None")                std::cout << "No index built." << std::endl;
}


template <typename T>
void call(T const & text, RunSACAOptions const & options)
{
    switch (options.shape) {
        case 0: call2(text, CyclicShape<FixedShape<0,GappedShape<HardwiredShape<> >,0> >(), options); break;
        case 1: call2(text, RunSACAOptions::TShape_1(), options); break;
        case 2: call2(text, RunSACAOptions::TShape_2(), options); break;
        case 3: call2(text, RunSACAOptions::TShape_3(), options); break;
        case 4: call2(text, RunSACAOptions::TShape_4(), options); break;
        case 5: call2(text, RunSACAOptions::TShape_5(), options); break;
        case 6: call2(text, RunSACAOptions::TShape_6(), options); break;
        case 7: call2(text, RunSACAOptions::TShape_7(), options); break;
    }
}

// --------------------------------------------------------------------------
// Function cast_and_run()
// --------------------------------------------------------------------------

template <typename T>
void cast_and_run(T & text, RunSACAOptions const & options)
{
    static bool arr[256]; // 0-initialized
    typedef typename Iterator<T, Standard>::Type TIter;
    for(TIter it=begin(text); it != end(text); ++it)
        arr[ordValue(*it)] = true;
    unsigned sigma = 0;
    for (unsigned i = 0; i< 256; ++i)
        if (arr[i]) ++sigma;

    std::cout << "1 String: total length " << length(text) << std::endl;
    std::cout << "alphabet size: " << sigma << " (assuming 4: DNA, 5: DNA5, <=24: AminoAcid, >24: char)" << std::endl;
    
    if (sigma ==4) { DnaString str; move(str,text); call(str, options); }
    if (sigma ==5) { Dna5String str; move(str,text); call(str, options); }
    if (sigma >5 && sigma < 25) { Peptide str; move(str,text); call(str, options); }
    if (sigma > 24) { CharString str; move(str,text); call(str, options); }
}

template <typename T, typename TSpec>
void cast_and_run(StringSet<T, TSpec> & text, RunSACAOptions const & options)
{
    static bool arr[256]; // 0-initialized
    typedef typename Iterator<typename Concatenator<StringSet<T, TSpec> >::Type, Standard>::Type TIter;
    for(TIter it=begin(concat(text)); it != end(concat(text)); ++it)
        arr[ordValue(*it)] = true;
    unsigned sigma = 0;
    for (unsigned i = 0; i< 256; ++i)
        if (arr[i]) ++sigma;
    
    std::cout << ordValue('N') << "," << ordValue('n') << "," << ordValue('A') << std::endl;
    std::cout << "StringSet " << length(text) << " with length sum " << lengthSum(text) << std::endl;
    std::cout << "alphabet size: " << sigma << " (assuming 4: DNA, 5: DNA5, <=24: AminoAcid, >24: char)" << std::endl;
    
    if (sigma ==4) { StringSet<DnaString> str; move(str,text); call(str, options); }
    if (sigma ==5) { StringSet<Dna5String> str; move(str,text); call(str, options); }
    if (sigma >5 && sigma < 25) { StringSet<Peptide> str; move(str,text); call(str, options); }
    if (sigma > 24) { StringSet<CharString> str; move(str,text); call(str, options); }
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    RunSACAOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    options.print();
    
    
    // Read file
    std::ifstream f(toCString(options.file));
    if (!f.is_open())
    {
        std::cout << "Cannot open " << f << std::endl;
        return 1;
    }
    
    if (options.set)
    {
        StringSet<CharString> set;
        unsigned count = 0;
        std::string line;
        while ( getline (f,line) )
        {
            if (line[0] == '>') {
                ++count;
                resize(set, count);
                continue;
            }
            std::transform(line.begin(), line.end(), line.begin(), ::tolower);
            set[count-1] += static_cast<CharString>(line);
        }
        cast_and_run(set, options);
    }
    else
    {
        CharString str;
        std::string line;
        while ( getline (f,line) )
        {
            if (line[0] == '>')
                continue;
            std::transform(line.begin(), line.end(), line.begin(), ::tolower);
            str += static_cast<CharString>(line);
        }
        cast_and_run(str, options);
    }
    f.close();
    return 0;
}
