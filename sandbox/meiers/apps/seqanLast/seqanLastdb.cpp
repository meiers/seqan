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

#include "seqanLast_IO.h"
#include "seqanLast_core.h"

using namespace seqan;

// debug output

template<typename TSA, typename TDir, typename TSeqSet>
void printTables(TSA const & sa, TDir const & dir, TSeqSet const & seqs)
{
    std::cout << "dir      \tsuffix array" << std::endl;
    typedef typename Size<TDir>::Type TSize;
    TSize dirPos = 0;
    TSize saPos = 0;

    while(dirPos < length(dir)-1)
    {
        std::cout << dirPos << ":\t" << dir[dirPos];

        if (saPos < dir[dirPos+1]) {
            std::cout  << "\t\t" << saPos << "\t" << sa[saPos] << "\t" << suffix(seqs, sa[saPos]) << std::endl;
            ++saPos;
        } else {
            std::cout << std::endl;
        }
        while (saPos < dir[dirPos+1])
        {
            std::cout << "\t\t\t" << saPos << "\t" << sa[saPos]  << "\t" << suffix(seqs, sa[saPos])  << std::endl;
            ++saPos;
        }
        ++dirPos;
    }

}




// =============================================================================
// Classes
// =============================================================================

// -----------------------------------------------------------------------------
// Class SeqanLastDbOptions
// -----------------------------------------------------------------------------

struct SeqanLastDbOptions
{
    int verbosity; // 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    CharString databaseFile;
    CharString outputName;
    int shapeChoice;
    int k;

    SeqanLastDbOptions() : verbosity(1), shapeChoice(1)
    {}

    void print()
    {
        std::cout << "Files:" << std::endl;
        std::cout << "   database:    " << databaseFile << std::endl;
        std::cout << "   output name: " << outputName  << std::endl;
        std::cout << "Other:" << std::endl;
        std::cout << "   shape:       " << shapeChoice << std::endl;
        std::cout << "   k:           " << k  << std::endl;
    }
};

// =============================================================================
// Functions
// =============================================================================

template <typename TRealShape, typename TCyclicShape, typename TSize>
void _getKmerShape(TRealShape & shape, TCyclicShape const & cycShape, TSize k)
{
    CharString shapeStr;
    cyclicShapeToString(shapeStr, cycShape);

    for(unsigned i=static_cast<unsigned>(WEIGHT<TCyclicShape>::VALUE); i<static_cast<unsigned>(k); i+= static_cast<unsigned>(WEIGHT<TCyclicShape>::VALUE))
        append(shapeStr, shapeStr);

    unsigned i = 0;
    for (unsigned w=0; i <length(shapeStr) && w < static_cast<unsigned>(k); ++i)
        if (shapeStr[i] == '1') ++w;

    resize(shapeStr, i);
    stringToShape(shape, shapeStr);
}

template <typename TSize, typename TType>
int _writePropertyFile(CharString const & fileName, TSize k, TType shapeChoice)
{
    std::fstream file(toCString(fileName), std::ios::binary | std::ios::out);
    if (!file.good()) {
        std::cerr << "Could not open " << fileName << " to write the propoerty file" << std::endl;
        return 3;
    }
    file << "k=" << k << std::endl;
    file << "shape=" << shapeChoice << std::endl;
    return 0;
}


// -----------------------------------------------------------------------------
// struct Lastdb()
// -----------------------------------------------------------------------------

template <typename TShape, unsigned K, typename TSeqSet, typename TIdSet>
struct Lastdb
{
    typedef typename Value<TSeqSet>::Type   TSeq;
    typedef typename Value<TSeq>::Type      TAlph;

    Lastdb() {}

    int build(TSeqSet const &databases, TIdSet const & ids, SeqanLastDbOptions &options)
    {
        if (options.verbosity > 1) std::cout << "Building Suffix array... " << std::endl;

        typedef Index<TSeqSet const, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;
        TIndex index(databases, TShape());
        indexCreate(index, FibreSA(), SAQSort() ); // TODO(meiers): choose algorithm!
        save(index, toCString(options.outputName));


        if (options.verbosity > 1) std::cout << "Building Look up table for K=" << K << "... " << std::endl;

        typedef Index<TSeqSet const, IndexQGram<Shape<TAlph, GenericShape> > > THashTable;
        THashTable hashTab(databases);
        Shape<TAlph, GenericShape> shape;
        _getKmerShape(shape, TShape(), K);
        indexShape(hashTab) = shape;
        resize(indexDir(hashTab), _fullDirLength(hashTab));
        adaptedCreateQGramIndexDirOnly(indexDir(hashTab), indexBucketMap(hashTab), databases, indexShape(hashTab));
        CharString fileName=options.outputName;
        append(fileName, ".dir");
        save(indexDir(hashTab), toCString(fileName));
        //printTables(indexSA(index), indexDir(hashTab), databases);

        // write the sequence Ids into file
        fileName = options.outputName;
        append(fileName, ".ids");
        std::ofstream file(toCString(fileName), std::ios::out);
        if (file.good())
        {
            for(typename Iterator<TIdSet const>::Type it=begin(ids); it!= end(ids); ++it)
                file << *it << std::endl;
        }
        else
        {
            std::cerr << "Cannot write sequence Ids to " << fileName << std::endl;
            return 4;
        }

        if (options.verbosity > 1) std::cout << "Writing property file... " << std::endl;

        fileName = options.outputName;
        append(fileName, ".prt");
        return _writePropertyFile(fileName, K, options.shapeChoice);
    }
};



// -----------------------------------------------------------------------------
// Function paramaterChoice
// -----------------------------------------------------------------------------

// 3: Run build
template <unsigned K, typename TSeqSet, typename TIdSet, typename TShape>
int paramaterChoice3(TSeqSet const &databases, TIdSet const &ids, SeqanLastDbOptions &options, TShape const &)
{
    Lastdb<TShape, K, TSeqSet, TIdSet> lastdb;
    return lastdb.build(databases, ids, options);
}

// 2: Choose k
template <typename TSeqSet, typename TIdSet, typename TShape>
int paramaterChoice2(TSeqSet const &databases, TIdSet const &ids, SeqanLastDbOptions &options, TShape const &)
{
    switch (options.k)
    {
        case 2:  return paramaterChoice3 <2> (databases, ids, options, TShape());
        case 3:  return paramaterChoice3 <3> (databases, ids, options, TShape());
        case 4:  return paramaterChoice3 <4> (databases, ids, options, TShape());
        case 5:  return paramaterChoice3 <5> (databases, ids, options, TShape());
        case 6:  return paramaterChoice3 <6> (databases, ids, options, TShape());
        case 7:  return paramaterChoice3 <7> (databases, ids, options, TShape());
        case 8:  return paramaterChoice3 <8> (databases, ids, options, TShape());
        case 9:  return paramaterChoice3 <9> (databases, ids, options, TShape());
        case 10: return paramaterChoice3 <10>(databases, ids, options, TShape());
        case 11: return paramaterChoice3 <11>(databases, ids, options, TShape());
        case 12: return paramaterChoice3 <12>(databases, ids, options, TShape());
        default:
            std::cerr << "No valid k-mer size chosen. Exit" << std::endl;
            return 2;
    }
}

// 1: Choose Shape
template <typename TSeqSet, typename TIdSet>
int paramaterChoice1(TSeqSet const &databases, TIdSet const &ids, SeqanLastDbOptions &options)
{
    switch (options.shapeChoice)
    {
        case 1: return paramaterChoice2(databases, ids, options, Shape1());
        case 2: return paramaterChoice2(databases, ids, options, Shape2());
        default:
            std::cerr << "No valid shape chosen. Exit" << std::endl;
            return 1;
    }
}


// -----------------------------------------------------------------------------
// Function main()
// -----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{

    // only Dna5 supported so far
    typedef Dna5 TAlphabet;
    typedef String<TAlphabet>                            TSeq;
    typedef StringSet<TSeq, Owner<> >                    TSeqSet; // TODO: MAke this concat direct ??


    // set option parser
    ArgumentParser parser;
    setShortDescription(parser, "seqanLast: build index");
    setDate(parser, "February 2014");
    setVersion(parser, "0.1");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIFASTA FILE\\fP> <\\fIOUTPUT NAME\\fP>");
    addDescription(parser, "Builds the index for the local aligner SeqanLast (see seqanLast -h for more). ");
    addDescription(parser, "Choose OUTPUT NAME carefully, as many files can be generated. "
                           "We recommend to use a subfolder for them.");
    addDescription(parser, "Note: Only Dna5 supported! Only up to 255 sequences (e.g. chromosomes) in database file allowed");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE"));
    setValidValues(parser, 0, "fa fasta");
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE, "OUTPUT NAME"));
    addOption(parser, ArgParseOption("v", "verbose", "Set verbosity mode."));
    addOption(parser, ArgParseOption("V", "very-verbose", "Set stronger verbosity mode."));
    addOption(parser, ArgParseOption("Q", "quiet", "No output, please."));
    addOption(parser, ArgParseOption("s", "shape", "shape used for the suffix array", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "s", "1");
    setMinValue(parser, "s", "0");
    setMaxValue(parser, "s", "2");
    addOption(parser, ArgParseOption("k", "k-mer", "k-mer size used in the hash table.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "k", "8");
    setMinValue(parser, "k", "2");
    setMaxValue(parser, "k", "12");


    // parse command line
    SeqanLastDbOptions options;
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    getArgumentValue(options.databaseFile, parser, 0);
    getArgumentValue(options.outputName, parser, 1);
    getOptionValue(options.shapeChoice, parser, "shape");
    getOptionValue(options.k, parser, "k-mer");
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    // import database sequence
    TSeqSet databases;
    StringSet<CharString> ids;
    if (!_importSequences(databases, ids, options.databaseFile, options.verbosity))
        return 1;

    // Output Options:
    if(options.verbosity > 1)
        options.print();

    // convert runtime parameters to compileTime parameter
    paramaterChoice1(databases, ids, options);
    
    return 0;
}
