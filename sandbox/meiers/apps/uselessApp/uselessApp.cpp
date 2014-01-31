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
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    /*
    // use the sequenceStream
    SequenceStream seqStream(argv[1]);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file." << std::endl;
        return 1;
    }
    if (readAll(ids, seqs, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from file" << std::endl;
        return 1;
    }
     */


    // variable typedefs:
    typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<2> >, 0> > TShape;
    typedef GappedTupler<TShape, false> GappedTupler_;

    CharString shapeString;
    cyclicShapeToString(shapeString, TShape());
    std::cout << "Shape: " << shapeString << std::endl << std::endl;

    {
        // StringSet
        std::cout << "StringSet" << std::endl;
        StringSet<CharString> set;
        appendValue(set, "0123456789");
        appendValue(set, "012");
        appendValue(set, "0");
        appendValue(set, "0123456789");
        appendValue(set, "012");



        typedef Concatenator<StringSet<CharString> >::Type         TConcat;
        typedef StringSetLimits<StringSet<CharString> >::Type      TLimits;
        typedef Multi< GappedTupler_, Pair<Size<CharString>::Type>, TLimits> TMultiGappedTupler;
        typedef Pipe< TConcat, Source<> >                                   TPipeSource;
        typedef Pipe< TPipeSource, TMultiGappedTupler >                     TPipeTupler;

        typedef _dislexTupleCompMulti<TypeOf_(TPipeTupler), TShape, TLimits> _tupleComparator;
        typedef Pool<TypeOf_(TPipeTupler),
                SorterSpec< SorterConfigSize<_tupleComparator, TSizeOf_(TPipeTupler) > > > TPoolSorter;
        typedef Pipe< TPoolSorter, Namer<_tupleComparator> > TPipeNamer;

        typedef _dislexMapMulti<TypeOf_(TPipeNamer), TLimits > _dislexMapper;
        typedef Pool< TypeOf_(TPipeNamer), MapperSpec< MapperConfigSize< _dislexMapper, TSizeOf_(TPipeNamer) > > > TPoolMapper;

        typedef Pipe< TPoolMapper, Filter< filterI2<TypeOf_(TPoolMapper)> > > TFilterI2;

        typedef Pipe<TFilterI2, Skew7> TSkew;

        typedef _dislexReverseTransformMulti<TypeOf_(TSkew), TLimits> _dislexReverse;
        // TODO(meiers): Specify TResult (third tmp arg) to be what the final SA wants.
        typedef Pipe<TSkew, Filter<_dislexReverse> > TPipeReverseTransform;


        
        // begin pipeline

        // Generate Tuples
        TPipeSource src(concat(set));
        TPipeTupler	tupler(src, stringSetLimits(set));

        // Sort Tuples
        _tupleComparator comp(TShape(), stringSetLimits(set));
        TPoolSorter sorter(tupler, comp);
        sorter << tupler;

        // Name sorted Tuples
        TPipeNamer namer(sorter, comp);

        // Map to LexText
        _dislexMapper map(TShape::span, stringSetLimits(set));
        TPoolMapper mapper(namer, map);
        mapper << namer;

        // Discard the original position
        TFilterI2 filter(mapper);

        // Run Skew
        TSkew skew(filter);

        // Reverse Map
        _dislexReverse mapReverse(TShape::span, stringSetLimits(set));
        TPipeReverseTransform reverse(skew, mapReverse);

//        reverseMapper << skew;
          std::cout << reverse;
    }

    {
        // String
        std::cout << "String" << std::endl;
        CharString s = "banaaanaaa";

        typedef Pipe< CharString, Source<> >                    TPipeSource;
        typedef Pipe< TPipeSource, GappedTupler_ >              TPipeTupler;

        typedef _dislexTupleComp<TypeOf_(TPipeTupler), TShape> _tupleComparator;
        typedef Pool<TypeOf_(TPipeTupler),
            SorterSpec< SorterConfigSize<_tupleComparator, TSizeOf_(TPipeTupler) > > > TPoolSorter;

        typedef Pipe< TPoolSorter, Namer<_tupleComparator> > TPipeNamer;

        typedef _dislexMap<TypeOf_(TPipeNamer) > _dislexMapper;
        typedef Pool< TypeOf_(TPipeNamer), MapperSpec< MapperConfigSize< _dislexMapper, TSizeOf_(TPipeNamer) > > > TPoolMapper;

        typedef Pipe< TPoolMapper, Filter< filterI2<TypeOf_(TPoolMapper)> > > TFilterI2;

        typedef Pipe<TFilterI2, Skew7> TSkew;

        typedef _dislexReverseTransform<TypeOf_(TSkew)> _dislexReverse;
        // TODO(meiers): Specify TResult (third tmp arg) to be what the final SA wants.
        typedef Pipe< TSkew, Filter<_dislexReverse> > TPipeReverseTransform;
        
        



        // begin pipeline

        // Generate Tuples
        TPipeSource src(s);
        TPipeTupler	tupler(src);

        // Sort Tuples
        _tupleComparator comp(TShape(), length(s));
        TPoolSorter sorter(tupler, comp);
        sorter << tupler;

        // Name sorted Tuples
        TPipeNamer namer(sorter, comp);

        // Map to LexText
        _dislexMapper map(TShape::span, length(s));
        TPoolMapper mapper(namer, map);
        mapper << namer;



        // Discard the original position
        TFilterI2 filter(mapper);

        // Run Skew
        TSkew skew(filter);

// up to here correct

        // Reverse Map
        _dislexReverse mapReverse(TShape::span, length(s));
        TPipeReverseTransform reverse(skew, mapReverse);

std::cout << reverse << std::endl;
        
    }

    
}

