// ==========================================================================
//                             inplace_radix_sort
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

#ifndef CORE_TESTS_INPLACE_RADIX_SORT_TEST_INPLACE_RADIX_SORT_H_
#define CORE_TESTS_INPLACE_RADIX_SORT_TEST_INPLACE_RADIX_SORT_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
using namespace seqan;

// A test for strings.
void test()
{
    typedef SimpleType<unsigned, Finite<50> > TAlph;
    typedef StringSet<String<TAlph> > TSet;

    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,2,1,2,1,1> >, 0> > TShape; // 0110110111
    typedef IndexSa<Gapped<ModCyclicShape<TShape> > > TIndexSpec;

    TSet s;
    generateText(s, 1000, 20000);
    std::cout << "StringSet with " << length(s) << " strings (" << lengthSum(s) << " characters)" << std::endl;

    Index<TSet, TIndexSpec> i1(s);
    Index<TSet, TIndexSpec> i2(s);
    std::cout << "1." << std::endl;
//indexCreate(i1, FibreSA(), Dislex<>());
    std::cout << "2." << std::endl;

    double x = cpuTime();
    indexCreate(i2, FibreSA(), InplaceRadixSort());
    std::cout << "radixTime: " << cpuTime() - x << std::endl;

    std::cout << "3." << std::endl;
    SEQAN_ASSERT_EQ(indexSA(i1), indexSA(i2));

}



template <unsigned Q, typename TSize>
void benchmarkRadixSortOnRandomSequences(TSize size)
{
    String<SimpleType<unsigned, Finite<Q> > > str;
    generateText(str, size);

    double time = 0;
    for (unsigned i=0; i<3; ++i)
    {
        Index<String<SimpleType<unsigned, Finite<Q> > > > index(str);
        double x = cpuTime();
        indexCreate(index, FibreSA(), InplaceRadixSort());
        time += cpuTime() - x;
    }
    std::cout << Q << "\t" << size << "   \t" << time/3 << std::endl;
}

void benchmarkRadixSortOnRandomSequences()
{
    benchmarkRadixSortOnRandomSequences<4>   (   200000);
    benchmarkRadixSortOnRandomSequences<4>   (   400000);
    benchmarkRadixSortOnRandomSequences<4>   (   800000);
    benchmarkRadixSortOnRandomSequences<4>   (  1600000);
    benchmarkRadixSortOnRandomSequences<4>   (  3200000);
    benchmarkRadixSortOnRandomSequences<4>   (  6400000);
    benchmarkRadixSortOnRandomSequences<4>   ( 12800000);
    benchmarkRadixSortOnRandomSequences<4>   ( 25600000);
    benchmarkRadixSortOnRandomSequences<4>   ( 51200000);
    benchmarkRadixSortOnRandomSequences<4>   (102400000);

    benchmarkRadixSortOnRandomSequences<8>   (   200000);
    benchmarkRadixSortOnRandomSequences<8>   (   400000);
    benchmarkRadixSortOnRandomSequences<8>   (   800000);
    benchmarkRadixSortOnRandomSequences<8>   (  1600000);
    benchmarkRadixSortOnRandomSequences<8>   (  3200000);
    benchmarkRadixSortOnRandomSequences<8>   (  6400000);
    benchmarkRadixSortOnRandomSequences<8>   ( 12800000);
    benchmarkRadixSortOnRandomSequences<8>   ( 25600000);
    benchmarkRadixSortOnRandomSequences<8>   ( 51200000);
    benchmarkRadixSortOnRandomSequences<8>   (102400000);

    benchmarkRadixSortOnRandomSequences<16>  (   200000);
    benchmarkRadixSortOnRandomSequences<16>  (   400000);
    benchmarkRadixSortOnRandomSequences<16>  (   800000);
    benchmarkRadixSortOnRandomSequences<16>  (  1600000);
    benchmarkRadixSortOnRandomSequences<16>  (  3200000);
    benchmarkRadixSortOnRandomSequences<16>  (  6400000);
    benchmarkRadixSortOnRandomSequences<16>  ( 12800000);
    benchmarkRadixSortOnRandomSequences<16>  ( 25600000);
    benchmarkRadixSortOnRandomSequences<16>  ( 51200000);
    benchmarkRadixSortOnRandomSequences<16>  (102400000);

    benchmarkRadixSortOnRandomSequences<32>  (   200000);
    benchmarkRadixSortOnRandomSequences<32>  (   400000);
    benchmarkRadixSortOnRandomSequences<32>  (   800000);
    benchmarkRadixSortOnRandomSequences<32>  (  1600000);
    benchmarkRadixSortOnRandomSequences<32>  (  3200000);
    benchmarkRadixSortOnRandomSequences<32>  (  6400000);
    benchmarkRadixSortOnRandomSequences<32>  ( 12800000);
    benchmarkRadixSortOnRandomSequences<32>  ( 25600000);
    benchmarkRadixSortOnRandomSequences<32>  ( 51200000);
    benchmarkRadixSortOnRandomSequences<32>  (102400000);

    benchmarkRadixSortOnRandomSequences<64>  (   200000);
    benchmarkRadixSortOnRandomSequences<64>  (   400000);
    benchmarkRadixSortOnRandomSequences<64>  (   800000);
    benchmarkRadixSortOnRandomSequences<64>  (  1600000);
    benchmarkRadixSortOnRandomSequences<64>  (  3200000);
    benchmarkRadixSortOnRandomSequences<64>  (  6400000);
    benchmarkRadixSortOnRandomSequences<64>  ( 12800000);
    benchmarkRadixSortOnRandomSequences<64>  ( 25600000);
    benchmarkRadixSortOnRandomSequences<64>  ( 51200000);
    benchmarkRadixSortOnRandomSequences<64>  (102400000);

    benchmarkRadixSortOnRandomSequences<128> (   200000);
    benchmarkRadixSortOnRandomSequences<128> (   400000);
    benchmarkRadixSortOnRandomSequences<128> (   800000);
    benchmarkRadixSortOnRandomSequences<128> (  1600000);
    benchmarkRadixSortOnRandomSequences<128> (  3200000);
    benchmarkRadixSortOnRandomSequences<128> (  6400000);
    benchmarkRadixSortOnRandomSequences<128> ( 12800000);
    benchmarkRadixSortOnRandomSequences<128> ( 25600000);
    benchmarkRadixSortOnRandomSequences<128> ( 51200000);
    benchmarkRadixSortOnRandomSequences<128> (102400000);

    benchmarkRadixSortOnRandomSequences<256> (   200000);
    benchmarkRadixSortOnRandomSequences<256> (   400000);
    benchmarkRadixSortOnRandomSequences<256> (   800000);
    benchmarkRadixSortOnRandomSequences<256> (  1600000);
    benchmarkRadixSortOnRandomSequences<256> (  3200000);
    benchmarkRadixSortOnRandomSequences<256> (  6400000);
    benchmarkRadixSortOnRandomSequences<256> ( 12800000);
    benchmarkRadixSortOnRandomSequences<256> ( 25600000);
    benchmarkRadixSortOnRandomSequences<256> ( 51200000);
    benchmarkRadixSortOnRandomSequences<256> (102400000);

    benchmarkRadixSortOnRandomSequences<512> (   200000);
    benchmarkRadixSortOnRandomSequences<512> (   400000);
    benchmarkRadixSortOnRandomSequences<512> (   800000);
    benchmarkRadixSortOnRandomSequences<512> (  1600000);
    benchmarkRadixSortOnRandomSequences<512> (  3200000);
    benchmarkRadixSortOnRandomSequences<512> (  6400000);
    benchmarkRadixSortOnRandomSequences<512> ( 12800000);
    benchmarkRadixSortOnRandomSequences<512> ( 25600000);
    benchmarkRadixSortOnRandomSequences<512> ( 51200000);
    benchmarkRadixSortOnRandomSequences<512> (102400000);
}

#endif  // CORE_TESTS_INPLACE_RADIX_SORT_TEST_INPLACE_RADIX_SORT_H_
