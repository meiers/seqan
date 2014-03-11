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
    typedef SimpleType<unsigned, Finite<22> > TAlph;
    typedef StringSet<String<TAlph> > TSet;

    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,2,1,2,1,1> >, 0> > TShape; // 0110110111
    typedef IndexSa<Gapped<ModCyclicShape<TShape> > > TIndexSpec;

    TSet s;
    generateText(s, 100, 1000);
    std::cout << "StringSet with " << length(s) << " strings (" << lengthSum(s) << " characters)" << std::endl;

    Index<TSet, TIndexSpec> i1(s);
    Index<TSet, TIndexSpec> i2(s);
    std::cout << "1." << std::endl;
    indexCreate(i1, FibreSA(), Dislex<>());
    std::cout << "2." << std::endl;
    indexCreate(i2, FibreSA(), InplaceRadixSort());

    std::cout << "3." << std::endl;
    SEQAN_ASSERT_EQ(indexSA(i1), indexSA(i2));

}

#endif  // CORE_TESTS_INPLACE_RADIX_SORT_TEST_INPLACE_RADIX_SORT_H_
