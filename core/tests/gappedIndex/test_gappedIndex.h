// ==========================================================================
//                                gappedIndex
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

#ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_H_
#define CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_H_

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

// A test for strings.
SEQAN_DEFINE_TEST(test_gappedIndex_strings_example1)
{
    using namespace seqan;

    typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,2> >, 1> > TShape; // 11010
    typedef Index<CharString, IndexSa<Gapped<ModCyclicShape<TShape> > > > TIndex;

    // Define some constant test data for comparison...
    CharString str = "01234567890123456789012345678901234567890123456789"; // 50

    TIndex index(str);
    indexRequire(index, FibreSA());

}

#endif  // CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_H_
