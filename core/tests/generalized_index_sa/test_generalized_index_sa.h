// ==========================================================================
//                            generalized_index_sa
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

#ifndef CORE_TESTS_GENERALIZED_INDEX_SA_TEST_GENERALIZED_INDEX_SA_H_
#define CORE_TESTS_GENERALIZED_INDEX_SA_TEST_GENERALIZED_INDEX_SA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

SEQAN_DEFINE_TEST(test_generalized_index_sa_strings_example1)
{
    using namespace seqan;
    
    typedef CyclicShape<GenericShape> TShape;
    typedef ModCyclicShape<TShape>    TModifier;
    typedef Index<CharString, IndexSa<Generalized<TModifier> > > TIndex;

    TShape shape;
    stringToCyclicShape(shape, "1001010");
    
    CharString str = "hallallalalhal al al";
    TIndex index1(str, shape);
    TIndex index2(index1);

    CharString s1; cyclicShapeToString(s1, suffixModifier(index1));
    std::cout << s1 << ": " << indexText(index1) << std::endl;

    CharString s2; cyclicShapeToString(s2, suffixModifier(index2));
    std::cout << s2 << ": " << indexText(index2) << std::endl;
    
}

#endif  // CORE_TESTS_GENERALIZED_INDEX_SA_TEST_GENERALIZED_INDEX_SA_H_
