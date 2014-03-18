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

#include <seqan/basic.h>
#include <seqan/file.h>

#include "../index/test_index_helpers.h"
#include "test_gappedIndex_find.h"
#include "test_gappedIndex_stree.h"



SEQAN_BEGIN_TESTSUITE(test_gappedIndex)
{
    SEQAN_CALL_TEST(test_gappedIndex_top_down_traversal);

	SEQAN_CALL_TEST(test_gappedIndex_find_10_DnaString);
    SEQAN_CALL_TEST(test_gappedIndex_find_10_DnaString_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_10_Dna5String);
	SEQAN_CALL_TEST(test_gappedIndex_find_10_Dna5String_Set);
    SEQAN_CALL_TEST(test_gappedIndex_find_10_Peptide);
    SEQAN_CALL_TEST(test_gappedIndex_find_10_Peptide_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_10_CharString);
	SEQAN_CALL_TEST(test_gappedIndex_find_10_CharString_Set);

    //special
	SEQAN_CALL_TEST(test_gappedIndex_find_10_Finite);
	SEQAN_CALL_TEST(test_gappedIndex_find_10_Finite_Set);

	SEQAN_CALL_TEST(test_gappedIndex_find_11010_DnaString);
    SEQAN_CALL_TEST(test_gappedIndex_find_11010_DnaString_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_11010_Dna5String);
	SEQAN_CALL_TEST(test_gappedIndex_find_11010_Dna5String_Set);
    SEQAN_CALL_TEST(test_gappedIndex_find_11010_Peptide);
    SEQAN_CALL_TEST(test_gappedIndex_find_11010_Peptide_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_11010_CharString);
	SEQAN_CALL_TEST(test_gappedIndex_find_11010_CharString_Set);

	SEQAN_CALL_TEST(test_gappedIndex_find_111100_DnaString);
    SEQAN_CALL_TEST(test_gappedIndex_find_111100_DnaString_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_111100_Dna5String);
	SEQAN_CALL_TEST(test_gappedIndex_find_111100_Dna5String_Set);
    SEQAN_CALL_TEST(test_gappedIndex_find_111100_Peptide);
    SEQAN_CALL_TEST(test_gappedIndex_find_111100_Peptide_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_111100_CharString);
	SEQAN_CALL_TEST(test_gappedIndex_find_111100_CharString_Set);

	SEQAN_CALL_TEST(test_gappedIndex_find_10001_DnaString);
    SEQAN_CALL_TEST(test_gappedIndex_find_10001_DnaString_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_10001_Dna5String);
	SEQAN_CALL_TEST(test_gappedIndex_find_10001_Dna5String_Set);
    SEQAN_CALL_TEST(test_gappedIndex_find_10001_Peptide);
    SEQAN_CALL_TEST(test_gappedIndex_find_10001_Peptide_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_10001_CharString);
	SEQAN_CALL_TEST(test_gappedIndex_find_10001_CharString_Set);

    SEQAN_CALL_TEST(test_gappedIndex_find_01_DnaString);
    SEQAN_CALL_TEST(test_gappedIndex_find_01_DnaString_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_01_Dna5String);
	SEQAN_CALL_TEST(test_gappedIndex_find_01_Dna5String_Set);
    SEQAN_CALL_TEST(test_gappedIndex_find_01_Peptide);
    SEQAN_CALL_TEST(test_gappedIndex_find_01_Peptide_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_01_CharString);
	SEQAN_CALL_TEST(test_gappedIndex_find_01_CharString_Set);

    SEQAN_CALL_TEST(test_gappedIndex_find_0011_DnaString);
    SEQAN_CALL_TEST(test_gappedIndex_find_0011_DnaString_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_0011_Dna5String);
	SEQAN_CALL_TEST(test_gappedIndex_find_0011_Dna5String_Set);
    SEQAN_CALL_TEST(test_gappedIndex_find_0011_Peptide);
    SEQAN_CALL_TEST(test_gappedIndex_find_0011_Peptide_Set);
	SEQAN_CALL_TEST(test_gappedIndex_find_0011_CharString);
	SEQAN_CALL_TEST(test_gappedIndex_find_0011_CharString_Set);

}
SEQAN_END_TESTSUITE
