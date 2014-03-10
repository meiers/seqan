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
#include "../index/test_index_helpers.h"


// TODO(meiers): Completely renew this!


SEQAN_DEFINE_TEST(test_generalized_index_sa_Demo)
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

    CharString s1; cyclicShapeToString(s1, shape);
    std::cout << s1 << ": " << indexText(index1) << std::endl;

    CharString s2; cyclicShapeToString(s2, index2.modifierCargo);
    std::cout << s2 << ": " << indexText(index2) << std::endl;


    //indexRequire(index1, FibreSA());
}


SEQAN_DEFINE_TEST(test_generalized_index_sa_qsort_modcyclicshape)
{
    using namespace seqan;

    typedef StringSet<CharString> TSet;
    typedef CyclicShape<GenericShape> TShape;
    typedef ModCyclicShape<TShape>    TModifier;
    typedef Index<TSet, IndexSa<Generalized<TModifier> > > TGappedIndex;

    TSet set;
    appendValue(set, "012341");
    appendValue(set, "12341");
    appendValue(set, "01011");

    unsigned correct[16][2]={   {2,4}, {1,4}, {0,5}, {2,1}, {2,3}, {2,2}, {1,3}, {0,4},
                                {2,0}, {0,0}, {1,0}, {0,1}, {1,1}, {0,2}, {1,2}, {0,3} };

    TShape shape;
    stringToCyclicShape(shape, "01");

    TGappedIndex gIndex(set, shape);
    resize(indexSA(gIndex), lengthSum(set));
    createGeneralizedSuffixArray(indexSA(gIndex), set, shape, TModifier(), SAQSort());

    for (unsigned i=0; i< 16; ++i)
    {
        SEQAN_ASSERT_EQ(correct[i][0], indexSA(gIndex)[i].i1);
        SEQAN_ASSERT_EQ(correct[i][1], indexSA(gIndex)[i].i2);
    }
}

SEQAN_DEFINE_TEST(test_generalized_index_sa_qsort_withoutmod)
{
    using namespace seqan;

    typedef StringSet<CharString> TSet;
    typedef Index<TSet, IndexSa<> > TNormalIndex;

    // String Set
    TSet set;
    appendValue(set, "012345");
    appendValue(set, "12345");
    appendValue(set, "0123456789");
    appendValue(set, "12345");

    TNormalIndex index1(set);
    resize(indexSA(index1), lengthSum(set));
    createSuffixArray(indexSA(index1), set, SAQSort());

    TNormalIndex index2(set);
    resize(indexSA(index2), lengthSum(set));
    createSuffixArray(indexSA(index2), set, Skew7());

    SEQAN_ASSERT_EQ(indexSA(index1), indexSA(index2));

    // String
    CharString str = "aklsdjfjklfjklsdjfjklaghkladshgkladjsfjladjsfjklf";
    Index<CharString, IndexSa<> > index3(str);
    resize(indexSA(index3), length(str));
    createSuffixArray(indexSA(index3), str, SAQSort());
    SEQAN_ASSERT( isSuffixArray(indexSA(index3), str) );
}



SEQAN_DEFINE_TEST(test_generalized_index_sa_qsort_modreverse)
{
    using namespace seqan;

    typedef StringSet<CharString> TSet;
    typedef Index<TSet, IndexSa<Generalized<ModReverse> > > TRevIndex;

    TSet set;
    appendValue(set, "010101");
    appendValue(set, "101");
    appendValue(set, "000");
    appendValue(set, "1000");

    TRevIndex revIndex(set);
    typename Cargo<TRevIndex>::Type xxx;

    resize(indexSA(revIndex), lengthSum(set));
    createGeneralizedSuffixArray (indexSA(revIndex), set, xxx, ModReverse(), SAQSort());

    unsigned correct[16][2]={   {3,3}, {2,2}, {3,2}, {2,1}, {3,1}, {2,0}, {3,0}, {1,2},
                                {0,5}, {1,1}, {0,4}, {1,0}, {0,3}, {0,2}, {0,1}, {0,0} };
    for (unsigned i=0; i< 16; ++i)
    {
        SEQAN_ASSERT_EQ(correct[i][0], indexSA(revIndex)[i].i1);
        SEQAN_ASSERT_EQ(correct[i][1], indexSA(revIndex)[i].i2);
    }
}

SEQAN_DEFINE_TEST(test_generalized_index_sa_qsort_modview)
{
    using namespace seqan;

    typedef ModView<FunctorComplement<Dna> > TMod;
    typedef Index<DnaString, IndexSa<Generalized<TMod > > > TViewIndex;

    DnaString dna = "ACGTAACGC";
    Cargo<Suffix<TViewIndex>::Type>::Type cargo;

    TViewIndex index(dna, cargo);
    indexRequire(index, FibreSA());


    // correct Index:
    typedef ModifiedString<DnaString,TMod> TModString;
    typedef Index<DnaString, IndexSa<> > TIndex;
    DnaString newStr = TModString(dna);
    TIndex corrIndex(newStr);
    indexRequire(corrIndex, FibreSA());

    for (unsigned i =0; i<length(dna);++i)
        SEQAN_ASSERT_EQ(indexSA(index), indexSA(corrIndex));
    
}







// find

SEQAN_DEFINE_TEST(test_generalized_index_sa_find_modreverse)
{
    using namespace seqan;

    typedef StringSet<CharString> TSet;
    typedef Index<TSet, IndexSa<Generalized<ModReverse> > > TRevIndex;

    TSet set;
    appendValue(set, "010101");
    appendValue(set, "101");
    appendValue(set, "000");
    appendValue(set, "1000");

    TRevIndex revIndex(set);
    typename Cargo<TRevIndex>::Type xxx;

    Finder<TRevIndex> finder(revIndex);
    while (find(finder, "0001"))
        std::cout << beginPosition(finder) << std::endl;

    }


SEQAN_DEFINE_TEST(test_generalized_index_sa_find_modcyclicshape)
{
    using namespace seqan;

    typedef StringSet<CharString> TSet;
    typedef CyclicShape<GenericShape> TShape;
    typedef ModCyclicShape<TShape>    TModifier;
    typedef Index<TSet, IndexSa<Generalized<TModifier> > > TGappedIndex;

    TSet set;
    appendValue(set, "012341");
    appendValue(set, "12341");
    appendValue(set, "01011");

    TShape shape;
    stringToCyclicShape(shape, "01");

    TGappedIndex gIndex(set, shape);

    Finder<TGappedIndex> finder(gIndex);
    while(find(finder, "24"))
        std::cout << beginPosition(finder) << std::endl;

}






// stree interface

SEQAN_DEFINE_TEST(test_generalized_index_sa_stree_modcyclicshape)
{
    using namespace seqan;

    typedef StringSet<CharString>                       TSet;
    typedef CyclicShape<GenericShape>                   TShape;
    typedef ModCyclicShape<TShape>                      TModifier;
    typedef Index<TSet, IndexSa<Generalized<TModifier> > > TGappedIndex;
    typedef Iterator<CharString, Standard>::Type        TQueryIter;
    typedef Iterator<TGappedIndex, TopDown<> >::Type    TTreeIter;

    TSet set;
    appendValue(set, "012341");
    appendValue(set, "12341");
    appendValue(set, "01011");
    CharString  query = "24";

    TShape shape;
    stringToCyclicShape(shape, "01");

    TGappedIndex index(set, shape);
    TTreeIter treeIter(index);
    TQueryIter qry = begin(query, Standard());
    TQueryIter qryEnd = end(query, Standard());
    while(qry < qryEnd)
    {
        std::cout << repLength(treeIter) << " : " << countOccurrences(treeIter) << "\t\"" << representative(treeIter) << "\"" << std::endl;
        if(!goDown(treeIter, *(qry++)))
            break;
    }

    SEQAN_ASSERT_EQ(representative(treeIter), "24");
    SEQAN_ASSERT_EQ(countOccurrences(treeIter), 2u);
}


#endif  // CORE_TESTS_GENERALIZED_INDEX_SA_TEST_GENERALIZED_INDEX_SA_H_
