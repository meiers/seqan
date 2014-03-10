// ==========================================================================
//                          test_gappedIndex_stree.h
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

#ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_STREE_H_
#define CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_STREE_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


void testManualExample()
{
    typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<3> >, 0> >      TShape; //1001
    typedef Index<Peptide, IndexEsa< > >  TIndex;        // Gapped<ModCyclicShape<TShape> >

    Peptide s = "KeinerKleinerFeinerReinerHeiligerWeicherNeider";
    TIndex index(s);            // , TShape()
    indexRequire(index, FibreSA());

    typedef Suffix<TIndex>::Type TModSuff;

    // Suffix array (ouput of:)
    //for (unsigned i=0; i<length(indexSA(index)); ++i)
    //    std::cout << "//  " << i << "\t" << indexSA(index)[i] << "\t" << suffix(index, indexSA(index)[i]) << std::endl;
    //
    //  0	45	R                       *
    //  1	19	RNEEIGEEIERID
    //  2	5	REIRFNEEIRHLIRWCHNEER
    //  3	12	RINREERILERICRNDE
    //  4	18	RINHEIGWEHEEIR
    //  5	39	RID
    //  6	32	RICRNDE
    //  7	24	RILERICRNDE
    //  8	16	NREERILERICRNDE         *
    //  9	40	NDE
    //  10	22	NHEIGWEHEEIR
    //  11	3	NKLNEEIRRNEEIGEEIERID
    //  12	10	NFEERINHEIGWEHEEIR
    //  13	43	D                       *
    //  14	36	CRNDE                   *
    //  15	44	E                       *---+
    //  16	41	EER                         *---+---+
    //  17	1	EEREIRFNEEIRHLIRWCHNEER             *
    //  18	8	EERINREERILERICRNDE                 *---+---+
    //  19	14	EERINHEIGWEHEEIR
    //  20	20	EERILERICRNDE                               *
    //  21	38	EEIR                            *---+---+---+
    //  22	11	EEIRRNEEIGEEIERID
    //  23	17	EEIRHLIRWCHNEER
    //  24	31	EEIERID
    //  25	23	EEIGEEIERID
    //  26	34	EHEEIR                      *---+
    //  27	26	EIGWEHEEIR                  *
    //  28	4	ELEERINREERILERICRNDE       *
    //  29	30	GWEHEEIR                *---+
    //  30	37	HNEER                   *
    //  31	25	HLIRWCHNEER
    //  32	42	IR                      *
    //  33	15	IRRNEEIGEEIERID
    //  34	21	IRHLIRWCHNEER
    //  35	2	IRKINFEERINHEIGWEHEEIR
    //  36	9	IRFNEEIRHLIRWCHNEER
    //  37	29	IRWCHNEER
    //  38	35	IERID
    //  39	27	IGEEIERID
    //  40	7	LNEEIRRNEEIGEEIERID     *
    //  41	28	LERICRNDE
    //  42	0	KNELEERINREERILERICRNDE *
    //  43	6	KINFEERINHEIGWEHEEIR
    //  44	13	FNEEIRHLIRWCHNEER       *
    //  45	33	WCHNEER                 *


/*
    Iterator<TIndex, TopDown<> >::Type treeIter(index);

    goDown(treeIter, 'E');
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(15,29));

    goDown(treeIter, 'E');
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(16,26));

    SEQAN_ASSERT_EQ(representative(treeIter), "EE");
    SEQAN_ASSERT_EQ(repLength(treeIter), 2u);

    goDown(treeIter, 'R');
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(16,21));

    goRight(treeIter);
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(21,26));

    goRoot(treeIter);
    goDown(treeIter, "EEIR");
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(21,24));
    SEQAN_ASSERT_EQ(representative(treeIter), "EEIR");
    SEQAN_ASSERT_EQ(repLength(treeIter), 4u);
*/
    // TODO: test isFirstChild _isLeaf and other functions from index_sa_stree.h
    // what about goUp, nodeUp

    Iterator<TIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(index);

    while (!atEnd(myIterator))
    {
        std::cout << representative(myIterator) << std::endl;
        ++myIterator;
    }

    return 0;
}


#endif  // #ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_STREE_H_
