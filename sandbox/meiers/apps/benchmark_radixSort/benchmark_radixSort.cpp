// ==========================================================================
//                            benchmark_radixSort
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

#include <seqan/index.h>
#include "../tests/index/test_index_helpers.h"
using namespace seqan;

// ==========================================================================
// Functions
// ==========================================================================

template <typename TString, typename TSA, typename TSize>
void radixSort1(TString const & text, TSA & sa, TSize depth)
{
    TSA saX;
    resize(saX, length(sa));    // swap memory
    unsigned const K = ValueSize<typename Value<TString>::Type>::VALUE; // alphabet size

    String<TSize, Alloc<> > cnt;
    resize(cnt, K, Exact());	// counter array

    for (unsigned level = depth; level >0; --level)
    {
        swap(sa, saX);
        radixPass(sa, saX,  text,    cnt, K, level);
    }
    swap(sa, saX);
    radixPass(sa, saX, text, cnt, K);
}


template <typename TString, typename TSA, typename TSize>
void radixSort2(TString const & text, TSA & sa, TSize depth)
{
    inplaceRadixSort(sa, text, depth);
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{

    if (argc != 4)
    {
        std::cerr << "Provide  (1) String Size  (2) depth  (3) radix version [1|2]" << std::endl;
        return 1;
    }

    unsigned len = atoi(argv[1]);
    unsigned depth = atoi(argv[2]);
    unsigned version = atoi(argv[3]);

    typedef Dna TAlph;


    String<TAlph> str;
    generateText(str, len);
    String<Size<String<TAlph> >::Type> sa;
    resize(sa, len);
    _initializeSA(sa, str);

    if (version ==1)
    {
        std::cout << "LSD radix sort" << std::endl;
        radixSort1(str, sa, depth);
    }
    if (version == 2)
    {
        std::cout << "MSD radix sort" << std::endl;
        radixSort2(str, sa, depth);
    }
    if (version == 3)
    {
        std::cout << "Only String construction" << std::endl;
    }

    return 0;
}
