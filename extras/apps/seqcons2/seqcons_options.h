// ==========================================================================
//                                 SeqCons
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef EXTRAS_APPS_SEQCONS2_SEQCONS_OPTIONS_H_
#define EXTRAS_APPS_SEQCONS2_SEQCONS_OPTIONS_H_

#include <iosfwd>
#include <string>

#include <seqan/arg_parse.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeqConsOptions
// ----------------------------------------------------------------------------

struct SeqConsOptions
{
    // Verbosity level, 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose
    int verbosity { 1 };

    // The operation to perform.
    enum Operation
    {
        NOP = 0,            // do nothing, just convert files
        REALIGN = 1,        // just realign, only works with SAM input
        POS_CONSENSUS = 2,  // consensus, interpreting positions
        CTG_CONSENSUS = 3,  // consensus, interpreting contig information
        MSA_CONSENSUS = 4   // "de novo" consensus, just MSA and consensus
    };

    // The operation to perform.
    Operation operation { NOP };

    // -----------------------------------------------------------------------
    // I/O Options
    // -----------------------------------------------------------------------

    // Path to the input file.
    std::string inputFile;
    // Path to the output files for read alignments and resulting consensus.
    std::string outputFileAlignment;
    std::string outputFileConsensus;

    // -----------------------------------------------------------------------
    // Alignment Quality Filter Options
    // -----------------------------------------------------------------------

    // Minimal overlap length.
    int overlapMinLength { 20 };
    // Maximal error rate in percent for a read.
    double overlapMaxErrorRate { 5 };
    // Minimal number of overlap for a read.
    int overlapMinCount { 3 };
    // Window size to look for overlaps when positions are considered.
    int overlapWindowSize { 20 };

    // -----------------------------------------------------------------------
    // K-mer Filter Options
    // -----------------------------------------------------------------------

    // K-mer length to use for identifying overlap candidates.
    int kMerSize { 20 };
    // K-mers with a higher number of occurences are ignored.
    int kMerMaxOcc { 200 };

    // -----------------------------------------------------------------------
    // Realignment Options
    // -----------------------------------------------------------------------

    // The bandwidth to use in the banded alignment when realigning.
    int reAlignmentBandwidth { 10 };
    // The environment for profile extraction.
    int reAlignmentEnvironment { 20 };

    // -----------------------------------------------------------------------
    // Member Functions
    // -----------------------------------------------------------------------

    SeqConsOptions() = default;

    // Check whether the combination of operation and input file is valid.  Throws an exception in the case that it is
    // not.
    void checkConsistency();

    // Print options to a C++ stream.
    void print(std::ostream & out) const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(SeqConsOptions & options, int argc, char const ** argv);

#endif  // #ifndef EXTRAS_APPS_SEQCONS2_SEQCONS_OPTIONS_H_
