// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Code for reading Bam.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_BAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_READ_BAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup SamBamIO SAM/BAM I/O
 * @brief Tags for identifying SAM/BAM format.
 */

/*!
 * @tag SamBamIO#Bam
 * @brief Identify the BAM format.
 *
 * @tag SamBamIO#Sam
 * @brief Identify the SAM format.
 */

/**
.Tag.Bam
..cat:BAM I/O
..signature:Bam
..summary:Tag for identifying the BAM format.
..include:seqan/bam_io.h
..see:Tag.Sam
*/

struct Bam_;
typedef Tag<Bam_> Bam;

template <typename T>
struct FileFormatExtensions<Bam, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileFormatExtensions<Bam, T>::VALUE[1] = {
    ".bam" };


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                              BamHeader
// ----------------------------------------------------------------------------

/*!
 * @fn SamBamIO#readRecord
 * @brief Read a record from a SAM/BAM file.
 *
 * @signature int readRecord(record, context, stream, tag);
 *
 * @param[out]    record  The @link BamAlignmentRecord @endlink object to read the information into.
 * @param[out]    header  The @link BamHeader @endlink object to read the header information into.
 * @param[in,out] context The @link BamIOContext @endlink object to use.
 * @param[in,out] stream  The @link StreamConcept Stream @endlink to read from.
 * @param[in]     tag     The format tag, one of <tt>Sam</tt> and <tt>Bam</tt>.
 *
 * @return int A status code, 0 on success, != 0 on failure.
 */

/**
.Function.readRecord
..signature:readRecord(headerRecord, context, stream, tag)
..param.header:@Class.BamHeader@ to read information into.
...type:Class.BamHeader
..param.context:The context to use for reading.
...type:Class.BamIOContext
..param.stream:The stream to read from (for BAM).
...remarks:BAM data can be read from any stream. For the proper decompression (from compressed BAM, the default) use @Spec.BGZF Stream@.
...type:Concept.StreamConcept
..param.tag:Format to read @Class.BamHeader@ from.
...type:Tag.Sam
...type:Tag.Bam
..include:seqan/bam_io.h
*/

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
void readRecord(BamHeader & header,
               BamIOContext<TNameStore, TNameStoreCache> & context,
               TForwardIter & iter,
               Bam const & /*tag*/)
{
    // Read BAM magic string.
    // char magic[5] = "\0\0\0\0";
    String<char> magic;// = "0000";
    clear(magic);
    readUntil(magic, iter, CountDownFunctor<>(4));
    if (magic != "BAM\1")
        throw std::runtime_error("Not in BAM format.");

    // Read header text, including null padding.
    __int32 lText;
    readRawByte(lText, iter);

    CharString samHeader;
    readUntil(samHeader, iter, CountDownFunctor<>(lText));

    // Truncate to first position of '\0'.
    Iterator<CharString, Rooted>::Type it = begin(samHeader);
    for (; !atEnd(it); ++it)
        if (*it == '\0')
            break;
    resize(samHeader, (it - begin(samHeader))); 

    // Parse out header records.
    BamHeaderRecord headerRecord;
    it = begin(samHeader);
    while (!atEnd(it))
    {
        clear(headerRecord);
        readRecord(headerRecord, context, it, Sam());
        appendValue(header.records, headerRecord);
    }

    // Read # reference sequences.
    __int32 nRef;
    readRawByte(nRef, iter);
    CharString name;

    clear(context.translateFile2GlobalRefId);
    resize(context.translateFile2GlobalRefId, nRef, -1);

    for (__int32 i = 0; i < nRef; ++i)
    {
        // Read length of the reference name.
        __int32 nName;
        readRawByte(nName, iter);
        clear(name);
        readUntil(name, iter, CountDownFunctor<>(nName));
        resize(name, nName - 1);
        // Read length of the reference sequence.
        __int32 lRef;
        readRawByte(lRef, iter);

        // Store sequence info.
        typedef typename BamHeader::TSequenceInfo TSequenceInfo;
        appendValue(header.sequenceInfos, TSequenceInfo(name, lRef));
        // Append contig name to name store, if not known already.
        typename Size<TNameStore>::Type globalRId = 0;
        if (!getIdByName(nameStore(context), name, globalRId, nameStoreCache(context)))
        {
            globalRId = length(nameStore(context));
            appendName(nameStore(context), name, nameStoreCache(context));
        }
        context.translateFile2GlobalRefId[i] = globalRId;
    }
}

// ----------------------------------------------------------------------------
// Function readRecord()                                     BamAlignmentRecord
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..signature:readRecord(alignmentRecord, context, stream, tag)
..param.alignmentRecord.type:Class.BamAlignmentRecord
*/

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
void readRecord(BamAlignmentRecord & record,
                BamIOContext<TNameStore, TNameStoreCache> & context,
                TForwardIter & iter,
                Bam const & /*tag*/)
{
    (void)context;  // Only used for assertions.

    // Read size of the remaining block.
    __int32 remainingBytes = 0;
    readRawByte(remainingBytes, iter);

    SEQAN_ASSERT_GT(remainingBytes, 32);
    readRawByte(static_cast<BamAlignmentRecordCore &>(record), iter);
    remainingBytes -= 32;

    // Translate file local rID into a global rID that is compatible with the context nameStore.
    if (record.rID >= 0 && !empty(context.translateFile2GlobalRefId))
        record.rID = context.translateFile2GlobalRefId[record.rID];

    if (record.rID >= 0)
        SEQAN_ASSERT_LT(static_cast<__uint64>(record.rID), length(nameStore(context)));

    // read name.
    SEQAN_ASSERT_GT(remainingBytes, (int)record._l_qname);
    //resize(record.qName);
    readUntil(record.qName, iter, CountDownFunctor<>(record._l_qname));
    resize(record.qName, record._l_qname - 1);
    remainingBytes -= record._l_qname;

    // cigar string.
    SEQAN_ASSERT_GT(remainingBytes, record._n_cigar * 4);
    resize(record.cigar, record._n_cigar, Exact());
    static char const * CIGAR_MAPPING = "MIDNSHP=";
    typedef typename Iterator<String<CigarElement<> >, Rooted>::Type TCigarIter;
    for (TCigarIter it = begin(record.cigar, Rooted()); !atEnd(it); goNext(it))
    {
        __uint32 ui = 0;
        readRawByte(ui, iter);
        it->operation = CIGAR_MAPPING[ui & 0x0007];
        it->count = ui >> 4;
    }
    remainingBytes -= record._n_cigar * 4;

    // sequence, 4-bit encoded "=ACMGRSVTWYHKDBN".
    SEQAN_ASSERT_GT(remainingBytes, (record._l_qseq + 2) / 2);
    resize(record.seq, record._l_qseq + 1, Exact());
    static char const * SEQ_MAPPING = "=ACMGRSVTWYHKDBN";

    typedef typename Iterator<CharString, Rooted>::Type TSeqIter;
    {
        // Note: Yes, we need separate index i and iterator.  The iterator allows the fast iteration and i is for
        // book-keeping since we potentially create too long seq records.
        TSeqIter it = begin(record.seq, Rooted());
        for (__int32 i = 0; i < record._l_qseq; i += 2)
        {
            __uint8 ui;
            readRawByte(ui, iter);
            *it++ = SEQ_MAPPING[ui >> 4];
            *it++ = SEQ_MAPPING[ui & 0x0f];
        }
    }
    resize(record.seq, record._l_qseq);  // Possibly trim last, overlap base.
    remainingBytes -= (record._l_qseq + 1) / 2;

    // phred quality
    SEQAN_ASSERT_GEQ(remainingBytes, record._l_qseq);
    if (record._l_qseq > 0)
    {
        clear(record.qual);
        readUntil(record.qual, iter, CountDownFunctor<>(record._l_qseq));
    }
    // If qual is a sequence of 0xff (heuristic same as samtools: Only look at first byte) then we clear it, to get the
    // representation of '*';
    if (!empty(record.qual) && record.qual[0] == '\xFF')
        clear(record.qual);
    typedef typename Iterator<CharString, Rooted>::Type TQualIter;
    for (TQualIter it = begin(record.qual, Rooted()); !atEnd(it); goNext(it))
        *it += '!';
    remainingBytes -= record._l_qseq;

    // tags
    if (remainingBytes > 0)
    {
        //resize(record.tags, remainingBytes);
        readUntil(record.tags, iter, CountDownFunctor<>(remainingBytes));
    }
    else
    {
        clear(record.tags);
    }
}
}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_BAM_H_
