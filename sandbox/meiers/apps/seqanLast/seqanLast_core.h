// ==========================================================================
//                              seqanLast_core.h
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

#ifndef SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_CORE_H_
#define SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_CORE_H_

// =============================================================================
// Global Definitions
// =============================================================================

// test: count vector of adaptive seeds
String<unsigned> globalCounts;

// StringSet
typedef StringSet<String<Dna5, MMap<> >, Owner<ConcatDirect<> > >         TMMapStringSet;
typedef StringSet<String<Dna5>,          Owner<ConcatDirect<> > >         TNormalStringSet;

// Shapes
typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1> >, 1> >       Shape1;   // 110
typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,1> >, 1> >     Shape2;   // 1110
typedef CyclicShape<FixedShape<0, ShapePatternHunter, 2> >                    Shape3;   // 11101001010011011100
typedef CyclicShape<FixedShape<0, ShapeIlieB3, 2> >                           Shape4;   // 11100010010000101011

// TODO(meiers): Define more shapes!

namespace SEQAN_NAMESPACE_MAIN
{
    // Define size type of Index<StringSet>
    template<>
    struct SAValue<TMMapStringSet>
    {
        typedef Pair<unsigned char, unsigned int, Pack> Type;
    };
    template<>
    struct SAValue<TMMapStringSet const>
    {
        typedef Pair<unsigned char, unsigned int, Pack> Type;
    };
    template<>
    struct SAValue<TNormalStringSet>
    {
        typedef Pair<unsigned char, unsigned int, Pack> Type;
    };
    template<>
    struct SAValue<TNormalStringSet const>
    {
        typedef Pair<unsigned char, unsigned int, Pack> Type;
    };
}

// =============================================================================
// Tags, Classes, Enums
// =============================================================================

/*!
 * @class SeqanLast SeqanLast App
 * @headerfile seqanLast/seqanLast_core.h
 *
 * @brief A local alignment search tool implemented in SeqAn.
 *
 * @signature struct SeqanLast;
 *
 * This is the documnetation of the local alignment search tool seqanLast,
 * which is a reimplementation of Last [Kielbasa et al., 2011].
 *
 * @section Overview
 *
 * Lalalala
 *
 */


// -----------------------------------------------------------------------------
// Struct DiagonalTable
// -----------------------------------------------------------------------------

template <typename TSigned, typename TSize, unsigned Q=256>
struct DiagonalTable
{
    typedef typename Iterator<String<Pair<TSigned,TSize> > >::Type TIter;
    String<String<Pair<TSigned, TSize> > > table;

    DiagonalTable()
    {
        resize(table, Q, Exact());
        for (unsigned i=0; i < Q; ++i)
            reserve(table[i], 100   );
    }

    void clear()
    {
        for (unsigned i=0; i < Q; ++i)
            ::seqan::clear(table[i]);
    }

    bool redundant(TSize pGenome, TSize pQuery)
    {
        TSigned d = static_cast<TSigned>(pGenome) - static_cast<TSigned>(pQuery);
        TSigned dd =(Q + d%Q)%Q; // hope the modulo is replaced by fast bit operations.
        TSize l=0;

        for(TIter it=begin(table[dd]); it != end(table[dd]); )
        {
            if(it->i2 >= pQuery) 
            {
                if(it->i1 == d) return true;
                ++it;
                ++l;
            }
            else // clean up old entries right away
            {
                erase(table[dd], l);
                //it = begin(table[dd])+l; // this happens automaticly
            }     
        }
        return false;
    }

    void add(TSize pGenome, TSize pQuery)
    {
        TSigned d = static_cast<TSigned>(pGenome) - static_cast<TSigned>(pQuery);
        TSigned dd =(Q + d%Q)%Q; // hope the modulo is replaced by fast bit operations.

/*
        for(TIter it=begin(table[dd]); it != end(table[dd]); ++it)
        {
            if(it->i1 != d) continue;
            if(it->i2 <= pQuery)
            {
                // remove smaller entries
                it->i2 = pQuery;
                break;
                //typename Size<String<Pair<TSigned, TSize> > >::Type pos = it - begin(table[dd]);
                //erase(table[dd], pos);
                // restore iterator (seems not even to be neccessary since erase does not realloc memory)
                //it = begin(table[dd]) + (pos-1);
            }
        }
*/
        appendValue(table[dd], Pair<TSigned,TSize>(d,pQuery));
    }

    void ___printSize()
    {
        for(unsigned d=0; d<Q-1; ++d)
            std::cout << length(table[d]) << ",";
        std::cout << length(table[Q-1]) << std::endl;
    }
};


// -----------------------------------------------------------------------------
// Struct SeqanLastMatch
// -----------------------------------------------------------------------------

template <typename TSize, typename TAlign, typename TScore=int>
struct SeqanLastMatch
{
    typedef TAlign Type;
    TAlign align;
    TSize  dbId;
    TSize  quId;
    TScore score;
};

// -----------------------------------------------------------------------------
// Struct MatchScoreLess
// -----------------------------------------------------------------------------

template <typename TMatch>
struct MatchScoreLess : std::binary_function<TMatch const &, TMatch const &, bool>
{
    inline bool operator()(TMatch const & a, TMatch const & b)
    {
        return a.score > b.score;
    }
};

// -----------------------------------------------------------------------------
// Struct LastParameters
// -----------------------------------------------------------------------------

template <typename TSize, typename TScoreMatrix>
struct LastParameters
{
    typedef typename Value<TScoreMatrix>::Type TScore;

    TSize                   maxFreq;
    TScoreMatrix const &    scoreMatrix;
    TScore                  Xgapless;
    TScore                  Xgapped;
    TScore                  Tgapless;
    TScore                  Tgapped;
    bool                    onlyUngappedAlignments;
    bool                    useHashTable;
    bool                    myUngappedExtend;
    int                     verbosity;

    LastParameters(TSize f, TScoreMatrix const & scoreMatrix, TScore Xgapless,
                   TScore Xgapped, TScore Tgapless, TScore Tgapped, bool ungAl,
                   bool useHash, bool myExtend, int v) :
        maxFreq(f), scoreMatrix(scoreMatrix), Xgapless(Xgapless), Xgapped(Xgapped),
        Tgapless(Tgapless), Tgapped(Tgapped), onlyUngappedAlignments(ungAl),
        useHashTable(useHash), myUngappedExtend(myExtend), verbosity(v)
    {}
};



// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
// Function adaptedCreateQGramIndexDirOnly()
// -----------------------------------------------------------------------------

template <typename TDir, typename TBucketMap, typename TText, typename TShape>
void _insertMissingQGrams(TDir &dir,
                          TBucketMap &bucketMap,
                          TText const &text,
                          TShape &shape)
{
    typedef typename Size<TText>::Type TPos;
    typedef typename Iterator<TText const,Standard>::Type TIter;

    if (!length(text)) return;

    std::cout << shape.span << std::endl;
    TPos start = length(text) < shape.span ? 0 : length(text)-shape.span+1;
    TIter it = begin(text, Standard()) + start;
    ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];

    for(++start, ++it; start < length(text); ++start, ++it)
    {
        ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];
    }
}

template <typename TDir, typename TBucketMap, typename TText, typename TSetSpec, typename TShape>
void _insertMissingQGrams(TDir &dir,
                          TBucketMap &bucketMap,
                          StringSet<TText, TSetSpec> const &textSet,
                          TShape &shape)
{
    typedef typename Size<TText>::Type TPos;
    typedef typename Iterator<TText const,Standard>::Type TIter;
    typedef typename Iterator<StringSet<TText, TSetSpec> const,Standard>::Type TSetIter;

    if (!length(textSet)) return;

    for (TSetIter set=begin(textSet,Standard()); set != end(textSet,Standard()); ++set)
    {
        TText const & text = *set;
        if (!length(text)) continue;

        TPos start = length(text) < shape.span ? 0 : length(text)-shape.span+1;
        TIter it = begin(text, Standard()) + start;
        ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];

        for(++start, ++it; start < length(text); ++start, ++it)
        {
            ++dir[requestBucket(bucketMap, hash(shape, it, length(text)-start))];
        }
    }
}

template <typename TDir, typename TBucketMap, typename TText, typename TShape>
void adaptedCreateQGramIndexDirOnly(TDir &dir,
                                    TBucketMap &bucketMap,
                                    TText const &text,
                                    TShape &shape)
{
    // 1. clear counters
    _qgramClearDir(dir, bucketMap);

    // 2. count q-grams
    _qgramCountQGrams(dir, bucketMap, text, shape, 1);

    // New part: Add Q-1 last q-grams (that are usually missed) to the count vector
    _insertMissingQGrams(dir, bucketMap, text, shape);
    
    // 3. cumulative sum (Step 4 is ommited)
    _qgramCummulativeSumAlt(dir, False());
}


// =============================================================================
// Functions of Last
// =============================================================================

// -----------------------------------------------------------------------------
// Function _fastTableLookup()
// -----------------------------------------------------------------------------

template <typename TTrieIt, typename TLookupTable, typename TQueryIt, typename TSize>
inline void _fastTableLookup(TTrieIt & trieIt,
                             TLookupTable & table,
                             TQueryIt & qryIt,      // will be modified
                             TQueryIt & qryEnd,
                             TSize maxFreq)
{
    typedef typename Fibre<TLookupTable, FibreShape>::Type TShape;
    typedef typename Value<TShape>::Type THashValue;
    typedef typename Size<TLookupTable>::Type TSaPos;

    TSaPos restLen = std::min(static_cast<TSaPos>(length(indexShape(table))), static_cast<TSaPos>(qryEnd - qryIt));
    if (restLen < 1) return;

    // determine hash values
    THashValue x = hash(indexShape(table), qryIt, restLen);
    THashValue y = x+1;
    if (restLen < length(indexShape(table)))
        y = hashUpper(indexShape(table), qryIt, restLen);

    
    // get range in SA
    TSaPos from  = indexDir(table)[x];
    TSaPos to    = indexDir(table)[y];


    // SA contains incomplete q-grams. These must be skipped:
    while( suffixLength(indexSA(container(trieIt))[from], container(trieIt)) < length(indexShape(table))
           && from + maxFreq < to)
        ++from;

    if (to - from > maxFreq)
    {
        value(trieIt).range.i1 = from;
        value(trieIt).range.i2 = to;
        value(trieIt).repLen = restLen;
        goFurther(qryIt, restLen-1);
        value(trieIt).lastChar = *qryIt++;
        // TODO: set parentRight? I think I don't need it because I won't goRight()
    }
    else
    {
        // TODO: Go Up in the table
    }
}

/*


// OR: make seed shorter
else
{
    std::cout << "Use new function!!" << std::endl;
    THashValue sigma, sigmaToQ;
    sigma = sigmaToQ = static_cast<THashValue>(ValueSize<typename Host<TShape>::Type>::VALUE);
    unsigned len = std::min(static_cast<TSaPos>(qryEnd - qryIt), weight(indexShape(table)));

    // for len = q-1 to 0
    for(--len; len > 0; --len, sigmaToQ *= sigma)
    {
        // hash values of the substring one shorter
        x = (x/sigmaToQ) * sigmaToQ;
        y = (x/sigmaToQ + 1) * sigmaToQ;

        // if the range is > maxFreq, the previous range was correct
        if (indexDir(table)[y] - indexDir(table)[x] > maxFreq)
        {
            value(trieIt).range.i1 = from;
            value(trieIt).range.i2 = to;
            value(trieIt).repLen = len+1;
            value(trieIt).lastChar = *(qryIt + len);
        }
        from = indexDir(table)[x];
        to   = indexDir(table)[y];
    }
}
*/


// -----------------------------------------------------------------------------
// Function _goDownTrie()
// -----------------------------------------------------------------------------

template <typename TTrieIt, typename TQueryIt, typename TSize>
inline void _goDownTrie(TTrieIt & trieIt,
                        TQueryIt qryIt,
                        TQueryIt qryEnd,
                        TSize maxFreq)
{
    // test: count vector of adaptive seeds
    unsigned x =0;


    while(qryIt < qryEnd)
    {
        if(!goDown(trieIt, *(qryIt++)))
            break;
        if(countOccurrences(trieIt) <= maxFreq)
            break;

        // test: count vector of adaptive seeds
        if (length(globalCounts)<= x)
        {
            resize(globalCounts, x+1);
            globalCounts[x] = 1;
        } else {
            globalCounts[x]++;
        }
        ++x;
    }
}

// -----------------------------------------------------------------------------
// Function adaptiveSeeds()                                           [ungapped]
// -----------------------------------------------------------------------------

/*!
 * @fn SeqanLast#adaptiveSeeds
 * @headerfile seqanLast/seqanLast_core.h
 *
 * @brief Find adaptive seeds to a query using a (gapped) suffix array
 *
 * @signature Pair<TSize> adaptiveSeeds(index, lookupTable, query, maxFreq);
 * @signature Pair<TSize> adaptiveSeeds(index, query, maxFreq);
 *
 * @param[in] index A trie-like index that can be traversed using goDown(). 
 *             Currently this parameter is limited to the IndexSa class
 *             specialization and its subclasses. For gapped adaptive seeds
 *             provide the IndexSa<Gapped<...> > index.
 * @param[in] lookupTable A hash table that supports the fast lookup of short
 *             strings in the suffix array. The q-gram index specialization
 *             is meant. Watch that the positions of the indexDir match those
 *             in the suffix array - usually the indexDir does not contain
 *             incomplete q-grams. This table significantly speeds up the
 *             determination of adaptive seeds.
 * @param[in] maxFreq maximum allowed frequency of a seed to match in the target
 *             sequence.
 * @return Range of occurrences in the suffix array. If adaptive seeds exist
 *             for this query and frequency, this range will be less than
 *             or equal maxFreq. If the range is zero or larger than maxFreq,
 *             no adaptive seeds exist according to the definition of adaptive
 *             seeds.
 */

// Ungapped seeds:

template <typename TTrieIndex, typename TLookupIndex, typename TQuery, typename TSize>
inline Pair<typename Size<TTrieIndex>::Type>
adaptiveSeeds(TTrieIndex   & index,
              TLookupIndex & table,
              TQuery const & query,
              TSize maxFreq)
{
    typedef typename Iterator<TTrieIndex, TopDown<> >::Type   TTreeIter;
    typedef typename Iterator<TQuery const, Standard>::Type   TQueryIter;

    TTreeIter   trieIt(index);
    TQueryIter  qryIt  = begin(query, Standard());
    TQueryIter  qryEnd = end(query, Standard());

    _fastTableLookup(trieIt, table, qryIt, qryEnd, maxFreq);
    _goDownTrie(trieIt, qryIt, qryEnd, maxFreq);
    return range(trieIt);
}

template <typename TTrieIndex, typename TQuery, typename TSize>
inline Pair<typename Size<TTrieIndex>::Type>
adaptiveSeeds(TTrieIndex   & index,
              TQuery const & query,
              TSize maxFreq)
{
    typedef typename Iterator<TTrieIndex, TopDown<> >::Type   TTreeIter;
    typedef typename Iterator<TQuery const, Standard>::Type   TQueryIter;

    TTreeIter   trieIt(index);
    TQueryIter  qryIt  = begin(query, Standard());
    TQueryIter  qryEnd = end(query, Standard());

    _goDownTrie(trieIt, qryIt, qryEnd, maxFreq);
    return range(trieIt);
}

// Gapped Seeds:

template <typename TIndexText, typename TMod, typename TLookupIndex, typename TQuery, typename TSize>
inline Pair<typename Size<Index<TIndexText, IndexSa<Gapped<TMod> > > >::Type>
adaptiveSeeds(Index<TIndexText, IndexSa<Gapped<TMod> > > & index,
              TLookupIndex & table,
              TQuery const & query,
              TSize maxFreq)
{
    typedef Index<TIndexText, IndexSa<Gapped<TMod> > >          TTrieIndex;
    typedef typename Iterator<TTrieIndex, TopDown<> >::Type     TTreeIter;
    typedef ModifiedString<TQuery const, TMod>                  TModStr;
    typedef typename Iterator<TModStr, Standard>::Type          TQueryIter;

    TTreeIter   trieIt(index);
    TModStr     modQuery(query);
    TQueryIter  qryIt  = begin(modQuery, Standard());
    TQueryIter  qryEnd = end(modQuery, Standard());

    _fastTableLookup(trieIt, table, qryIt, qryEnd, maxFreq);
    _goDownTrie(trieIt, qryIt, qryEnd, maxFreq);
    return range(trieIt);
}

template <typename TIndexText, typename TMod, typename TQuery, typename TSize>
inline Pair<typename Size<Index<TIndexText, IndexSa<Gapped<TMod> > > >::Type>
adaptiveSeeds(Index<TIndexText, IndexSa<Gapped<TMod> > > & index,
              TQuery const & query,
              TSize maxFreq)
{
    typedef Index<TIndexText, IndexSa<Gapped<TMod> > >          TTrieIndex;
    typedef typename Iterator<TTrieIndex, TopDown<> >::Type     TTreeIter;
    typedef ModifiedString<TQuery const, TMod>                  TModStr;
    typedef typename Iterator<TModStr, Standard>::Type          TQueryIter;

    TTreeIter   trieIt(index);
    TModStr     modQuery(query);
    TQueryIter  qryIt  = begin(modQuery, Standard());
    TQueryIter  qryEnd = end(modQuery, Standard());

    _goDownTrie(trieIt, qryIt, qryEnd, maxFreq);
    return range(trieIt);
}


// --------------------------------------------------------------------------
// Function _myUngapedExtendSeed()
// --------------------------------------------------------------------------

/*!
 * @fn SeqanLast#_myUngapedExtendSeed
 * @headerfile seqanLast/seqanLast_core.h
 *
 * @brief Ungapped xDrop Extension
 * 
 * write a documentation?
 */

template <typename TConfig,
    typename TDatabase,
    typename TQuery,
    typename TScoreValue,
    typename TScoreSpec>
inline void
_myUngapedExtendSeed(Seed<Simple, TConfig> & seed,
                    TDatabase const & database,
                    TQuery const & query,
                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                    TScoreValue scoreDropOff)
{
    // Horizontal = database
    // Vertical   = query

    typedef typename Position<Seed<Simple, TConfig> >::Type     TPosition;
    typedef typename Size<Seed<Simple, TConfig> >::Type         TSize;
    typedef typename Iterator<TDatabase const, Standard>::Type  TDbIter;
    typedef typename Iterator<TQuery const, Standard>::Type     TQuIter;

    TDbIter dbBeg = begin(database, Standard());
    TDbIter dbEnd = end(database, Standard());
    TDbIter quBeg = begin(query, Standard());
    TDbIter quEnd = end(query, Standard());

    TScoreValue tmpScore, maxScoreLeft, maxScoreRight;
    TPosition len, optLenLeft, optLenRight;
    TDbIter dbIt;
    TQuIter quIt;

   	// Extension to the left
    dbIt = dbBeg + beginPositionH(seed);
    quIt = quBeg + beginPositionV(seed);
    tmpScore = maxScoreLeft = score(seed);
    len = optLenLeft = 0;
    while (dbIt > dbBeg && quIt > quBeg && tmpScore >= maxScoreLeft - scoreDropOff)
    {
        --dbIt; --quIt; ++len;
        tmpScore += score(scoringScheme, *dbIt, *quIt);
        if (tmpScore > maxScoreLeft)
        {
            maxScoreLeft = tmpScore;
            optLenLeft = len;
        }
    }

    // Extension to the right
    dbIt = dbBeg + endPositionH(seed);
    quIt = quBeg + endPositionV(seed);
    tmpScore = maxScoreRight = score(seed);
    len = optLenRight = 0;
    while (dbIt < dbEnd && quIt < quEnd && tmpScore >= maxScoreRight - scoreDropOff)
    {
        tmpScore += score(scoringScheme, *dbIt, *quIt);
        ++dbIt; ++quIt; ++len;
        if (tmpScore > maxScoreRight)
        {
            maxScoreRight = tmpScore;
            optLenRight = len;
        }
    }

    // update seed
    setBeginPositionH(seed, beginPositionH(seed) - optLenLeft);
    setBeginPositionV(seed, beginPositionV(seed) - optLenLeft);
    setEndPositionH(seed, endPositionH(seed) + optLenRight);
    setEndPositionV(seed, endPositionV(seed) + optLenRight);
    setScore(seed, score(seed) + maxScoreLeft + maxScoreRight);
}


// --------------------------------------------------------------------------
// Function _seqanUngapedExtendSeed()
// --------------------------------------------------------------------------

template <typename TConfig,
typename TDatabase,
typename TQuery,
typename TScoreValue,
typename TScoreSpec>
inline void
_seqanUngappedExtendSeed(Seed<Simple, TConfig> & seed,
                    TDatabase const & database,
                    TQuery const & query,
                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                    TScoreValue scoreDropOff)
{
    // Horizontal = database
    // Vertical   = query

    typedef typename Position<Seed<Simple, TConfig> >::Type     TPosition;
    typedef typename Size<Seed<Simple, TConfig> >::Type         TSize;
    typedef typename Iterator<TDatabase const, Standard>::Type  TDbIter;
    typedef typename Iterator<TQuery const, Standard>::Type     TQuIter;

    TDbIter dbBeg = begin(database, Standard());
    TDbIter dbEnd = end(database, Standard());
    TDbIter quBeg = begin(query, Standard());
    TDbIter quEnd = end(query, Standard());

    TScoreValue tmpScore, maxScoreLeft, maxScoreRight;
    TPosition len, optLenLeft, optLenRight;
    TDbIter dbIt;
    TQuIter quIt;

    std::cout << "The Seqan module for the ungapped extension is not fully implemented" << std::endl;
    // TODO
}

// -----------------------------------------------------------------------------
// Function myExtendAlignment()
// -----------------------------------------------------------------------------

template <typename TAlignObject,
    typename TDatabase,
    typename TQuery,
    typename TScoreValue,
    typename TScoreSpec>
inline TScoreValue myExtendAlignment(
                                     TAlignObject &                  alignObj,
                                     TDatabase const &               database,
                                     TQuery const &                  query,
                                     Score<TScoreValue, TScoreSpec> const & scoreMatrix,
                                     TScoreValue                     gappedXDropScore)
{
    typedef typename Size<TDatabase>::Type                     TSize;

    // Run a local alignment first to get the "core" of the alignment
    TScoreValue localScore = localAlignment(alignObj, scoreMatrix);

    /* TODO
    // Now extend both ends
    Tuple<TSize, 4> positions;
    positions[0] = beginPositionH(seed) + beginPosition(row(alignObj, 0));
    positions[1] = beginPositionV(seed) + beginPosition(row(alignObj, 1));
    positions[2] = beginPositionH(seed) + endPosition(row(alignObj, 0));
    positions[3] = beginPositionV(seed) + endPosition(row(alignObj, 1));


    // TODO: extendAlignment mit AliExtContext damit die Matrizen nicht immer wieder allokiert werden müssen!

    TScoreValue finalScore = extendAlignment(alignObj, localScore, database, query, positions,
                                             EXTEND_BOTH,
                                             -25,       // lower Diag           // TODO(meiers): Choose band width
                                             +25,       // upper Diag
                                             gappedXDropScore, scoreMatrix);


    return finalScore;
     */
    return localScore;
}




// TODO:
//      - make sort outside only sort references, not objects
//      - enable hashTable
//      - write a bit of documentation


// -----------------------------------------------------------------------------
// Function _prepareMatchObject()
// -----------------------------------------------------------------------------

template <typename TMatch, typename TSeed, typename TString1, typename TString2, typename TSize1, typename TSize2>
inline void _prepareMatchObject(TMatch & match,
                                TSeed const & seed,
                                TString1 const & database,
                                TString2 const & query,
                                TSize1 const & dbId,
                                TSize2 const & quId)
{
    // create ungapped align Object from the seed
    resize(rows(match.align), 2);
    assignSource(row(match.align, 0), infix(database,   beginPositionH(seed), endPositionH(seed)));
    assignSource(row(match.align, 1), infix(query,      beginPositionV(seed), endPositionV(seed)));

    // Set other parameters of the match
    match.score = score(seed);
    match.dbId = dbId;
    match.quId = quId;
}

// -----------------------------------------------------------------------------
// Function linearLastal()
// -----------------------------------------------------------------------------

template <typename TMatches,
    typename TDbSet,
    typename TIndexSpec,
    typename TShape,
    typename TQuerySet,
    typename TSize2,
    typename TScoreMatrix>
void linearLastal(TMatches                                   & finalMatches,
                  Index<TDbSet, IndexSa<TIndexSpec> >        & index,
                  Index<TDbSet, IndexQGram<TShape> >         & table,
                  TQuerySet                            const & querySet,
                  LastParameters<TSize2, TScoreMatrix> const & params)
{
    typedef Index<TDbSet, IndexSa<TIndexSpec> >                     TIndex;
    typedef typename Size<TIndex>::Type                             TDbSize;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;
    typedef typename Iterator<TSA, Standard>::Type                  TSAIter;
    typedef typename Value<TDbSet const>::Type                      TDbString;
    typedef typename Value<TQuerySet const>::Type                   TQueryString;
    typedef typename Size<TQuerySet>::Type                          TQuSize;
    typedef typename Iterator<TQueryString const, Standard>::Type   TQueryIter;
    typedef DiagonalTable<typename Difference<TDbString>::Type,
                          typename Size<TDbString>::Type>           TDiagTable;
    typedef typename Value<TMatches>::Type                          TMatch;
    typedef typename Value<TScoreMatrix>::Type                      TScore;


    // Self-measurements
    double      _tASCalls = 0;
    unsigned    _cASCalls = 0;
    unsigned    _cSeeds = 0;
    double      _tglAlsCalls = 0;
    unsigned    _cglAls = 0;
    double      _tgpAlsCalls = 0;
    unsigned    _cgpAls = 0;


    resize(globalCounts, 0);

    // diagonal tables for the identification of redundant hits
    String<TDiagTable> diagTables;
    resize(diagTables, length(indexText(index)));

    for(TQuSize queryId=0; queryId < length(querySet); ++queryId)
    {
        TQueryString const &    query = querySet[queryId];
        TQuSize                 queryPos = 0;
        TQueryIter              queryBeg = begin(query, Standard());
        TQueryIter              queryEnd = end(query, Standard());

        // reset diagonal tables
        for(unsigned i=0; i<length(indexText(index)); ++i)
            diagTables[i].clear();

        // Sequential search over queries
        for(TQueryIter queryIt = queryBeg; queryIt != queryEnd; ++queryIt, ++queryPos)
        {
            // Lookup adaptive Seed
            //double xxx = cpuTime();
            Pair<TDbSize> range = (params.useHashTable ?
                                   adaptiveSeeds(index, table, suffix(query, queryPos), params.maxFreq) :
                                   adaptiveSeeds(index, suffix(query, queryPos), params.maxFreq));
            //_tASCalls += cpuTime() - xxx;
            ++_cASCalls;

            if(range.i2 <= range.i1) continue; // seed doesn't hit at all
            if(range.i2 - range.i1 > params.maxFreq) continue; // seed hits too often

            // Enumerate adaptive seeds
            TSAIter saIt = begin(indexSA(index), Standard()) + range.i1;
            TSAIter saEnd  = begin(indexSA(index), Standard()) + range.i2;

            for(; saIt != saEnd; ++saIt)
            {
                ++_cSeeds;
                TDiagTable        & diagTable = diagTables[getSeqNo(*saIt)];
                TDbString const   & database  = indexText(index)[getSeqNo(*saIt)];
                Seed<Simple>        seed( getSeqOffset(*saIt), queryIt - queryBeg, 0 );

                // Check whether seed is redundant
                if (diagTable.redundant(beginPositionH(seed), beginPositionV(seed)))
                    continue;

                // Gapless Alignment in both directions with a XDrop
                //double xxxx = cpuTime();
                if (params.myUngappedExtend)
                    _myUngapedExtendSeed(seed, database, query, params.scoreMatrix, params.Xgapless);
                else
                    _seqanUngappedExtendSeed(seed, database, query, params.scoreMatrix, params.Xgapless);
                //_tglAlsCalls += cpuTime() - xxxx;
                ++_cglAls;

                // gapLess alignment too weak
                if (score(seed) < params.Tgapless) continue;

                // Mark diagonal as already visited
                diagTable.add(endPositionH(seed), endPositionV(seed));

                // Prepare a match object with an align object inside
                TMatch matchObj;
                _prepareMatchObject(matchObj, seed, database, query, getSeqNo(*saIt), queryId);

                // If only ungapped matches are wanted
                if (params.onlyUngappedAlignments)
                {
                    appendValue(finalMatches, matchObj);
                    continue;
                }

                // Gapped alignment:
                //double xxxxx = cpuTime();
                TScore finalScore = myExtendAlignment(matchObj.align, database, query, params.scoreMatrix, params.Xgapped);
                //_tgpAlsCalls += cpuTime() - xxxxx; ++_cgpAls;

                if (finalScore > params.Tgapped) appendValue(finalMatches, matchObj);
                
            } // for(; saIt != saEnd; ++saIt)
        } //for(; queryIt != queryEnd; ++queryIt)
    }


    // test: count vector of adaptive seeds
    std::cout << "Counts:" << std::endl;
    for (unsigned x=0; x<length(globalCounts); ++x)
        std::cout << x << "\t" << globalCounts[x] << std::endl;


    if (params.verbosity > 1)
    {
        std::cout << "Time spend in adaptive seeding:  " << _tASCalls <<    "\t(" << _cASCalls << " calls)" << std::endl;
        std::cout << "Time spend in gapless alignment: " << _tglAlsCalls << "\t(" << _cglAls <<   " calls)" << std::endl;
        std::cout << "Time spend in gapped alignment:  " << _tgpAlsCalls << "\t(" << _cgpAls <<   " calls)" << std::endl;
        std::cout << std::endl;
        std::cout << " # adaptive seeds:     " << _cSeeds << std::endl;
    }
    if (params.verbosity)
        std::cout << " # final matches:      " << length(finalMatches) << std::endl;
}




#endif  // #ifndef SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_CORE_H_
