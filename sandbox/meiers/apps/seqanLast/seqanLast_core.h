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


// only output the ungapped extension!
#define SEQANLAST_ONLY_UNGAPPED_EXTENSION

// =============================================================================
// Global Definitions
// =============================================================================

// StringSet
typedef StringSet<String<Dna5, MMap<> >, Owner<ConcatDirect<> > >         TMMapStringSet;
typedef StringSet<String<Dna5>,          Owner<ConcatDirect<> > >         TNormalStringSet;

// Shapes
typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1> >, 1> >       Shape1;   // 110
typedef CyclicShape<FixedShape<0, GappedShape<HardwiredShape<1,1> >, 1> >     Shape2;   // 1110
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
    }

    bool redundant(TSize pGenome, TSize pQuery)
    {
        TSigned d = static_cast<TSigned>(pGenome) - static_cast<TSigned>(pQuery);
        TSigned dd =(Q + d%Q)%Q; // hope the modulo is replaced by fast bit operations.

        for(TIter it=begin(table[dd]); it != end(table[dd]); ++it)
        {
            if(it->i1 != d) continue;
            if(it->i2 > pQuery) return true;
        }
        return false;
    }

    void add(TSize pGenome, TSize pQuery)
    {
        TSigned d = static_cast<TSigned>(pGenome) - static_cast<TSigned>(pQuery);
        TSigned dd =(Q + d%Q)%Q; // hope the modulo is replaced by fast bit operations.
        String<Pair<TSigned, TSize> > newList;
        for(TIter it=begin(table[dd]); it != end(table[dd]); ++it)
        {
            if(it->i1 != d)
                appendValue(newList, *it);
            else
                if(it->i2 > pQuery) // there is already one further right
                    pQuery = it->i2;
        }
        appendValue(newList, Pair<TSigned,TSize>(d,pQuery));
        table[dd] = newList;
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
    bool onlyUngappedAlignments;
    int verbosity;

    LastParameters(TSize f, TScoreMatrix const & scoreMatrix, TScore Xgapless,
                   TScore Xgapped, TScore Tgapless, TScore Tgapped, bool ungAl, int v) :
        maxFreq(f), scoreMatrix(scoreMatrix), Xgapless(Xgapless), Xgapped(Xgapped),
        Tgapless(Tgapless), Tgapped(Tgapped), onlyUngappedAlignments(ungAl),
        verbosity(v)
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
// Function _goDownTrie()
// -----------------------------------------------------------------------------

template <typename TTrieIt, typename TLookupIndex, typename TQueryIt, typename TSize>
inline void _goDownTrie(TTrieIt & trieIt,
                        TLookupIndex  & table,
                        TQueryIt qryIt,
                        TQueryIt qryEnd,
                        TSize maxFreq)
{

    typedef typename Fibre<TLookupIndex, FibreShape>::Type TShape;
    typedef typename Value<TShape>::Type THashValue;
    typedef typename Size<TLookupIndex>::Type TSaPos;

    // determine hash values
    THashValue x = hash(indexShape(table), qryIt, qryEnd - qryIt);
    THashValue y = x+1;
    if (static_cast<TSaPos>(qryEnd - qryIt) < length(indexShape(table)))
        y = hashUpper(indexShape(table), qryIt, qryEnd - qryIt);

    // get range in SA
    TSaPos from  = indexDir(table)[x];
    TSaPos to    = indexDir(table)[y];

    // EITHER: make seed longer
    if ( to - from > maxFreq)
    {
        value(trieIt).range.i1 = from;
        value(trieIt).range.i2 = to;
        value(trieIt).repLen = weight(indexShape(table));
        goFurther(qryIt, weight(indexShape(table)) - 1);
        value(trieIt).lastChar = *qryIt++;
    }

    // OR: make seed shorter
    else
    {
        std::cout << "make shorter" << std::endl;
    }

    while(qryIt < qryEnd)
    {
        if(!goDown(trieIt, *(qryIt++)))
            break;
        if(countOccurrences(trieIt) <= maxFreq)
            break;
    }


}

// -----------------------------------------------------------------------------
// Function adaptiveSeeds()                                           [ungapped]
// -----------------------------------------------------------------------------

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

    _goDownTrie(trieIt, table, qryIt, qryEnd, maxFreq);
    return range(trieIt);
}

// -----------------------------------------------------------------------------
// Function adaptiveSeeds()                                             [Gapped]
// -----------------------------------------------------------------------------

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

    _goDownTrie(trieIt, table, qryIt, qryEnd, maxFreq);
    return range(trieIt);
}



// --------------------------------------------------------------------------
// Function myUngapedExtendSeed()
// --------------------------------------------------------------------------

template <typename TConfig,
    typename TDatabase,
    typename TQuery,
    typename TScoreValue,
    typename TScoreSpec>
inline void
myUngapedExtendSeed(Seed<Simple, TConfig> & seed,
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
    while (dbIt > dbBeg && quIt > quBeg && tmpScore > maxScoreLeft - scoreDropOff)
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
    while (dbIt < dbEnd && quIt < quEnd && tmpScore > maxScoreRight - scoreDropOff)
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


// -----------------------------------------------------------------------------
// Function myExtendAlignment()
// -----------------------------------------------------------------------------

template <
    typename TAlignObject,
    typename TConfig,
    typename TDatabase,
    typename TQuery,
    typename TScoreValue,
    typename TScoreSpec>
inline TScoreValue myExtendAlignment(
                                     TAlignObject &                  alignObj,
                                     Seed<Simple, TConfig> const &   seed,
                                     TDatabase const &               database,
                                     TQuery const &                  query,
                                     Score<TScoreValue, TScoreSpec> const & scoreMatrix,
                                     TScoreValue                     gappedXDropScore)
{
    typedef typename Size<TDatabase>::Type                     TSize;
    
    resize(rows(alignObj), 2);
    assignSource(row(alignObj, 0), infix(database, beginPositionH(seed), endPositionH(seed)));
    assignSource(row(alignObj, 1), infix(query, beginPositionV(seed), endPositionV(seed)));

    // Run a local alignment first to get the "core" of the alignment
    TScoreValue localScore = localAlignment(alignObj, scoreMatrix);

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
void linearLastal(TMatches & finalMatches,
                  Index<TDbSet, IndexSa<TIndexSpec> > & index,
                  Index<TDbSet, IndexQGram<TShape> >  & table,
                  TQuerySet      const & querySet,
                  LastParameters<TSize2, TScoreMatrix> const & params)
{
    typedef Index<TDbSet, IndexSa<TIndexSpec> >                     TIndex;
    typedef typename Size<TIndex>::Type                             TDbSize;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;
    typedef typename Iterator<TSA, Standard>::Type                  TSAIter;

    typedef typename Value<TDbSet const>::Type                      TDbString;

    typedef typename Value<TQuerySet const>::Type                   TQueryString;
    typedef typename Size<TQuerySet>::Type                          TQuSize;
    typedef typename Iterator<TQuerySet const, Standard>::Type      TQuerySetIter;
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


    // TODO(meiers): Only need |db| many diag tables, since only one query is passed at a time !!
    TDbSize L = length(indexText(index));
    String<TDiagTable> diagTables;
    resize(diagTables, L * length(querySet));

    // Sequential search over queries
    for(TQuSize queryId=0; queryId < length(querySet); ++queryId)
    {
        TQueryString const &    query = querySet[queryId];
        TQuSize                 queryPos = 0;
        TQueryIter              queryBeg = begin(query, Standard());
        TQueryIter              queryEnd = end(query, Standard());

        for(TQueryIter queryIt = queryBeg; queryIt != queryEnd; ++queryIt, ++queryPos)
        {
            // Lookup adaptive Seed
                    double xxx = cpuTime();
            Pair<TDbSize> range = adaptiveSeeds(index, table, suffix(query, queryPos), params.maxFreq);
                    _tASCalls += cpuTime() - xxx;
                    ++_cASCalls;

            if(range.i2 <= range.i1) continue; // seed doesn't hit at all
            if(range.i2 - range.i1 > params.maxFreq) continue; // seed hits too often

            // Enumerate adaptive seeds
            TSAIter saFrom = begin(indexSA(index), Standard()) + range.i1;
            TSAIter saEnd  = begin(indexSA(index), Standard()) + range.i2;

            for(; saFrom != saEnd; ++saFrom)
            {
                ++_cSeeds;
                TDiagTable          &diagTable = diagTables[getSeqNo(*saFrom) + L*queryId];
                TDbString const     &database  = indexText(index)[getSeqNo(*saFrom)];
                Seed<Simple>        seed( getSeqOffset(*saFrom), queryIt - queryBeg, 0 );

                // Check whether seed is redundant
                if (diagTable.redundant(beginPositionH(seed), beginPositionV(seed)))
                    continue;

                // Gapless Alignment in both directions with a XDrop
                        double xxxx = cpuTime();
                myUngapedExtendSeed(seed, database, query, params.scoreMatrix, params.Xgapless);
                        _tglAlsCalls += cpuTime() - xxxx;
                        ++_cglAls;


                // Mark diagonal as already
                diagTable.add(endPositionH(seed), endPositionV(seed));

                // gapLess alignment too weak
                if (score(seed) < params.Tgapless) continue;



                typename TMatch::Type alignObj;

                // TODO: make a switch for ungapped extension based on params
                //      - wrap function to prepare matchObject
                //      - handle diagTables only per query, not M x N many
                //      - make sort outside only sort references, not objects
                //      - use iterators in ungapped Extension
                //      - make a switch between my ungapped extension and the seqan version
                //      - enable hashTable
                //      - write a bit of documentation
#ifdef SEQANLAST_ONLY_UNGAPPED_EXTENSION
                resize(rows(alignObj), 2);
                assignSource(row(alignObj, 0), infix(database, beginPositionH(seed), endPositionH(seed)));
                assignSource(row(alignObj, 1), infix(query, beginPositionV(seed), endPositionV(seed)));
                TMatch matchObj;
                matchObj.quId = queryId;
                matchObj.dbId = getSeqNo(*saFrom);
                matchObj.align = alignObj;
                matchObj.score = score(seed);
                appendValue(finalMatches, matchObj);
                continue;
#endif

                // Gapped alignment:

                        double xxxxx = cpuTime();
                TScore finalScore = myExtendAlignment(alignObj, seed, database, query, params.scoreMatrix, params.Xgapped);
                        _tgpAlsCalls += cpuTime() - xxxxx;
                        ++_cgpAls;

                if (finalScore > params.Tgapped)
                {
                    TMatch matchObj;
                    matchObj.quId = queryId;
                    matchObj.dbId = getSeqNo(*saFrom);
                    matchObj.align = alignObj;
                    matchObj.score = finalScore;
                    appendValue(finalMatches, matchObj);
                }
                
            } // for(; saFrom != saEnd; ++saFrom)
            

        } //for(; queryIt != queryEnd; ++queryIt)
    }

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
