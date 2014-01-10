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
// Forwards
// =============================================================================

// =============================================================================
// Tags, Classes, Enums
// =============================================================================

// -----------------------------------------------------------------------------
// Struct DiagonalTable
// -----------------------------------------------------------------------------

template <typename TSigned, typename TSize, unsigned Q=64>
struct DiagonalTable
{
    typedef typename Iterator<String<Pair<TSigned,TSize> > >::Type TIter;
    String<Pair<TSigned, TSize> > table[Q];

    DiagonalTable()
    {
        for(unsigned i=0; i<Q; ++i)
            table[i] = String<Pair<TSigned, TSize> >();
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

template <typename TSize, typename TAlign>
struct SeqanLastMatch
{
    typedef TAlign Type;
    TAlign align;
    TSize  dbId;
    TSize  quId;
};

// =============================================================================
// Metafunctions
// =============================================================================

// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
// Function adaptiveSeeds()
// -----------------------------------------------------------------------------

template <typename TTrieIndex,
typename TLookupIndex,
typename TQuery,
typename TSize2 >
inline Pair<typename Size<TTrieIndex>::Type>
adaptiveSeeds(TTrieIndex &index, TLookupIndex & table,
              TQuery const &query, TSize2 maxFreq, int verbosity,
              TSize2 goDownAtOnce = 0)
{
    // NOTE: query and database modifiers are expected to match
    //       (i.e. be of the same subset/shape type)
    // NOTE: Does only return the SA range, not the seed length

    typedef typename Iterator<TTrieIndex, TopDown<> >::Type   TTreeIter;
    typedef typename Size<TTrieIndex>::Type                   TSize;
    typedef typename Iterator<TQuery const, Standard>::Type   TQueryIter;
    typedef typename Fibre<TLookupIndex, QGramShape>::Type    TShape;
    // TODO: Here Shape is always ungapped, but in the index it can be gapped !!
    typedef typename Value<TShape>::Type                      THash;
    typedef typename Host<TShape>::Type                       TShapeAlph;

    TShape qryHash;
    TQueryIter qry = begin(query, Standard());
    TQueryIter qryEnd = end(query, Standard());
    TTreeIter treeIter(index); // root

    // 1.
    // Lookup the hash table only for full length Shapes to avoid
    // possible problems with Open Adressing indiex.
    if(static_cast<TSize>(qryEnd - qry) >= weight(qryHash) )
    {
        hash(qryHash, qry);
        Pair<TSize> initialRange = range(table, qryHash);

        if(initialRange.i2 - initialRange.i1 > maxFreq)
        {
            treeIter.vDesc.range = initialRange;
            treeIter.vDesc.repLen = weight(qryHash);
            qry += weight(qryHash);
            treeIter.vDesc.lastChar = *qry;
        }
    }

    // 2.
    // Continue or restart search in the Suffix array.
    while(qry < qryEnd)
    {
        if(!goDown(treeIter, *(qry++)))
            break;

        if(countOccurrences(treeIter) <= maxFreq)
            break;
    }

    return range(treeIter);
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
    // TODO: Use Iterators instead of positional access

    // Horizontal = database
    // Vertical   = query

    typedef Seed<ChainedSeed, TConfig>      TSeed;
    typedef typename Position<TSeed>::Type  TPosition;
    typedef typename Size<TSeed>::Type      TSize;

	// Extension to the left
    TScoreValue tmpScore = score(seed);
    TScoreValue maxScore = score(seed);
    TPosition posH = beginPositionH(seed);
    TPosition posV = beginPositionV(seed);

    for (; posH >= 1 && posV >= 1 && tmpScore > maxScore - scoreDropOff; --posH, --posV)
    {
        tmpScore += score(scoringScheme, database[posH -1], query[posV -1]);
        if (tmpScore > maxScore)
        {
            maxScore = tmpScore;
            setBeginPositionH(seed, posH -1);
            setBeginPositionV(seed, posV -1);
        }
    }
    setScore(seed, maxScore);

    // Extension to the right
    tmpScore = score(seed);
    maxScore = score(seed);
    posH = endPositionH(seed);
    posV = endPositionV(seed);
    TSize lengthH = length(database);
    TSize lengthV = length(query);

    for (; posH < lengthH && posV < lengthV && tmpScore > maxScore - scoreDropOff; ++posH, ++posV)
    {
        tmpScore += score(scoringScheme, database[posH], query[posV]);
        if (tmpScore > maxScore)
        {
            maxScore = tmpScore;
            setEndPositionH(seed, posH + 1);
            setEndPositionV(seed, posV + 1);
        }

    }
    setScore(seed, maxScore);
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

    TScoreValue finalScore = extendAlignment(alignObj, localScore, database, query, positions, EXTEND_BOTH, gappedXDropScore, scoreMatrix);

    return finalScore;
}

// -----------------------------------------------------------------------------
// Function linearLastal()
// -----------------------------------------------------------------------------

template <typename TMatches,
    typename TDatabase,
    typename TDatabaseSetSpec,
    typename TIndexSpec,
    typename TShape,
    typename TQueryString,
    typename TQuerySetSpec,
    typename TSize2,
    typename TScoreMatrix,
    typename TScore>
void linearLastal(
                  TMatches & finalMatches,
                  Index<StringSet<TDatabase, TDatabaseSetSpec>, IndexSa<TIndexSpec> > &index,
                  Index<StringSet<TDatabase, TDatabaseSetSpec>, IndexQGram<TShape> > &table,
                  StringSet<TQueryString, TQuerySetSpec> const &querySet,
                  TSize2 maxFreq,
                  TScoreMatrix const & scoreMatrix,
                  TScore glXdrop,
                  TScore gpXdrop,
                  TScore glThr,
                  TScore gpThr,
                  int verbosity)
{
    typedef StringSet<TDatabase, TDatabaseSetSpec>                  TDatabaseSet;
    typedef Index<TDatabaseSet, IndexSa<TIndexSpec> >               TIndex;
    typedef typename Size<TIndex>::Type                             TSize;
    typedef StringSet<TQueryString, TQuerySetSpec> const            TQuerySet;
    typedef typename Iterator<TQuerySet, Standard>::Type            TQuerySetIter;
    typedef typename Iterator<TQueryString const, Standard>::Type   TQueryIter;
    typedef String<typename SAValue<TIndex>::Type>                  TSA;
    typedef typename Iterator<TSA, Standard>::Type                  TSAIter;
    typedef DiagonalTable<typename Difference<TDatabase>::Type, typename Size<TDatabase>::Type> TDiagTable;
    typedef typename Value<TMatches>::Type TMatch;

    // Self-measurements
    double      _tASCalls = 0;
    unsigned    _cASCalls = 0;
    unsigned    _cSeeds = 0;
    unsigned    _cgpAls = 0;
    unsigned    _cglAls = 0;

    TSize L = length(indexText(index));
    String<TDiagTable> diagTables;
    resize(diagTables, L * length(querySet));

    // Step 0:
    // Check that SA exists
    indexRequire(index, FibreSA());
    indexRequire(table, FibreDir());

    // Linear search on query
    for(typename Size<TQuerySet>::Type queryId=0; queryId < length(querySet); ++queryId)
    {
        TQueryString const & query = querySet[queryId];
        typename Size<TQueryString>::Type queryPos = 0;
        TQueryIter queryBeg = begin(query, Standard());
        TQueryIter queryEnd = end(query, Standard());

        for(TQueryIter queryIt = queryBeg; queryIt != queryEnd; ++queryIt, ++queryPos)
        {

            // Step 1:
            // Lookup adaptive Seed
            double xxx = cpuTime();
            Pair<TSize> range = adaptiveSeeds(index, table, suffix(query, queryPos), maxFreq, verbosity);
            _tASCalls += cpuTime() - xxx;
            ++_cASCalls;

            if (verbosity > 2)
                std::cout << "query [" << queryId << "," << queryPos << "] : \t\"" <<
                infix(query, queryPos, std::min(static_cast<TSize>(queryPos + 30),
                                                static_cast<TSize>(queryEnd-queryBeg))) << "...\", SA range: " << range <<
                (range.i2 > range.i1 && range.i2 <= range.i1 + maxFreq ? " good" : " BAD") << std::endl;

            if(range.i2 <= range.i1) continue; // seed doesn't hit at all
            if(range.i2 - range.i1 > maxFreq) continue; // seed hits too often

            TSAIter saFrom = begin(indexSA(index), Standard()) + range.i1;
            TSAIter saEnd  = begin(indexSA(index), Standard()) + range.i2;

            for(; saFrom != saEnd; ++saFrom)
            {
                ++_cSeeds;

                // Step 2:
                // Gapless Alignment in both directions with a XDrop
                TDatabase const & database = indexText(index)[getSeqNo(*saFrom)];
                Seed<Simple> seed(getSeqOffset(*saFrom), queryIt - queryBeg, 0);

                // Check whether seed is redundant
                TDiagTable &diagTable = diagTables[getSeqNo(*saFrom) * L + queryId];
                if (diagTable.redundant(beginPositionH(seed), beginPositionV(seed)))
                    continue;

                if (verbosity > 2)
                    std::cout << "      seed adaptive   database [" << getSeqNo(*saFrom) << "," << beginPositionH(seed) <<
                    "-" << endPositionH(seed) << "]\tquery [" << queryId << "," << beginPositionV(seed)
                    << "-" << endPositionV(seed) << "]"   << std::endl;

                myUngapedExtendSeed(seed, database, query, scoreMatrix, glXdrop);
                diagTable.add(endPositionH(seed), endPositionV(seed));
                ++_cglAls;

                if (verbosity > 2)
                    std::cout << "           extended   database [" << getSeqNo(*saFrom) << "," << beginPositionH(seed) <<
                    "-" << endPositionH(seed) << "]\tquery [" << queryPos << "," << beginPositionV(seed)
                    << "-" << endPositionV(seed) << "]\tscore: " << score(seed) << std::endl;

                // gapLess alignment too weak
                if (score(seed) < glThr) continue;

                // Step 3:
                // Gapped alignment:
                typename TMatch::Type alignObj;

                TScore finalScore = myExtendAlignment(alignObj, seed, database, query, scoreMatrix, gpXdrop);
                ++_cgpAls;

                if (verbosity > 2)
                    std::cout << "      * aligned    database [" << getSeqNo(*saFrom) << "," <<
                    beginPosition(row(alignObj,0)) << "-" << endPosition(row(alignObj,0)) << "]\tquery [" <<
                    queryPos << "," << beginPosition(row(alignObj,1)) << "-" <<
                    endPosition(row(alignObj,1)) << "]\tscore: " << finalScore << std::endl;

                if (finalScore > gpThr)
                {
                    TMatch matchObj;
                    matchObj.quId = queryId;
                    matchObj.dbId = getSeqNo(*saFrom);
                    matchObj.align = alignObj;
                    appendValue(finalMatches, matchObj);
                }
                
            } // for(; saFrom != saEnd; ++saFrom)
            
        } //for(; queryIt != queryEnd; ++queryIt)
    }
    if (verbosity > 1)
        std::cout << "Time spend in adaptive seeds: " << _tASCalls << "/" << _cASCalls << " calls" << std::endl;
    if (verbosity > 1)
    {
        std::cout << " # adaptive seeds:     " << _cSeeds << std::endl;
        std::cout << " # gapless alignments: " << _cglAls << std::endl;
        std::cout << " # gapped alignments:  " << _cgpAls << std::endl;
        std::cout << " # final matches:      " << length(finalMatches) << std::endl;
    }
}




#endif  // #ifndef SANDBOX_MEIERS_APPS_SEQANLAST_SEQANLAST_CORE_H_
