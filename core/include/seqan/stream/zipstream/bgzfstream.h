/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003   (original zlib stream)
Author: David Weese, dave.weese@gmail.com, 2014             (extension to parallel block-wise compression in bgzf format)
*/

#ifndef BGZFSTREAM_HPP
#define BGZFSTREAM_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include "zutil.h"

namespace seqan {

/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

const unsigned BGZF_MAX_BLOCK_SIZE = 64 * 1024;
const unsigned BGZF_BLOCK_HEADER_LENGTH = 18;
const unsigned BGZF_BLOCK_FOOTER_LENGTH = 8;
const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

// Reduce the maximal input size, such that the compressed data 
// always fits in one block even for level Z_NO_COMPRESSION.
const unsigned BGZF_BLOCK_SIZE = BGZF_MAX_BLOCK_SIZE - BGZF_BLOCK_HEADER_LENGTH - BGZF_BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;

/// Compression strategy, see bgzf doc.
enum EStrategy
{
	StrategyFiltered = 1,
	StrategyHuffmanOnly = 2,
	DefaultStrategy = 0
};


// new Event class

struct Event2
{
    enum { Infinite = LONG_MAX };

    pthread_cond_t  data_cond;
    Mutex &         mutex;

    explicit
    Event2(Mutex &mutex) :
        mutex(mutex)
    {
        int result = pthread_cond_init(&data_cond, NULL);
        ignoreUnusedVariableWarning(result);
        SEQAN_ASSERT_EQ(result, 0);
    }

    ~Event2()
    {
        int result = pthread_cond_destroy(&data_cond);
        ignoreUnusedVariableWarning(result);
        SEQAN_ASSERT_EQ(result, 0);
    }
};

inline void
waitFor(Event2 &event)
{
    int result = pthread_cond_wait(&event.data_cond, event.mutex.hMutex);
    ignoreUnusedVariableWarning(result);
    SEQAN_ASSERT_EQ(result, 0);
}

inline void
waitFor(Event2 &event, long timeoutMilliSec, bool & inProgress)
{
    if (timeoutMilliSec != Event2::Infinite)
    {
        timespec ts;
        ts.tv_sec = timeoutMilliSec / 1000;
        ts.tv_nsec = (timeoutMilliSec % 1000) * 1000;
        int result = pthread_cond_timedwait(&event.data_cond, event.mutex.hMutex, &ts);
        inProgress = (result == ETIMEDOUT);
        SEQAN_ASSERT(result == 0 || inProgress);
    }
    else
    {
        inProgress = false;
        waitFor(event);
    }
}

inline void
signal(Event2 &event)
{
    int result = pthread_cond_broadcast(&event.data_cond);
    ignoreUnusedVariableWarning(result);
    SEQAN_ASSERT_EQ(result, 0);
}



template <typename TSpec = void>
struct Suspendable;

template <typename TValue, typename TSpec>
class ConcurrentQueue<TValue, Suspendable<TSpec> >
{
public:
    typedef typename Host<ConcurrentQueue>::Type    TString;
    typedef typename Size<TString>::Type            TSize;

    size_t      readerCount;
    size_t      writerCount;

    TString     data;
    TSize       occupied;
    TSize       nextIn;
    TSize       nextOut;

    Mutex       mutex;
    Event2      more;

    bool        virgin;

    ConcurrentQueue():
        readerCount(0),
        writerCount(0),
        occupied(0),
        nextIn(0),
        nextOut(0),
        mutex(false),
        more(mutex),
        virgin(true)
    {
    }

    ~ConcurrentQueue()
    {
        TValue tmp;

        // wait for all pending readers to finish
        while (readerCount != 0u)
        {}

        while (popFront(tmp, *this))
        {}

        _setLength(data, 0);
    }
};

template <typename TValue>
class ConcurrentQueue<TValue, Suspendable<Limit> >:
    public ConcurrentQueue<TValue, Suspendable<> >
{
public:
    typedef ConcurrentQueue<TValue, Suspendable<> > TBase;
    typedef typename Host<ConcurrentQueue>::Type    TString;
    typedef typename Size<TString>::Type            TSize;

    Event2      less;

    ConcurrentQueue(TSize maxSize):
        TBase(),
        less(this->mutex)
    {
        reserve(this->data, maxSize, Exact());
    }
};

template <typename TValue, typename TSpec>
inline void
lockReading(ConcurrentQueue<TValue, Suspendable<TSpec> > &)
{}

template <typename TValue, typename TSpec>
inline void
unlockReading(ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    if (--me.readerCount == 0u)
        signal(me.less);
}

template <typename TValue, typename TSpec>
inline void
lockWriting(ConcurrentQueue<TValue, Suspendable<TSpec> > &)
{}

template <typename TValue, typename TSpec>
inline void
unlockWriting(ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    if (--me.writerCount == 0u)
        signal(me.more);
}

template <typename TValue, typename TSize, typename TSpec>
inline void
setReaderCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize readerCount)
{
    ScopedLock<Mutex> lock(me.mutex);
    me.readerCount = readerCount;
}

template <typename TValue, typename TSize, typename TSpec>
inline void
setWriterCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize writerCount)
{
    ScopedLock<Mutex> lock(me.mutex);
    me.writerCount = writerCount;
}

template <typename TValue, typename TSize1, typename TSize2, typename TSpec>
inline void
setReaderWriterCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize1 readerCount, TSize2 writerCount)
{
    ScopedLock<Mutex> lock(me.mutex);
    me.readerCount = readerCount;
    me.writerCount = writerCount;
}

template <typename TValue, typename TSpec>
inline bool
_popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;
    typedef typename Iterator<TString, Standard>::Type      TIter;

    TSize cap = capacity(me.data);

    while (me.occupied == 0u && me.writerCount > 0u)
        waitFor(me.more);

    if (me.occupied == 0u)
        return false;

    SEQAN_ASSERT_NEQ(me.occupied, 0u);

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + me.nextOut;
    std::swap(result, *it);
    valueDestruct(it);

    me.nextOut = (me.nextOut + 1) % cap;
    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    return true;
}

template <typename TValue, typename TSpec>
inline bool
_popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;
    typedef typename Iterator<TString, Standard>::Type      TIter;

    TSize cap = capacity(me.data);

    while (me.occupied == 0u && me.writerCount > 0u)
        waitFor(me.more);

    if (me.occupied == 0u)
        return false;

    SEQAN_ASSERT_NEQ(me.occupied, 0u);

    me.nextIn = (me.nextIn + cap - 1) % cap;

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + me.nextIn;
    std::swap(result, *it);
    valueDestruct(it);

    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    return true;
}

template <typename TValue, typename TSpec>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    return _popFront(result, me);
}

template <typename TValue>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<Limit> > & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    if (_popFront(result, me))
    {
        signal(me.less);
        return true;
    }
    return false;
}

template <typename TValue, typename TSpec>
inline bool
popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    return _popBack(result, me);
}

template <typename TValue>
inline bool
popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<Limit> > & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    if (_popBack(result, me))
    {
        signal(me.less);
        return true;
    }
    return false;
}


template <typename TValue, typename TValue2, typename TSpec>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
            TValue2 SEQAN_FORWARD_CARG val)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;

    ScopedLock<Mutex> lock(me.mutex);
    TSize cap = capacity(me.data);

    if (me.occupied >= cap)
    {
        if (me.nextOut >= me.nextIn)
            ++me.nextOut;
        
        _setLength(me.data, cap);
        insertValue(me.data, me.nextIn, val);
        ++me.nextIn;
    }
    else
    {
        valueConstruct(begin(me.data, Standard()) + me.nextIn, val);
        me.nextIn = (me.nextIn + 1) % cap;
    }

    me.occupied++;

    /* now: either me.occupied < BSIZE and me.nextin is the index
       of the next empty slot in the buffer, or
       me.occupied == BSIZE and me.nextin is the index of the
       next (occupied) slot that will be emptied by a consumer
       (such as me.nextin == me.nextout) */

    signal(me.more);
    return true;
}

template <typename TValue, typename TValue2>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<Limit> > & me,
            TValue2 SEQAN_FORWARD_CARG val)
{
    typedef ConcurrentQueue<TValue, Suspendable<Limit> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;

    ScopedLock<Mutex> lock(me.mutex);
    TSize cap = capacity(me.data);

    while (me.occupied >= cap && me.readerCount > 0u)
        waitFor(me.less);

    if (me.occupied >= cap)
        return false;

    SEQAN_ASSERT_LT(me.occupied, cap);

    valueConstruct(begin(me.data, Standard()) + me.nextIn, val);
    me.nextIn = (me.nextIn + 1) % cap;
    me.occupied++;

    /* now: either me.occupied < BSIZE and me.nextin is the index
       of the next empty slot in the buffer, or
       me.occupied == BSIZE and me.nextin is the index of the
       next (occupied) slot that will be emptied by a consumer
       (such as me.nextin == me.nextout) */

    signal(me.more);
    return true;
}

template <typename TValue, typename TSize, typename TSpec>
inline bool
waitForMinSize(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
               TSize minSize)
{
    ScopedLock<Mutex> lock(me.mutex);
    while (me.occupied < minSize && me.writerCount > 0u)
        waitFor(me.more);
    return me.occupied >= minSize;
}




template <typename TValue>
struct ResourcePool
{
    typedef ConcurrentQueue<TValue *, Suspendable<> >   TStack;
    typedef typename Size<TStack>::Type                 TSize;

    TStack recycled;

    ResourcePool(TSize maxSize)
    {
        setWriterCount(recycled, 1);
        for (; maxSize != 0; --maxSize)
            appendValue(recycled, (TValue *)NULL);
    }

    ~ResourcePool()
    {
        unlockWriting(recycled);
        TValue *ptr;
        unsigned count =0;
        while (popBack(ptr, recycled))
        {
            if (ptr != NULL)
            count++;
            delete ptr;
        }
        std::cerr<<"RESOURCES:"<<count<<std::endl;
    }
};

template <typename TValue>
inline TValue *
aquireValue(ResourcePool<TValue> & me)
{
    TValue *ptr;
    if (!popBack(ptr, me.recycled))
        return NULL;

    if (ptr == NULL)
        ptr = new TValue();
    return ptr;
}

template <typename TValue>
inline void
releaseValue(ResourcePool<TValue> & me, TValue *ptr)
{
    appendValue(me.recycled, ptr);
}



template <typename TValue>
struct SerializerItem
{
    TValue          val;
    SerializerItem  *next;
    bool            ready;
};

template <typename TValue, typename TWorker>
struct Serializer
{
    typedef SerializerItem<TValue>  TItem;
    typedef TItem *                 TItemPtr;

    Mutex               mutex;
    TWorker             worker;
    TItemPtr            first;
    TItemPtr            last;
    ResourcePool<TItem> pool;
    bool                stop;

    Serializer():
        mutex(false),
        first(NULL),
        last(NULL),
        pool(1024),
        stop(false)
    {
        clear(*this);
    }

    template <typename TArg>
    Serializer(TArg &arg):
        mutex(false),
        worker(arg),
        first(NULL),
        last(NULL),
        pool(1024),
        stop(false)
    {}

    ~Serializer()
    {
        while (first != NULL)
        {
            TItemPtr item = first;
            first = first->next;
            delete item;
            std::cout << "free"<<std::endl;
        }
        std::cerr<<"Serializer:"<<capacity(pool.recycled.data)<<std::endl;
    }

    operator bool()
    {
        return !stop;
    }
};

template <typename TValue, typename TWorker>
inline void
clear(Serializer<TValue, TWorker> & me)
{
    me.stop = false;
    while (me.first != NULL)
    {
        TValue *item = me.first;
        me.first = me.first->next;
        releaseValue(me.recycled, item);
    }
    me.last = NULL;
}

template <typename TValue, typename TWorker>
inline TValue *
aquireValue(Serializer<TValue, TWorker> & me)
{
    typedef SerializerItem<TValue> TItem;

    TItem *item = aquireValue(me.pool);
    item->next = NULL;
    item->ready = false;

    // add item to the end of our linked list
    {
        ScopedLock<Mutex> lock(me.mutex);

        if (me.first != NULL)
        {
            SEQAN_ASSERT(me.last != NULL);
            me.last->next = item;
        }
        else
        {
            me.first = item;
            me.last = item;
        }
    }

//    std::cerr<<"aquire\t"<<(size_t)item<<std::endl;
    return &item->val;
}

template <typename TValue, typename TWorker>
inline bool
releaseValue(Serializer<TValue, TWorker> & me, TValue *ptr)
{
    typedef SerializerItem<TValue> TItem;

    TItem *item = reinterpret_cast<TItem *>(ptr);
    SEQAN_ASSERT_NOT(item->ready);

//    {
//            ScopedLock<Mutex> lock(me.mutex);
//            // recycle released items
//            appendValue(me.recycled, item);
//    }
//    {
//        ScopedLock<Mutex> lock(me.mutex);
//        item->ready = true;
//        if (item != me.first)
//            return true;
//    }
//
//    return true;

//    std::cerr<<"RELEASE\t"<<(size_t)ptr<<std::endl;


    // changing me.first or the ready flag must be done synchronized (me.mutex)
    // the thread who changed me.first->ready to be true has to write it.

    // change our ready flag and test if me.first->ready became true
    {
        ScopedLock<Mutex> lock(me.mutex);
        item->ready = true;
        if (item != me.first)
            return true;
    }

    // ok, if we are here it seems that we are responsible for writing the buffer

    SEQAN_ASSERT(me.first != NULL);

    bool success;
    do
    {
        success = me.worker(item->val);
//        std::cerr<<"WRITE"<<std::endl;
//success=true;

        {
            ScopedLock<Mutex> lock(me.mutex);
            me.first = item->next;
            // recycle released items
            releaseValue(me.pool, item);
//    std::cerr<<"FREE\t"<<(size_t)ptr<<std::endl;

            item = me.first;

            // can we leave?
            if (item == NULL || !item->ready)
                return success;
        }
        // we continue to write the next buffer
    }
    while (success);

    return false;
}

/** \brief A stream decorator that takes raw input and zips it to a ostream.

The class wraps up the inflate method of the bgzf library 1.1.4 http://www.gzip.org/bgzf/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_streambuf : public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef typename Tr::char_type char_type;
	typedef typename Tr::int_type int_type;

    typedef ConcurrentQueue<size_t, Suspendable<Limit> > TJobQueue;

    struct OutputBuffer
    {
        char    buffer[BGZF_MAX_BLOCK_SIZE];
        size_t  size;
    };

    struct BufferWriter
    {
        ostream_reference ostream;

        BufferWriter(ostream_reference ostream) :
            ostream(ostream)
        {}

        bool operator() (OutputBuffer const & outputBuffer)
        {
//    std::cerr<<"work\t"<<(size_t)&outputBuffer<<std::endl;
            ostream.write(outputBuffer.buffer, outputBuffer.size);
            return ostream.good();
        }
    };

    struct CompressionJob
    {
        typedef std::vector<char_type, char_allocator_type> TBuffer;

        TBuffer         buffer;
        size_t          size;
        OutputBuffer    *outputBuffer;

        CompressionJob() :
            buffer(BGZF_BLOCK_SIZE / sizeof(char_type), 0)
        {}
    };

    // string of recycable jobs
    size_t                  numThreads;
    size_t                  numJobs;
    String<CompressionJob>  jobs;
    TJobQueue               jobQueue;
    TJobQueue               idleQueue;
    Serializer<
        OutputBuffer,
        BufferWriter>       serializer;

    size_t                  currentJobId;
    bool                    currentJobAvail;


    struct CompressionThread
    {
        basic_bgzf_streambuf            *streamBuf;
        CompressionContext<BgzfFile>    compressionCtx;

        void operator()()
        {
            ScopedReadLock<TJobQueue> readLock(streamBuf->jobQueue);
            ScopedWriteLock<TJobQueue> writeLock(streamBuf->idleQueue);

            // wait for a new job to become available
            bool success = true;
            while (success)
            {
                size_t jobId;
                if (!popFront(jobId, streamBuf->jobQueue))
                    return;

                CompressionJob &job = streamBuf->jobs[jobId];
                OutputBuffer *outputBuffer = aquireValue(streamBuf->serializer);

                // compress block with zlib
                outputBuffer->size = _compressBlock(
                    outputBuffer->buffer, sizeof(outputBuffer->buffer),
                    &job.buffer[0], job.size, compressionCtx);

                success = releaseValue(streamBuf->serializer, outputBuffer);
                appendValue(streamBuf->idleQueue, jobId);
            }
            std::cerr<<"DAMN"<<std::endl;
        }
    };

    // array of worker threads
    Thread<CompressionThread>   *threads;

    /** Construct a zip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
    basic_bgzf_streambuf(ostream_reference ostream_,
                         size_t numThreads = 16,
                         size_t jobsPerThread = 8) :
        numThreads(numThreads),
        numJobs(numThreads * jobsPerThread),
        jobQueue(numJobs),
        idleQueue(numJobs),
		serializer(ostream_)
    {
        resize(jobs, numJobs, Exact());
        currentJobId = 0;

        lockWriting(jobQueue);
        lockReading(idleQueue);
        setReaderWriterCount(jobQueue, numThreads, 1);
        setReaderWriterCount(idleQueue, 1, numThreads);

        for (unsigned i = 0; i < numJobs; ++i)
        {
            bool success = appendValue(idleQueue, i);
            ignoreUnusedVariableWarning(success);
            SEQAN_ASSERT(success);
        }

        threads = new Thread<CompressionThread>[numThreads];
        for (unsigned i = 0; i < numThreads; ++i)
        {
            threads[i].worker.streamBuf = this;
            run(threads[i]);
        }

        currentJobAvail = popFront(currentJobId, idleQueue);
        SEQAN_ASSERT(currentJobAvail);

        CompressionJob &job = jobs[currentJobId];
		this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
    }

    ~basic_bgzf_streambuf()
    {
        // the buffer is now (after addFooter()) and flush will append the empty EOF marker
        flush(true);

        unlockWriting(jobQueue);
        unlockReading(idleQueue);

        for (unsigned i = 0; i < numThreads; ++i)
            waitFor(threads[i]);
        delete[] threads;
    }

    bool compressBuffer(size_t size)
    {
        // submit current job
        if (currentJobAvail)
        {
            jobs[currentJobId].size = size;
            appendValue(jobQueue, currentJobId);
        }

        // recycle existing idle job
        if (!(currentJobAvail = popFront(currentJobId, idleQueue)))
            return false;

        return serializer;
    }

    int_type overflow(int_type c)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if (c != EOF)
        {
            *this->pptr() = c;
            ++w;
        }
        if (compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
            return c;
        }
        else
        {
            return EOF;
        }
    }

	/** flushes the zip buffer and output buffer.

	This method should be called at the end of the compression. Calling flush multiple times, will lower the
	compression ratio.
	*/
	std::streamsize flush(bool flushEmptyBuffer = false)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if ((w != 0 || flushEmptyBuffer) && compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
        }
        else
        {
            w = 0;
        }

        // wait for running compressor threads
        waitForMinSize(idleQueue, numJobs - 1);

		serializer.worker.ostream.flush();
		return w;
    }

	int sync()
    {
		if (this->pptr() != this->pbase())
		{
			int c = overflow(EOF);
			if (c == EOF)
				return -1;
        }
        return 0;
    }

    void addFooter()
    {
        // we flush the filled buffer here, so that an empty (EOF) buffer is flushed in the d'tor
        if (this->pptr() != this->pbase())
            overflow(EOF);
    }

	/// returns a reference to the output stream
	ostream_reference get_ostream() const	{ return serializer.worker.ostream; };
};

/** \brief A stream decorator that takes compressed input and unzips it to a istream.

The class wraps up the deflate method of the bgzf library 1.1.4 http://www.gzip.org/bgzf/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_unbgzf_streambuf : 
	public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef typename Tr::char_type char_type;
	typedef typename Tr::int_type int_type;
	typedef typename Tr::off_type off_type;
	typedef typename Tr::pos_type pos_type;

    typedef std::vector<char_type, char_allocator_type> TBuffer;
    typedef ConcurrentQueue<size_t, Suspendable<Limit> >  TJobQueue;

    static const size_t MAX_PUTBACK = 4;

    struct Serializer
    {
        istream_reference   istream;
        Mutex               lock;
        std::exception      *error;
        off_t               fileOfs;

        Serializer(istream_reference istream) :
            istream(istream),
            lock(false),
            error(NULL),
            fileOfs(0u)
        {}

        ~Serializer()
        {
            delete error;
        }
    };

    Serializer serializer;

    struct DecompressionJob
    {
        typedef std::vector<byte_type, byte_allocator_type> TInputBuffer;

        TInputBuffer    inputBuffer;
        TBuffer         buffer;
        off_t           fileOfs;
        size_t          size;

        Mutex           mutex;
        Event2          readyEvent;
        bool            ready;

        DecompressionJob() :
            inputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(MAX_PUTBACK + BGZF_MAX_BLOCK_SIZE / sizeof(char_type), 0),
            size(0),
            mutex(false),
            readyEvent(mutex),
            ready(false)
        {}

        DecompressionJob(DecompressionJob const &other) :
            inputBuffer(other.inputBuffer),
            buffer(other.buffer),
            size(other.size),
            mutex(false),
            readyEvent(mutex),
            ready(other.ready)
        {}
    };

    // string of recycable jobs
    size_t                      numThreads;
    size_t                      numJobs;
    String<DecompressionJob>    jobs;
    TJobQueue                   runningQueue;
    TJobQueue                   idleQueue;
    size_t                      currentJobId;
    bool                        currentJobAvail;

    struct DecompressionThread
    {
        basic_unbgzf_streambuf          *streamBuf;
        CompressionContext<BgzfFile>    compressionCtx;

        void operator()()
        {
            ScopedReadLock<TJobQueue> readLock(streamBuf->idleQueue);
            ScopedWriteLock<TJobQueue> writeLock(streamBuf->runningQueue);

            // wait for a new job to become available
            while (true)
            {
                size_t jobId;
                if (!popFront(jobId, streamBuf->idleQueue))
                    return;

                DecompressionJob &job = streamBuf->jobs[jobId];
                size_t tailLen;

                {
                    ScopedLock<Mutex> scopedLock(streamBuf->serializer.lock);

                    if (streamBuf->serializer.error != NULL)
                        return;

                    // remember start offset (for tellg later)
                    job.fileOfs = streamBuf->serializer.fileOfs;

                    // read header
                    streamBuf->serializer.istream.read(
                        (char*)&(job.inputBuffer[0]),
                        BGZF_BLOCK_HEADER_LENGTH);

                    if (!streamBuf->serializer.istream.good())
                    {
                        if (!streamBuf->serializer.istream.eof())
                            streamBuf->serializer.error = new IOException("Stream read error.");
                        return;
                    }

                    // check header
                    if (!_bgzfCheckHeader(&(job.inputBuffer[0])))
                    {
                        streamBuf->serializer.error = new IOException("Invalid BGZF block header.");
                        return;
                    }

                    // extract length of compressed data
                    tailLen = _bgzfUnpack16(&(job.inputBuffer[16])) + 1u - BGZF_BLOCK_HEADER_LENGTH;

                    // read compressed data and tail
                    streamBuf->serializer.istream.read(
                        (char*)&(job.inputBuffer[BGZF_BLOCK_HEADER_LENGTH]),
                        tailLen);
                    
                    if (!streamBuf->serializer.istream.good())
                    {
                        if (!streamBuf->serializer.istream.eof())
                            streamBuf->serializer.error = new IOException("Stream read error.");
                        return;
                    }

                    streamBuf->serializer.fileOfs += BGZF_BLOCK_HEADER_LENGTH + tailLen;
                    job.ready = false;

                    if (!appendValue(streamBuf->runningQueue, jobId))
                        return;
                }

                // decompress block
                job.size = _decompressBlock(
                    &job.buffer[MAX_PUTBACK], capacity(job.buffer),
                    &job.inputBuffer[0], BGZF_BLOCK_HEADER_LENGTH + tailLen, compressionCtx);

                // signal that job is ready
                {
                    ScopedLock<Mutex> lock(job.mutex);
                    job.ready = true;
                    signal(job.readyEvent);
                }
            }
        }
    };

    // array of worker threads
    Thread<DecompressionThread> *threads;
    TBuffer                     putbackBuffer;

    /** Construct a unzip stream
    * More info on the following parameters can be found in the bgzf documentation.
    */
    basic_unbgzf_streambuf(istream_reference istream_,
                           size_t numThreads = 16,
                           size_t jobsPerThread = 8) :
		serializer(istream_),
        numThreads(numThreads),
        numJobs(numThreads * jobsPerThread),
        runningQueue(numJobs),
        idleQueue(numJobs),
        putbackBuffer(MAX_PUTBACK)
    {
        resize(jobs, numJobs, Exact());
        currentJobId = 0;
        currentJobAvail = false;

        lockReading(runningQueue);
        lockWriting(idleQueue);
        setReaderWriterCount(runningQueue, 1, numThreads);
        setReaderWriterCount(idleQueue, numThreads, 1);

        for (unsigned i = 0; i < numJobs; ++i)
        {
            bool success = appendValue(idleQueue, i);
            ignoreUnusedVariableWarning(success);
            SEQAN_ASSERT(success);
        }

        threads = new Thread<DecompressionThread>[numThreads];
        for (unsigned i = 0; i < numThreads; ++i)
        {
            threads[i].worker.streamBuf = this;
            run(threads[i]);
        }
    }

	~basic_unbgzf_streambuf()
    {
        unlockWriting(idleQueue);
        unlockReading(runningQueue);

        for (unsigned i = 0; i < numThreads; ++i)
            waitFor(threads[i]);
        delete[] threads;
        std::cerr << "runnning queue:\t" << (__int64)&runningQueue << std::endl;
        std::cerr << "idle queue:    \t" << (__int64)&idleQueue << std::endl;
    }

    int_type underflow()
    {
        // no need to use the next buffer?
        if (this->gptr() && this->gptr() < this->egptr())
            return Tr::to_int_type(*this->gptr());

        size_t putback = this->gptr() - this->eback();
        if (putback > MAX_PUTBACK)
            putback = MAX_PUTBACK;

        // save at most MAX_PUTBACK characters from previous page to putback buffer
        if (putback != 0)
            std::copy(
                this->gptr() - putback,
                this->gptr(),
                &putbackBuffer[0]);

        if (currentJobAvail)
            appendValue(idleQueue, currentJobId);

        while (true)
        {
            if (!(currentJobAvail = popFront(currentJobId, runningQueue)))
            {
                if (serializer.error != NULL)
                    throw *serializer.error;
                return EOF;
            }

            DecompressionJob &job = jobs[currentJobId];

            // restore putback buffer
            this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
            if (putback != 0)
                std::copy(
                    &putbackBuffer[0],
                    &putbackBuffer[putback],
                    &job.buffer[MAX_PUTBACK - putback]);

            // wait for the end of decompression
            {
                ScopedLock<Mutex> lock(job.mutex);
                if (!job.ready)
                    waitFor(job.readyEvent);
            }

            // reset buffer pointers
            this->setg( 
                  &job.buffer[MAX_PUTBACK - putback],      // beginning of putback area
                  &job.buffer[MAX_PUTBACK],                // read position
                  &job.buffer[MAX_PUTBACK + job.size]);    // end of buffer

            if (job.size != 0)
                break;
        }

        // return next character
        return Tr::to_int_type(*this->gptr());
    }

    pos_type seekoff(off_type ofs, std::ios_base::seekdir dir, std::ios_base::openmode openMode)
    {
        DecompressionJob &job = jobs[currentJobId];
        if ((openMode & (std::ios_base::in | std::ios_base::out)) == std::ios_base::in)
        {
            if (dir == std::ios_base::cur && ofs == 0)
                return pos_type((job.fileOfs << 16) + (this->gptr() - this->eback()));
        }
        return pos_type(off_type(-1));
    }

    pos_type seekpos(pos_type pos, std::ios_base::openmode openMode)
    {
        return seekoff(off_type(pos), std::ios_base::beg, openMode);
    }

	/// returns the compressed input istream
	istream_reference get_istream()	{ return serializer.istream;};
};

/*! \brief Base class for zip ostreams

Contains a basic_bgzf_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
	typedef basic_bgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > bgzf_streambuf_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
	basic_bgzf_ostreambase(ostream_reference ostream_)
		: m_buf(ostream_)
	{
		this->init(&m_buf );
	};
	
	/// returns the underlying zip ostream object
	bgzf_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the bgzf error state
	int get_zerr() const					{	return m_buf.get_err();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the compressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the uncompressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	bgzf_streambuf_type m_buf;
};

/*! \brief Base class for unzip istreams

Contains a basic_unbgzf_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
	typedef basic_unbgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > unbgzf_streambuf_type;

	basic_bgzf_istreambase(istream_reference ostream_)
		: m_buf(ostream_)
	{
		this->init(&m_buf );
	};
	
	/// returns the underlying unzip istream object
	unbgzf_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the bgzf error state
	int get_zerr() const					{	return m_buf.get_zerr();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the uncompressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the compressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	unbgzf_streambuf_type m_buf;
};

/*! \brief A zipper ostream

This class is a ostream decorator that behaves 'almost' like any other ostream.

At construction, it takes any ostream that shall be used to output of the compressed data.

When finished, you need to call the special method zflush or call the destructor 
to flush all the intermidiate streams.

Example:
\code
// creating the target zip string, could be a fstream
ostringstream ostringstream_;
// creating the zip layer
bgzf_ostream zipper(ostringstream_);

	
// writing data	
zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// zip ostream needs special flushing...
zipper.zflush();
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_ostream : 
	public basic_bgzf_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_ostream<Elem,Tr>
{
public:
	typedef basic_bgzf_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_ostreambase_type;
	typedef std::basic_ostream<Elem,Tr> ostream_type;
    typedef ostream_type& ostream_reference;

	/** Constructs a zipper ostream decorator
	 *
	 * \param ostream_ ostream where the compressed output is written

	 When is_gbgzf_ is true, a gzip header and footer is automatically added.
	 */
	basic_bgzf_ostream(ostream_reference ostream_)
	: 
		bgzf_ostreambase_type(ostream_),
		ostream_type(this->rdbuf())
	{}

	/// flush inner buffer and zipper buffer
	basic_bgzf_ostream<Elem,Tr>& zflush()	
	{	
		this->flush(); this->rdbuf()->flush(); return *this;
	};

    ~basic_bgzf_ostream()
    {
        this->rdbuf()->addFooter();
    }


private:
    static void put_long(ostream_reference out_, unsigned long x_);
};

/*! \brief A zipper istream

This class is a istream decorator that behaves 'almost' like any other ostream.

At construction, it takes any istream that shall be used to input of the compressed data.

Simlpe example:
\code
// create a stream on zip string
istringstream istringstream_( ostringstream_.str());
// create unzipper istream
bgzf_istream unzipper( istringstream_);

// read and unzip
unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istream : 
	public basic_bgzf_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_istream<Elem,Tr>
{
public:
	typedef basic_bgzf_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_istreambase_type;
	typedef std::basic_istream<Elem,Tr> istream_type;
    typedef istream_type& istream_reference;
	typedef char byte_type;

	/** Construct a unzipper stream
	 *
	 * \param istream_ input buffer
	 */
	basic_bgzf_istream(istream_reference istream_)
	  : 
		bgzf_istreambase_type(istream_),
		istream_type(this->rdbuf()),
		m_is_gzip(false),
		m_gbgzf_data_size(0)
	{};

	/// returns true if it is a gzip file
	bool is_gzip() const				{	return m_is_gzip;};
	/// return data size check
	bool check_data_size() const		{	return this->get_out_size() == m_gbgzf_data_size;};

	/// return the data size in the file 
	long get_gbgzf_data_size() const		{	return m_gbgzf_data_size;};
protected:
    static void read_long(istream_reference in_, unsigned long& x_);

	int check_header();
	bool m_is_gzip;
	unsigned long m_gbgzf_data_size;
};

/// A typedef for basic_bgzf_ostream<char>
typedef basic_bgzf_ostream<char> bgzf_ostream;
/// A typedef for basic_bgzf_ostream<wchar_t>
typedef basic_bgzf_ostream<wchar_t> bgzf_wostream;
/// A typedef for basic_bgzf_istream<char>
typedef basic_bgzf_istream<char> bgzf_istream;
/// A typedef for basic_bgzf_istream<wchart>
typedef basic_bgzf_istream<wchar_t> bgzf_wistream;

}  // namespace seqan

#include "bgzfstream_impl.h"

#endif
