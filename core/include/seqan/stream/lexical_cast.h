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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// String <=> Numerical conversions
// ==========================================================================

#ifndef SEQAN_STREAM_LEXICAL_CAST_H
#define SEQAN_STREAM_LEXICAL_CAST_H

namespace seqan {

// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Exception BadLexicalCast
// ----------------------------------------------------------------------------

struct BadLexicalCast : ParseError
{
    template <typename TTarget, typename TSource>
    BadLexicalCast(TTarget const & target, TSource const & source) :
        ParseError(std::string("Unable to convert '") + toCString(source) + "' into " + toCString(Demangler<TTarget const>(target)) + ".")
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IntegerFormatString_
// ----------------------------------------------------------------------------
// Return the format string for numbers.

template <typename TUnsigned, unsigned SIZE, typename T = void>
struct IntegerFormatString_;


template <typename T>
struct IntegerFormatString_<False, 1, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 1, T>::VALUE[] = "%hhi%n";


template <typename T>
struct IntegerFormatString_<True, 1, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 1, T>::VALUE[] = "%hhu%n";


template <typename T>
struct IntegerFormatString_<False, 2, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 2, T>::VALUE[] = "%hi%n";


template <typename T>
struct IntegerFormatString_<True, 2, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 2, T>::VALUE[] = "%hu%n";


template <typename T>
struct IntegerFormatString_<False, 4, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 4, T>::VALUE[] = "%i%n";


template <typename T>
struct IntegerFormatString_<True, 4, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 4, T>::VALUE[] = "%u%n";


template <typename T>
struct IntegerFormatString_<False, 8, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 8, T>::VALUE[] = "%lli%n";


template <typename T>
struct IntegerFormatString_<True, 8, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 8, T>::VALUE[] = "%llu%n";

// ============================================================================
// Functions
// ============================================================================

// TODO(dox): lexicalCast2
/*!
 * @fn lexicalCast2
 * @headerfile seqan/stream.h
 * @brief Cast from a String-type to a numerical type
 * 
 * @signature bool lexicalCast2(target, source);
 * 
 * @param[out] target Object to hold result of cast.
 * @param[in]  source The string to be read from.  Type: @link SequenceConcept @endlink.

 * @return bool <tt>true</tt> if cast was successful, <tt>false</tt> otherwise.
 * 
 * @section Remarks
 * 
 * Uses istringstream internally, so right now "123foobar" will be succesfully cast to an int of 123.
 * 
 * @section Examples
 * 
 * Using lexicalCast2 is straightforward and we can detect errors.
 * 
 * @code{.cpp}
 * unsigned u = 0;
 * int i = 0;
 * double = 0;
 * bool success = false;
 *  
 * success = lexicalCast2(u, "3");      // => success is true, u is 3.
 * success = lexicalCast2(u, "-3");     // => success is false, u is undefined.
 * success = lexicalCast2(i, "-3");     // => success is true, i is -3.
 * success = lexicalCast2(d, "-3.99");  // => success is true, d is -3.99.
 * @endcode
 *
 * @see lexicalCast
 */

/**
.Function.lexicalCast2
..cat:Input/Output
..summary:Cast from a String-type to a numerical type
..signature:lexicalCast2(TTarget & target, TSource const & source)
..signature:lexicalCast2(TTarget & target, String<TValue, TSpec> const & source)
..param.target:Object to hold result of cast
...type:nolink:$int$
...type:nolink:$unsigned int$
...type:nolink:$double$
...type:nolink:or similar
..param.source:The string to be read from
...type:Shortcut.CharString
...type:nolink:char[]
...type:nolink:std::string
...type:nolink:or similar
..returns:$true$ if cast was successful, $false$ otherwise
...type:nolink:$bool$
..remarks:uses istringstream internally, so right now "123foobar" will be
succesfully cast to an int of 123
..include:seqan/stream.h
..see:Function.lexicalCast
..example.text:Using lexicalCast2 is straightforward and we can detect errors.
..example.code:
unsigned u = 0;
int i = 0;
double = 0;
bool success = false;

success = lexicalCast2(u, "3");      // => success is true, u is 3.
success = lexicalCast2(u, "-3");     // => success is false, u is undefined.
success = lexicalCast2(i, "-3");     // => success is true, i is -3.
success = lexicalCast2(d, "-3.99");  // => success is true, d is -3.99.
 */

// ----------------------------------------------------------------------------
// Function lexicalCast()
// ----------------------------------------------------------------------------
// Generic version for integers.

template <typename TInteger, typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<SignedIntegerConcept<TInteger> >, bool)
lexicalCast(TInteger & target, TSource const & source)
{
    int offset;
    return (sscanf(toCString(source), IntegerFormatString_<False, sizeof(TInteger)>::VALUE, &target, &offset) == 1) &&
           (static_cast<typename Size<TSource>::Type>(offset) == length(source));
}

template <typename TInteger, typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<UnsignedIntegerConcept<TInteger> >, bool)
lexicalCast(TInteger & target, TSource const & source)
{
    if (!empty(source) && front(source) == '-') return false;

    int offset;
    return (sscanf(toCString(source), IntegerFormatString_<True, sizeof(TInteger)>::VALUE, &target, &offset) == 1) &&
           (static_cast<typename Size<TSource>::Type>(offset) == length(source));
}

// ----------------------------------------------------------------------------
// Function lexicalCast(float)
// ----------------------------------------------------------------------------

template <typename TSource>
inline bool lexicalCast(float & target, TSource const & source)
{
    int offset;
    return (sscanf(toCString(source), "%g%n", &target, &offset) == 1) &&
           (static_cast<typename Size<TSource>::Type>(offset) == length(source));
}

// ----------------------------------------------------------------------------
// Function lexicalCast(double)
// ----------------------------------------------------------------------------

template <typename TSource>
inline bool lexicalCast(double & target, TSource const & source)
{
    int offset;
    return (sscanf(toCString(source), "%lg%n", &target, &offset) == 1) &&
           (static_cast<typename Size<TSource>::Type>(offset) == length(source));
}

// TODO(dox): lexicalCast
/*!
 * @fn lexicalCast
 * @headerfile <seqan/stream.h>
 * @brief Cast from a String-type to a numerical type
 * 
 * @signature template <typename TTarget>
 *            TTarget lexicalCast<TTarget>(source);
 * 
 * @tparam TTarget Target type to cast to.
 *
 * @param[in] source The string to be read from.  Type: @link SequenceConcept @endlink.
 * 
 * @return TTarget Value of Type TTarget with cast contents of source.
 * 
 * @section Remarks
 * 
 * Return value undefined if casting fails, see @link lexicalCast2 @endlink for a more robust variant.
 * 
 * This function uses <tt>std::istringstream</tt> internally, so right now "123foobar" will be succesfully cast to an
 * int of 123.
 * 
 * @section Examples
 * 
 * Using <tt>lexicalCast<>()</tt> is easy but not as robust as @link lexicalCast2 @endlink: We cannot detect parsing or
 * conversion errors.
 * 
 * @code{.cpp}
 * unsigned u = 0;
 * int i = 0;
 * double = 0;
 * bool success = false;
 *  
 * u = lexicalCast<unsigned>( "3");   // => u is 3.
 * u = lexicalCast<unsigned>("-3");   // => u is undefined.
 * i = lexicalCast<int>("-3");        // => i is -3.
 * d = lexicalCast<double>("-3.99");  // => d is -3.99.
 * @endcode
 *
 * @see lexicalCast2
 */

/**
.Function.lexicalCast
..cat:Input/Output
..summary:Cast from a String-type to a numerical type
..signature:lexicalCast<TTarget>(TSource const & source)
..signature:lexicalCast<TTarget>(String<TValue, TSpec> const & source)
..param.source:The string to be read from
...type:Shortcut.CharString
...type:nolink:char[]
...type:nolink:std::string
...type:nolink:or similar
..param.TTarget:Type to be casted to
...type:nolink:$int$
...type:nolink:$unsigned int$
...type:nolink:$double$
...type:nolink:or similar
..returns:Value of Type TTarget with casted contents of source
...type:nolink:TTarget
..remarks:Return value undefined if casting fails, see @Function.lexicalCast2@ for a more robust variant.
..remarks:This function uses $std::istringstream$ internally, so right now "123foobar" will be
succesfully cast to an int of 123.
..include:seqan/stream.h
..see:Function.lexicalCast2
..example.text:Using $lexicalCast<>()$ is easy but not as robust as @Function.lexicalCast2@: We cannot detect parsing or conversion errors.
..example.code:
unsigned u = 0;
int i = 0;
double = 0;
bool success = false;

u = lexicalCast<unsigned>( "3");   // => u is 3.
u = lexicalCast<unsigned>("-3");   // => u is undefined.
i = lexicalCast<int>("-3");        // => i is -3.
d = lexicalCast<double>("-3.99");  // => d is -3.99.
 */

template <typename TTarget, typename TSource>
inline TTarget lexicalCast(TSource const & source)
{
    TTarget target;
    if (!lexicalCast(target, source))
        throw BadLexicalCast(target, source);
    return target;
}

// ----------------------------------------------------------------------------
// Function appendNumber()
// ----------------------------------------------------------------------------
// Generic version for integers.

template <typename TTarget, typename TInteger>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TInteger> >, typename Size<TTarget>::Type)
appendNumber(TTarget & target, TInteger i)
{
    // 1 byte has at most 3 decimal digits (plus 1 for the NULL character)
    char buffer[sizeof(TInteger) * 3 + 2];
    int offset;
    size_t len = snprintf(buffer, sizeof(buffer),
                          IntegerFormatString_<typename Is<UnsignedIntegerConcept<TInteger> >::Type,
                          sizeof(TInteger)>::VALUE, i, &offset);

    Range<char *> range = toRange(buffer + 0, buffer + len);
    write(target, range);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(float)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, float source)
{
    char buffer[32];
    int offset;
    size_t len = snprintf(buffer, 32, "%g%n", source, &offset);
    Range<char *> range = toRange(buffer + 0, buffer + len);
    write(target, range);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(double)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, double source)
{
    char buffer[32];
    int offset;
    size_t len = snprintf(buffer, 32, "%lg%n", source, &offset);
    Range<char *> range = toRange(buffer + 0, buffer + len);
    write(target, range);
    return len;
}

}

#endif //def SEQAN_STREAM_LEXICAL_CAST_H
