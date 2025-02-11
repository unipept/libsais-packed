/*--

This file is a part of libsais, a library for linear time suffix array,
longest common prefix array and burrows wheeler transform construction.

   Copyright (c) 2021-2024 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information.

--*/

/*--
Modifications made by team Unipept as of 18/10/2024 and ongoing.
--*/

#include "libsais16x64.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>

typedef int64_t                         sa_sint_t;
typedef uint64_t                        sa_uint_t;
typedef int64_t                         fast_sint_t;
typedef uint64_t                        fast_uint_t;

#define SAINT_BIT                       (64)
#define SAINT_MAX                       INT64_MAX
#define SAINT_MIN                       INT64_MIN

#define ALPHABET_SIZE                   (1 << CHAR_BIT << CHAR_BIT)
#define UNBWT_FASTBITS                  (17)

#define SUFFIX_GROUP_BIT                (SAINT_BIT - 1)
#define SUFFIX_GROUP_MARKER             (((sa_sint_t)1) << (SUFFIX_GROUP_BIT - 1))

#define BUCKETS_INDEX2(_c, _s)          ((((fast_sint_t)_c) << 1) + (fast_sint_t)(_s))
#define BUCKETS_INDEX4(_c, _s)          ((((fast_sint_t)_c) << 2) + (fast_sint_t)(_s))

#define LIBSAIS_LOCAL_BUFFER_SIZE       (1024)
#define LIBSAIS_PER_THREAD_CACHE_SIZE   (2097184)

typedef struct LIBSAIS_THREAD_CACHE
{
        sa_sint_t                       symbol;
        sa_sint_t                       index;
} LIBSAIS_THREAD_CACHE;

typedef union LIBSAIS_THREAD_STATE
{
    struct
    {
        fast_sint_t                     position;
        fast_sint_t                     count;

        fast_sint_t                     m;
        fast_sint_t                     last_lms_suffix;

        sa_sint_t *                     buckets;
        LIBSAIS_THREAD_CACHE *          cache;
    } state;

    uint8_t padding[64];
} LIBSAIS_THREAD_STATE;

typedef struct LIBSAIS_CONTEXT
{
    sa_sint_t *                         buckets;
} LIBSAIS_CONTEXT;

typedef struct LIBSAIS_UNBWT_CONTEXT
{
    sa_uint_t *                         bucket2;
    uint16_t *                          fastbits;
    sa_uint_t *                         buckets;
} LIBSAIS_UNBWT_CONTEXT;

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
    #define RESTRICT __restrict
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if defined(__has_builtin)
    #if __has_builtin(__builtin_prefetch)
        #define HAS_BUILTIN_PREFETCH
    #endif
#elif defined(__GNUC__) && (((__GNUC__ == 3) && (__GNUC_MINOR__ >= 2)) || (__GNUC__ >= 4))
    #define HAS_BUILTIN_PREFETCH
#endif 

#if defined(HAS_BUILTIN_PREFETCH)
    #define libsais16x64_prefetchr(address) __builtin_prefetch((const void *)(address), 0, 3)
    #define libsais16x64_prefetchw(address) __builtin_prefetch((const void *)(address), 1, 3)
#elif defined (_M_IX86) || defined (_M_AMD64)
    #include <intrin.h>
    #define libsais16x64_prefetchr(address) _mm_prefetch((const void *)(address), _MM_HINT_T0)
    #define libsais16x64_prefetchw(address) _m_prefetchw((const void *)(address))
#elif defined (_M_ARM)
    #include <intrin.h>
    #define libsais16x64_prefetchr(address) __prefetch((const void *)(address))
    #define libsais16x64_prefetchw(address) __prefetchw((const void *)(address))
#elif defined (_M_ARM64)
    #include <intrin.h>
    #define libsais16x64_prefetchr(address) __prefetch2((const void *)(address), 0)
    #define libsais16x64_prefetchw(address) __prefetch2((const void *)(address), 16)
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if !defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
    #if defined(_LITTLE_ENDIAN) \
            || (defined(BYTE_ORDER) && defined(LITTLE_ENDIAN) && BYTE_ORDER == LITTLE_ENDIAN) \
            || (defined(_BYTE_ORDER) && defined(_LITTLE_ENDIAN) && _BYTE_ORDER == _LITTLE_ENDIAN) \
            || (defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && __BYTE_ORDER == __LITTLE_ENDIAN) \
            || (defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
        #define __LITTLE_ENDIAN__
    #elif defined(_BIG_ENDIAN) \
            || (defined(BYTE_ORDER) && defined(BIG_ENDIAN) && BYTE_ORDER == BIG_ENDIAN) \
            || (defined(_BYTE_ORDER) && defined(_BIG_ENDIAN) && _BYTE_ORDER == _BIG_ENDIAN) \
            || (defined(__BYTE_ORDER) && defined(__BIG_ENDIAN) && __BYTE_ORDER == __BIG_ENDIAN) \
            || (defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
        #define __BIG_ENDIAN__
    #elif defined(_WIN32)
        #define __LITTLE_ENDIAN__
    #endif
#endif

static void * libsais16x64_align_up(const void * address, size_t alignment)
{
    return (void *)((((ptrdiff_t)address) + ((ptrdiff_t)alignment) - 1) & (-((ptrdiff_t)alignment)));
}

static void * libsais16x64_alloc_aligned(size_t size, size_t alignment)
{
    void * address = malloc(size + sizeof(short) + alignment - 1);
    if (address != NULL)
    {
        void * aligned_address = libsais16x64_align_up((void *)((ptrdiff_t)address + (ptrdiff_t)(sizeof(short))), alignment);
        ((short *)aligned_address)[-1] = (short)((ptrdiff_t)aligned_address - (ptrdiff_t)address);

        return aligned_address;
    }

    return NULL;
}

static void libsais16x64_free_aligned(void * aligned_address)
{
    if (aligned_address != NULL)
    {
        free((void *)((ptrdiff_t)aligned_address - ((short *)aligned_address)[-1]));
    }
}

static void libsais16x64_gather_lms_suffixes_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, fast_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 128;

        fast_sint_t i, j = omp_block_start + omp_block_size, c0 = T[omp_block_start + omp_block_size - 1], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t s = c0 >= c1;

        for (i = omp_block_start + omp_block_size - 2, j = omp_block_start + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchr(&T[i - prefetch_distance]);

            c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 0); m -= ((s & 3) == 1);
            c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 1); m -= ((s & 3) == 1);
            c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 2); m -= ((s & 3) == 1);
        }

        for (j -= 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
        }

        SA[m] = (sa_sint_t)(i + 1);
    }
}

static void libsais16x64_gather_lms_suffixes_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = n - omp_block_start;

    libsais16x64_gather_lms_suffixes_16u(T, SA, n, (fast_sint_t)n - 1, omp_block_start, omp_block_size);
}

static sa_sint_t libsais16x64_gather_lms_suffixes_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t             i   = n - 2;
    sa_sint_t             m   = n - 1;
    fast_uint_t           s   = 1;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

    for (; i >= 3; i -= 4)
    {
        libsais16x64_prefetchr(&T[i - prefetch_distance]);

        c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = i - 0; m -= ((s & 3) == 1);
        c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = i - 1; m -= ((s & 3) == 1);
        c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = i - 2; m -= ((s & 3) == 1);
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
    }

    return n - 1 - m;
}

static sa_sint_t libsais16x64_gather_compacted_lms_suffixes_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t             i   = n - 2;
    sa_sint_t             m   = n - 1;
    fast_uint_t           s   = 1;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

    for (; i >= 3; i -= 4)
    {
        libsais16x64_prefetchr(&T[i - prefetch_distance]);

        c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = i + 1; m -= ((fast_sint_t)(s & 3) == (c0 >= 0));
        c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = i - 0; m -= ((fast_sint_t)(s & 3) == (c1 >= 0));
        c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = i - 1; m -= ((fast_sint_t)(s & 3) == (c0 >= 0));
        c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = i - 2; m -= ((fast_sint_t)(s & 3) == (c1 >= 0));
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = i + 1; m -= ((fast_sint_t)(s & 3) == (c1 >= 0));
    }

    return n - 1 - m;
}

static void libsais16x64_count_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 32;

    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    sa_sint_t             i   = n - 2;
    fast_uint_t           s   = 1;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais16x64_prefetchr(&T[i - 2 * prefetch_distance]);

        libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0], 0)]);
        libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1], 0)]);
        libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2], 0)]);
        libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3], 0)]);

        c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1)));
        buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1)));
        buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;

        c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1)));
        buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1)));
        buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1)));
        buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;
    }

    buckets[BUCKETS_INDEX2((fast_uint_t)c0, 0)]++;
}

static sa_sint_t libsais16x64_count_and_gather_lms_suffixes_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, (size_t)4 * ALPHABET_SIZE * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 128;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t s = c0 >= c1;

        for (i = m - 1, j = omp_block_start + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchr(&T[i - prefetch_distance]);

            c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, s & 3)]++;

            c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 0); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]++;

            c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, s & 3)]++;

            c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 2); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]++;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]++;
        }

        c1 = (i >= 0) ? T[i] : -1; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((fast_uint_t)c0, s & 3)]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais16x64_count_and_gather_lms_suffixes_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t m = 0;

    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = n - omp_block_start;

    m = libsais16x64_count_and_gather_lms_suffixes_16u(T, SA, n, buckets, omp_block_start, omp_block_size);

    return m;
}

static sa_sint_t libsais16x64_count_and_gather_lms_suffixes_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 4 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 32;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t s = c0 >= c1;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 0], 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 1], 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 2], 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 3], 0)]);

            c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, s & 3)]++;

            c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 0); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]++;

            c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, s & 3)]++;

            c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 2); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]++;
        }

        c1 = (i >= 0) ? T[i] : -1; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((fast_uint_t)c0, s & 3)]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais16x64_count_and_gather_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 32;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t s = c0 >= c1;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0], 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1], 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2], 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3], 0)]);

            c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;

            c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 0); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;

            c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;

            c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 2); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;
        }

        c1 = (i >= 0) ? T[i] : -1; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 32;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t s = c0 >= c1;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0] & SAINT_MAX, 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1] & SAINT_MAX, 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2] & SAINT_MAX, 0)]);
            libsais16x64_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3] & SAINT_MAX, 0)]);

            c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((fast_sint_t)(s & 3) == (c0 >= 0));
            c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;

            c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 0); m -= ((fast_sint_t)(s & 3) == (c1 >= 0));
            c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;

            c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 1); m -= ((fast_sint_t)(s & 3) == (c0 >= 0));
            c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;

            c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i - 2); m -= ((fast_sint_t)(s & 3) == (c1 >= 0));
            c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((fast_sint_t)(s & 3) == (c1 >= 0));
            c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (s & 3) == 1)]++;
        }

        c1 = (i >= 0) ? T[i] : -1; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); SA[m] = (sa_sint_t)(i + 1); m -= ((fast_sint_t)(s & 3) == (c0 >= 0));
        c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (s & 3) == 1)]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais16x64_count_and_gather_lms_suffixes_32s_4k_nofs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t m = libsais16x64_count_and_gather_lms_suffixes_32s_4k(T, SA, n, k, buckets, 0, n);

    return m;
}

static sa_sint_t libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_nofs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t m = 0;

    m = libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k(T, SA, n, k, buckets, 0, n);

    return m;
}

static sa_sint_t libsais16x64_count_and_gather_lms_suffixes_32s_4k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t m = libsais16x64_count_and_gather_lms_suffixes_32s_4k_nofs_omp(T, SA, n, k, buckets);

    return m;
}

static void libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_nofs_omp(T, SA, n, k, buckets);
}

static void libsais16x64_count_suffixes_32s(const sa_sint_t * RESTRICT T, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 32;

    memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

    fast_sint_t i, j;
    for (i = 0, j = (fast_sint_t)n - 7; i < j; i += 8)
    {
        libsais16x64_prefetchr(&T[i + prefetch_distance]);

        buckets[T[i + 0]]++;
        buckets[T[i + 1]]++;
        buckets[T[i + 2]]++;
        buckets[T[i + 3]]++;
        buckets[T[i + 4]]++;
        buckets[T[i + 5]]++;
        buckets[T[i + 6]]++;
        buckets[T[i + 7]]++;
    }

    for (j += 7; i < j; i += 1)
    {
        buckets[T[i]]++;
    }
}

static sa_sint_t libsais16x64_initialize_buckets_start_and_end_16u(sa_sint_t * RESTRICT buckets, sa_sint_t * RESTRICT freq)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[6 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT bucket_end   = &buckets[7 * ALPHABET_SIZE];

    fast_sint_t k = -1;

    if (freq != NULL)
    {
        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
        {
            sa_sint_t total = buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];

            bucket_start[j] = sum; sum += total; bucket_end[j] = sum; k = total > 0 ? j : k; freq[j] = total;
        }
    }
    else
    {
        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
        {
            sa_sint_t total = buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];

            bucket_start[j] = sum; sum += total; bucket_end[j] = sum; k = total > 0 ? j : k;
        }
    }

    return (sa_sint_t)(k + 1);
}

static void libsais16x64_initialize_buckets_start_and_end_32s_6k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[4 * (fast_sint_t)k];
    sa_sint_t * RESTRICT bucket_end   = &buckets[5 * (fast_sint_t)k];

    fast_sint_t i, j; sa_sint_t sum = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
    {
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];
        bucket_end[j] = sum;
    }
}

static void libsais16x64_initialize_buckets_start_and_end_32s_4k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[2 * (fast_sint_t)k];
    sa_sint_t * RESTRICT bucket_end   = &buckets[3 * (fast_sint_t)k];

    fast_sint_t i, j; sa_sint_t sum = 0;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    { 
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        bucket_end[j] = sum;
    }
}

static void libsais16x64_initialize_buckets_start_32s_1k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    fast_sint_t i; sa_sint_t sum = 0;
    for (i = 0; i <= (fast_sint_t)k - 1; i += 1) { sa_sint_t tmp = buckets[i]; buckets[i] = sum; sum += tmp; }
}

static void libsais16x64_initialize_buckets_end_32s_1k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    fast_sint_t i; sa_sint_t sum = 0;
    for (i = 0; i <= (fast_sint_t)k - 1; i += 1) { sum += buckets[i]; buckets[i] = sum; }
}

static sa_sint_t libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix)
{
    {
        fast_uint_t     s = 0;
        fast_sint_t     c0 = T[first_lms_suffix];
        fast_sint_t     c1 = 0;

        for (; --first_lms_suffix >= 0; )
        {
            c1 = c0; c0 = T[first_lms_suffix]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1)));
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]--;
        }

        buckets[BUCKETS_INDEX4((fast_uint_t)c0, (s << 1) & 3)]--;
    }

    {
        sa_sint_t * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
        { 
            temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum; sum += buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 3)]; temp_bucket[j] = sum;
        }

        return sum;
    }
}

static sa_sint_t libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix)
{
    {
        fast_uint_t     s = 0;
        fast_sint_t     c0 = T[first_lms_suffix];
        fast_sint_t     c1 = 0;

        for (; --first_lms_suffix >= 0; )
        {
            c1 = c0; c0 = T[first_lms_suffix]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1)));
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, s & 3)]--;
        }

        buckets[BUCKETS_INDEX4((fast_uint_t)c0, (s << 1) & 3)]--;
    }

    {
        sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];

        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
        { 
            sum += buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 3)]; temp_bucket[j] = sum;
        }

        return sum;
    }
}

static void libsais16x64_radix_sort_lms_suffixes_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais16x64_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 0]]);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 1]]);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 2]]);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 3]]);

        sa_sint_t p0 = SA[i - 0]; SA[--induction_bucket[BUCKETS_INDEX2(T[p0], 0)]] = p0;
        sa_sint_t p1 = SA[i - 1]; SA[--induction_bucket[BUCKETS_INDEX2(T[p1], 0)]] = p1;
        sa_sint_t p2 = SA[i - 2]; SA[--induction_bucket[BUCKETS_INDEX2(T[p2], 0)]] = p2;
        sa_sint_t p3 = SA[i - 3]; SA[--induction_bucket[BUCKETS_INDEX2(T[p3], 0)]] = p3;
    }

    for (j -= prefetch_distance + 3; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[--induction_bucket[BUCKETS_INDEX2(T[p], 0)]] = p;
    }
}

static void libsais16x64_radix_sort_lms_suffixes_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_radix_sort_lms_suffixes_16u(T, SA, &buckets[4 * ALPHABET_SIZE], (fast_sint_t)n - (fast_sint_t)m + 1, (fast_sint_t)m - 1);
}

static void libsais16x64_radix_sort_lms_suffixes_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais16x64_prefetchr(&SA[i - 3 * prefetch_distance]);
        
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 0]]);
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 1]]);
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 2]]);
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 3]]);

        libsais16x64_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 0]]]);
        libsais16x64_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 1]]]);
        libsais16x64_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 2]]]);
        libsais16x64_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 3]]]);

        sa_sint_t p0 = SA[i - 0]; SA[--induction_bucket[T[p0]]] = p0;
        sa_sint_t p1 = SA[i - 1]; SA[--induction_bucket[T[p1]]] = p1;
        sa_sint_t p2 = SA[i - 2]; SA[--induction_bucket[T[p2]]] = p2;
        sa_sint_t p3 = SA[i - 3]; SA[--induction_bucket[T[p3]]] = p3;
    }

    for (j -= 2 * prefetch_distance + 3; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[--induction_bucket[T[p]]] = p;
    }
}

static void libsais16x64_radix_sort_lms_suffixes_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t * RESTRICT induction_bucket)
{
    libsais16x64_radix_sort_lms_suffixes_32s_6k(T, SA, induction_bucket, (fast_sint_t)n - (fast_sint_t)m + 1, (fast_sint_t)m - 1);
}

static sa_sint_t libsais16x64_radix_sort_lms_suffixes_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t             i = n - 2;
    sa_sint_t             m = 0;
    fast_uint_t           s = 1;
    fast_sint_t           c0 = T[n - 1];
    fast_sint_t           c1 = 0;
    fast_sint_t           c2 = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais16x64_prefetchr(&T[i - 2 * prefetch_distance]);

        libsais16x64_prefetchw(&buckets[T[i - prefetch_distance - 0]]);
        libsais16x64_prefetchw(&buckets[T[i - prefetch_distance - 1]]);
        libsais16x64_prefetchw(&buckets[T[i - prefetch_distance - 2]]);
        libsais16x64_prefetchw(&buckets[T[i - prefetch_distance - 3]]);

        c1 = T[i - 0]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c0]] = i + 1; m++; }
        
        c0 = T[i - 1]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c1]] = i - 0; m++; }

        c1 = T[i - 2]; s = (s << 1) + (fast_uint_t)(c1 > (c0 - (fast_sint_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c0]] = i - 1; m++; }

        c0 = T[i - 3]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c1]] = i - 2; m++; }
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (fast_uint_t)(c0 > (c1 - (fast_sint_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c1]] = i + 1; m++; }
    }

    if (m > 1)
    {
        SA[buckets[c2]] = 0;
    }

    return m;
}

static void libsais16x64_radix_sort_set_markers_32s_6k(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais16x64_prefetchr(&induction_bucket[i + 2 * prefetch_distance]);

        libsais16x64_prefetchw(&SA[induction_bucket[i + prefetch_distance + 0]]);
        libsais16x64_prefetchw(&SA[induction_bucket[i + prefetch_distance + 1]]);
        libsais16x64_prefetchw(&SA[induction_bucket[i + prefetch_distance + 2]]);
        libsais16x64_prefetchw(&SA[induction_bucket[i + prefetch_distance + 3]]);

        SA[induction_bucket[i + 0]] |= SAINT_MIN;
        SA[induction_bucket[i + 1]] |= SAINT_MIN;
        SA[induction_bucket[i + 2]] |= SAINT_MIN;
        SA[induction_bucket[i + 3]] |= SAINT_MIN;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[induction_bucket[i]] |= SAINT_MIN;
    }
}

static void libsais16x64_radix_sort_set_markers_32s_6k_omp(sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = (fast_sint_t)k - 1;

    libsais16x64_radix_sort_set_markers_32s_6k(SA, induction_bucket, omp_block_start, omp_block_size);
}

static void libsais16x64_initialize_buckets_for_partial_sorting_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count)
{
    sa_sint_t * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

    buckets[BUCKETS_INDEX4((fast_uint_t)T[first_lms_suffix], 1)]++;

    fast_sint_t i, j; sa_sint_t sum0 = left_suffixes_count + 1, sum1 = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    { 
        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;

        sum0 += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 2)];
        sum1 += buckets[i + BUCKETS_INDEX4(0, 1)];

        buckets[j + BUCKETS_INDEX2(0, 0)] = sum0;
        buckets[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static void libsais16x64_initialize_buckets_for_partial_sorting_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count)
{
    sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];

    fast_sint_t i, j; sa_sint_t sum0 = left_suffixes_count + 1, sum1 = 0, sum2 = 0;
    for (first_lms_suffix = T[first_lms_suffix], i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4((fast_sint_t)first_lms_suffix - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    {
        sa_sint_t SS = buckets[i + BUCKETS_INDEX4(0, 0)];
        sa_sint_t LS = buckets[i + BUCKETS_INDEX4(0, 1)];
        sa_sint_t SL = buckets[i + BUCKETS_INDEX4(0, 2)];
        sa_sint_t LL = buckets[i + BUCKETS_INDEX4(0, 3)];

        buckets[i + BUCKETS_INDEX4(0, 0)] = sum0;
        buckets[i + BUCKETS_INDEX4(0, 1)] = sum2;
        buckets[i + BUCKETS_INDEX4(0, 2)] = 0;
        buckets[i + BUCKETS_INDEX4(0, 3)] = 0;

        sum0 += SS + SL; sum1 += LS; sum2 += LS + LL;

        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;
        temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }

    for (sum1 += 1; i <= BUCKETS_INDEX4((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    { 
        sa_sint_t SS = buckets[i + BUCKETS_INDEX4(0, 0)];
        sa_sint_t LS = buckets[i + BUCKETS_INDEX4(0, 1)];
        sa_sint_t SL = buckets[i + BUCKETS_INDEX4(0, 2)];
        sa_sint_t LL = buckets[i + BUCKETS_INDEX4(0, 3)];

        buckets[i + BUCKETS_INDEX4(0, 0)] = sum0;
        buckets[i + BUCKETS_INDEX4(0, 1)] = sum2;
        buckets[i + BUCKETS_INDEX4(0, 2)] = 0;
        buckets[i + BUCKETS_INDEX4(0, 3)] = 0;

        sum0 += SS + SL; sum1 += LS; sum2 += LS + LL;

        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;
        temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static sa_sint_t libsais16x64_partial_sorting_scan_left_to_right_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais16x64_prefetchr(&T[SA[i + prefetch_distance + 0] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i + prefetch_distance + 0] & SAINT_MAX] - 2);
        libsais16x64_prefetchr(&T[SA[i + prefetch_distance + 1] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i + prefetch_distance + 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i + 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] >= T[p0 - 1]);
        SA[induction_bucket[v0]++] = (p0 - 1) | ((sa_sint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d;

        sa_sint_t p1 = SA[i + 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] >= T[p1 - 1]);
        SA[induction_bucket[v1]++] = (p1 - 1) | ((sa_sint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] >= T[p - 1]);
        SA[induction_bucket[v]++] = (p - 1) | ((sa_sint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static sa_sint_t libsais16x64_partial_sorting_scan_left_to_right_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t left_suffixes_count, sa_sint_t d)
{
    sa_sint_t * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    SA[induction_bucket[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])]++] = (n - 1) | SAINT_MIN;
    distinct_names[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])] = ++d;

    d = libsais16x64_partial_sorting_scan_left_to_right_16u(T, SA, buckets, d, 0, left_suffixes_count);

    return d;
}

static sa_sint_t libsais16x64_partial_sorting_scan_left_to_right_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchr(&SA[i + 3 * prefetch_distance]);

        libsais16x64_prefetchr(&T[SA[i + 2 * prefetch_distance + 0] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i + 2 * prefetch_distance + 0] & SAINT_MAX] - 2);
        libsais16x64_prefetchr(&T[SA[i + 2 * prefetch_distance + 1] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i + 2 * prefetch_distance + 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i + prefetch_distance + 0] & SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX4(T[p0 - (p0 > 0)], 0); libsais16x64_prefetchw(&buckets[v0]);
        sa_sint_t p1 = SA[i + prefetch_distance + 1] & SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX4(T[p1 - (p1 > 0)], 0); libsais16x64_prefetchw(&buckets[v1]);

        sa_sint_t p2 = SA[i + 0]; d += (p2 < 0); p2 &= SAINT_MAX; sa_sint_t v2 = BUCKETS_INDEX4(T[p2 - 1], T[p2 - 2] >= T[p2 - 1]);
        SA[buckets[v2]++] = (p2 - 1) | ((sa_sint_t)(buckets[2 + v2] != d) << (SAINT_BIT - 1)); buckets[2 + v2] = d;

        sa_sint_t p3 = SA[i + 1]; d += (p3 < 0); p3 &= SAINT_MAX; sa_sint_t v3 = BUCKETS_INDEX4(T[p3 - 1], T[p3 - 2] >= T[p3 - 1]);
        SA[buckets[v3]++] = (p3 - 1) | ((sa_sint_t)(buckets[2 + v3] != d) << (SAINT_BIT - 1)); buckets[2 + v3] = d;
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX4(T[p - 1], T[p - 2] >= T[p - 1]);
        SA[buckets[v]++] = (p - 1) | ((sa_sint_t)(buckets[2 + v] != d) << (SAINT_BIT - 1)); buckets[2 + v] = d;
    }

    return d;
}

static void libsais16x64_partial_sorting_scan_left_to_right_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchw(&SA[i + 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i + 2 * prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i + 2 * prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);
        sa_sint_t s2 = SA[i + 1 * prefetch_distance + 0]; if (s2 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s2 - 1]]); libsais16x64_prefetchr(&T[s2] - 2); }
        sa_sint_t s3 = SA[i + 1 * prefetch_distance + 1]; if (s3 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s3 - 1]]); libsais16x64_prefetchr(&T[s3] - 2); }

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 & SAINT_MAX; if (p0 > 0) { SA[i + 0] = 0; SA[induction_bucket[T[p0 - 1]]++] = (p0 - 1) | ((sa_sint_t)(T[p0 - 2] < T[p0 - 1]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 & SAINT_MAX; if (p1 > 0) { SA[i + 1] = 0; SA[induction_bucket[T[p1 - 1]]++] = (p1 - 1) | ((sa_sint_t)(T[p1 - 2] < T[p1 - 1]) << (SAINT_BIT - 1)); }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p & SAINT_MAX; if (p > 0) { SA[i] = 0; SA[induction_bucket[T[p - 1]]++] = (p - 1) | ((sa_sint_t)(T[p - 2] < T[p - 1]) << (SAINT_BIT - 1)); }
    }
}

static sa_sint_t libsais16x64_partial_sorting_scan_left_to_right_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t left_suffixes_count, sa_sint_t d)
{
    SA[buckets[BUCKETS_INDEX4(T[n - 1], T[n - 2] >= T[n - 1])]++] = (n - 1) | SAINT_MIN;
    buckets[2 + BUCKETS_INDEX4(T[n - 1], T[n - 2] >= T[n - 1])] = ++d;

    d = libsais16x64_partial_sorting_scan_left_to_right_32s_6k(T, SA, buckets, d, 0, left_suffixes_count);

    return d;
}

static void libsais16x64_partial_sorting_scan_left_to_right_32s_1k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets)
{
    SA[buckets[T[n - 1]]++] = (n - 1) | ((sa_sint_t)(T[n - 2] < T[n - 1]) << (SAINT_BIT - 1));

    libsais16x64_partial_sorting_scan_left_to_right_32s_1k(T, SA, buckets, 0, n);
}

static void libsais16x64_partial_sorting_shift_markers_16u_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, const sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 32;

    const sa_sint_t * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

    fast_sint_t c;

    for (c = BUCKETS_INDEX2(ALPHABET_SIZE - 1, 0); c >= BUCKETS_INDEX2(1, 0); c -= BUCKETS_INDEX2(1, 0))
    {
        fast_sint_t i, j; sa_sint_t s = SAINT_MIN;
        for (i = (fast_sint_t)temp_bucket[c] - 1, j = (fast_sint_t)buckets[c - BUCKETS_INDEX2(1, 0)] + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchw(&SA[i - prefetch_distance]);

            sa_sint_t p0 = SA[i - 0], q0 = (p0 & SAINT_MIN) ^ s; s = s ^ q0; SA[i - 0] = p0 ^ q0;
            sa_sint_t p1 = SA[i - 1], q1 = (p1 & SAINT_MIN) ^ s; s = s ^ q1; SA[i - 1] = p1 ^ q1;
            sa_sint_t p2 = SA[i - 2], q2 = (p2 & SAINT_MIN) ^ s; s = s ^ q2; SA[i - 2] = p2 ^ q2;
            sa_sint_t p3 = SA[i - 3], q3 = (p3 & SAINT_MIN) ^ s; s = s ^ q3; SA[i - 3] = p3 ^ q3;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            sa_sint_t p = SA[i], q = (p & SAINT_MIN) ^ s; s = s ^ q; SA[i] = p ^ q;
        }
    }
}

static void libsais16x64_partial_sorting_shift_markers_32s_6k_omp(sa_sint_t * RESTRICT SA, sa_sint_t k, const sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 32;

    const sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];
    
    fast_sint_t c;

    for (c = (fast_sint_t)k - 1; c >= 1; c -= 1)
    {
        fast_sint_t i, j; sa_sint_t s = SAINT_MIN;
        for (i = (fast_sint_t)buckets[BUCKETS_INDEX4(c, 0)] - 1, j = (fast_sint_t)temp_bucket[BUCKETS_INDEX2(c - 1, 0)] + 3; i >= j; i -= 4)
        {
            libsais16x64_prefetchw(&SA[i - prefetch_distance]);

            sa_sint_t p0 = SA[i - 0], q0 = (p0 & SAINT_MIN) ^ s; s = s ^ q0; SA[i - 0] = p0 ^ q0;
            sa_sint_t p1 = SA[i - 1], q1 = (p1 & SAINT_MIN) ^ s; s = s ^ q1; SA[i - 1] = p1 ^ q1;
            sa_sint_t p2 = SA[i - 2], q2 = (p2 & SAINT_MIN) ^ s; s = s ^ q2; SA[i - 2] = p2 ^ q2;
            sa_sint_t p3 = SA[i - 3], q3 = (p3 & SAINT_MIN) ^ s; s = s ^ q3; SA[i - 3] = p3 ^ q3;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            sa_sint_t p = SA[i], q = (p & SAINT_MIN) ^ s; s = s ^ q; SA[i] = p ^ q;
        }
    }
}

static void libsais16x64_partial_sorting_shift_buckets_32s_6k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];

    fast_sint_t i;
    for (i = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0))
    {
        buckets[2 * i + BUCKETS_INDEX4(0, 0)] = temp_bucket[i + BUCKETS_INDEX2(0, 0)];
        buckets[2 * i + BUCKETS_INDEX4(0, 1)] = temp_bucket[i + BUCKETS_INDEX2(0, 1)];
    }
}

static sa_sint_t libsais16x64_partial_sorting_scan_right_to_left_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 2);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i - 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]);
        SA[--induction_bucket[v0]] = (p0 - 1) | ((sa_sint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d;

        sa_sint_t p1 = SA[i - 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]);
        SA[--induction_bucket[v1]] = (p1 - 1) | ((sa_sint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
        SA[--induction_bucket[v]] = (p - 1) | ((sa_sint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static void libsais16x64_partial_sorting_scan_right_to_left_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t d)
{
    fast_sint_t scan_start    = (fast_sint_t)left_suffixes_count + 1;
    fast_sint_t scan_end      = (fast_sint_t)n - (fast_sint_t)first_lms_suffix;

    libsais16x64_partial_sorting_scan_right_to_left_16u(T, SA, buckets, d, scan_start, scan_end - scan_start);
}

static sa_sint_t libsais16x64_partial_sorting_scan_right_to_left_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchr(&SA[i - 3 * prefetch_distance]);

        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 0] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 0] & SAINT_MAX] - 2);
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 1] & SAINT_MAX] - 1);
        libsais16x64_prefetchr(&T[SA[i - 2 * prefetch_distance - 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i - prefetch_distance - 0] & SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX4(T[p0 - (p0 > 0)], 0); libsais16x64_prefetchw(&buckets[v0]);
        sa_sint_t p1 = SA[i - prefetch_distance - 1] & SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX4(T[p1 - (p1 > 0)], 0); libsais16x64_prefetchw(&buckets[v1]);

        sa_sint_t p2 = SA[i - 0]; d += (p2 < 0); p2 &= SAINT_MAX; sa_sint_t v2 = BUCKETS_INDEX4(T[p2 - 1], T[p2 - 2] > T[p2 - 1]);
        SA[--buckets[v2]] = (p2 - 1) | ((sa_sint_t)(buckets[2 + v2] != d) << (SAINT_BIT - 1)); buckets[2 + v2] = d;

        sa_sint_t p3 = SA[i - 1]; d += (p3 < 0); p3 &= SAINT_MAX; sa_sint_t v3 = BUCKETS_INDEX4(T[p3 - 1], T[p3 - 2] > T[p3 - 1]);
        SA[--buckets[v3]] = (p3 - 1) | ((sa_sint_t)(buckets[2 + v3] != d) << (SAINT_BIT - 1)); buckets[2 + v3] = d;
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX4(T[p - 1], T[p - 2] > T[p - 1]);
        SA[--buckets[v]] = (p - 1) | ((sa_sint_t)(buckets[2 + v] != d) << (SAINT_BIT - 1)); buckets[2 + v] = d;
    }

    return d;
}

static void libsais16x64_partial_sorting_scan_right_to_left_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchw(&SA[i - 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i - 2 * prefetch_distance - 0]; const sa_sint_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i - 2 * prefetch_distance - 1]; const sa_sint_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);
        sa_sint_t s2 = SA[i - 1 * prefetch_distance - 0]; if (s2 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s2 - 1]]); libsais16x64_prefetchr(&T[s2] - 2); }
        sa_sint_t s3 = SA[i - 1 * prefetch_distance - 1]; if (s3 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s3 - 1]]); libsais16x64_prefetchr(&T[s3] - 2); }

        sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { SA[i - 0] = 0; SA[--induction_bucket[T[p0 - 1]]] = (p0 - 1) | ((sa_sint_t)(T[p0 - 2] > T[p0 - 1]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { SA[i - 1] = 0; SA[--induction_bucket[T[p1 - 1]]] = (p1 - 1) | ((sa_sint_t)(T[p1 - 2] > T[p1 - 1]) << (SAINT_BIT - 1)); }
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { SA[i] = 0; SA[--induction_bucket[T[p - 1]]] = (p - 1) | ((sa_sint_t)(T[p - 2] > T[p - 1]) << (SAINT_BIT - 1)); }
    }
}

static sa_sint_t libsais16x64_partial_sorting_scan_right_to_left_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t d)
{
    fast_sint_t scan_start    = (fast_sint_t)left_suffixes_count + 1;
    fast_sint_t scan_end      = (fast_sint_t)n - (fast_sint_t)first_lms_suffix;

    d = libsais16x64_partial_sorting_scan_right_to_left_32s_6k(T, SA, buckets, d, scan_start, scan_end - scan_start);

    return d;
}

static void libsais16x64_partial_sorting_scan_right_to_left_32s_1k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_partial_sorting_scan_right_to_left_32s_1k(T, SA, buckets, 0, n);
}

static fast_sint_t libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k(sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j, l;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3, l = omp_block_start; i < j; i += 4)
    {
        libsais16x64_prefetchr(&SA[i + prefetch_distance]);

        sa_sint_t s0 = SA[i + 0]; SA[l] = s0 & SAINT_MAX; l += (s0 < 0);
        sa_sint_t s1 = SA[i + 1]; SA[l] = s1 & SAINT_MAX; l += (s1 < 0);
        sa_sint_t s2 = SA[i + 2]; SA[l] = s2 & SAINT_MAX; l += (s2 < 0);
        sa_sint_t s3 = SA[i + 3]; SA[l] = s3 & SAINT_MAX; l += (s3 < 0);
    }

    for (j += 3; i < j; i += 1)
    {
        sa_sint_t s = SA[i]; SA[l] = s & SAINT_MAX; l += (s < 0);
    }

    return l;
}

static void libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k_omp(sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = n - omp_block_start;

    libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k(SA, omp_block_start, omp_block_size);
}

static void libsais16x64_induce_partial_order_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count)
{
    memset(&buckets[2 * ALPHABET_SIZE], 0, (size_t)2 * ALPHABET_SIZE * sizeof(sa_sint_t));

    sa_sint_t d = libsais16x64_partial_sorting_scan_left_to_right_16u_omp(T, SA, n, k, buckets, left_suffixes_count, 0);
    libsais16x64_partial_sorting_shift_markers_16u_omp(SA, n, buckets);
    libsais16x64_partial_sorting_scan_right_to_left_16u_omp(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count, d);
}

static void libsais16x64_induce_partial_order_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count)
{
    sa_sint_t d = libsais16x64_partial_sorting_scan_left_to_right_32s_6k_omp(T, SA, n, buckets, left_suffixes_count, 0);
    libsais16x64_partial_sorting_shift_markers_32s_6k_omp(SA, k, buckets);
    libsais16x64_partial_sorting_shift_buckets_32s_6k(k, buckets);
    libsais16x64_partial_sorting_scan_right_to_left_32s_6k_omp(T, SA, n, buckets, first_lms_suffix, left_suffixes_count, d);
}

static void libsais16x64_induce_partial_order_32s_1k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_start_32s_1k(k, buckets);
    libsais16x64_partial_sorting_scan_left_to_right_32s_1k_omp(T, SA, n, buckets);

    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_end_32s_1k(k, buckets);
    libsais16x64_partial_sorting_scan_right_to_left_32s_1k_omp(T, SA, n, buckets);

    libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k_omp(SA, n);
}

static sa_sint_t libsais16x64_renumber_lms_suffixes_16u(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t name, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais16x64_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 0] & SAINT_MAX) >> 1]);
        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 1] & SAINT_MAX) >> 1]);
        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 2] & SAINT_MAX) >> 1]);
        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 3] & SAINT_MAX) >> 1]);

        sa_sint_t p0 = SA[i + 0]; SAm[(p0 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p0 < 0;
        sa_sint_t p1 = SA[i + 1]; SAm[(p1 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p1 < 0;
        sa_sint_t p2 = SA[i + 2]; SAm[(p2 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p2 < 0;
        sa_sint_t p3 = SA[i + 3]; SAm[(p3 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p3 < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SAm[(p & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p < 0;
    }

    return name;
}

static fast_sint_t libsais16x64_gather_marked_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t l, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    l -= 1;

    fast_sint_t i, j;
    for (i = (fast_sint_t)m + omp_block_start + omp_block_size - 1, j = (fast_sint_t)m + omp_block_start + 3; i >= j; i -= 4)
    {
        libsais16x64_prefetchr(&SA[i - prefetch_distance]);

        sa_sint_t s0 = SA[i - 0]; SA[l] = s0 & SAINT_MAX; l -= s0 < 0;
        sa_sint_t s1 = SA[i - 1]; SA[l] = s1 & SAINT_MAX; l -= s1 < 0;
        sa_sint_t s2 = SA[i - 2]; SA[l] = s2 & SAINT_MAX; l -= s2 < 0;
        sa_sint_t s3 = SA[i - 3]; SA[l] = s3 & SAINT_MAX; l -= s3 < 0;
    }

    for (j -= 3; i >= j; i -= 1)
    {
        sa_sint_t s = SA[i]; SA[l] = s & SAINT_MAX; l -= s < 0;
    }

    l += 1;

    return l;
}

static sa_sint_t libsais16x64_renumber_lms_suffixes_16u_omp(sa_sint_t * RESTRICT SA, sa_sint_t m)
{
    sa_sint_t name = 0;

    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = m - omp_block_start;

    name = libsais16x64_renumber_lms_suffixes_16u(SA, m, 0, omp_block_start, omp_block_size);

    return name;
}

static void libsais16x64_gather_marked_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = ((fast_sint_t)n >> 1) - omp_block_start;

    libsais16x64_gather_marked_lms_suffixes(SA, m, (fast_sint_t)n + (fast_sint_t)fs, omp_block_start, omp_block_size);
}

static sa_sint_t libsais16x64_renumber_and_gather_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs)
{
    memset(&SA[m], 0, ((size_t)n >> 1) * sizeof(sa_sint_t));

    sa_sint_t name = libsais16x64_renumber_lms_suffixes_16u_omp(SA, m);
    if (name < m)
    {
        libsais16x64_gather_marked_lms_suffixes_omp(SA, n, m, fs);
    }
    else
    {
        fast_sint_t i; for (i = 0; i < m; i += 1) { SA[i] &= SAINT_MAX; }
    }

    return name;
}

static sa_sint_t libsais16x64_renumber_distinct_lms_suffixes_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t name, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j; sa_sint_t p0, p1, p2, p3 = 0;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais16x64_prefetchw(&SA[i + 2 * prefetch_distance]);

        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 0] & SAINT_MAX) >> 1]);
        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 1] & SAINT_MAX) >> 1]);
        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 2] & SAINT_MAX) >> 1]);
        libsais16x64_prefetchw(&SAm[(SA[i + prefetch_distance + 3] & SAINT_MAX) >> 1]);

        p0 = SA[i + 0]; SAm[(SA[i + 0] = p0 & SAINT_MAX) >> 1] = name | (p0 & p3 & SAINT_MIN); name += p0 < 0;
        p1 = SA[i + 1]; SAm[(SA[i + 1] = p1 & SAINT_MAX) >> 1] = name | (p1 & p0 & SAINT_MIN); name += p1 < 0;
        p2 = SA[i + 2]; SAm[(SA[i + 2] = p2 & SAINT_MAX) >> 1] = name | (p2 & p1 & SAINT_MIN); name += p2 < 0;
        p3 = SA[i + 3]; SAm[(SA[i + 3] = p3 & SAINT_MAX) >> 1] = name | (p3 & p2 & SAINT_MIN); name += p3 < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        p2 = p3; p3 = SA[i]; SAm[(SA[i] = p3 & SAINT_MAX) >> 1] = name | (p3 & p2 & SAINT_MIN); name += p3 < 0;
    }

    return name;
}

static void libsais16x64_mark_distinct_lms_suffixes_32s(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j; sa_sint_t p0, p1, p2, p3 = 0;
    for (i = (fast_sint_t)m + omp_block_start, j = (fast_sint_t)m + omp_block_start + omp_block_size - 3; i < j; i += 4)
    {
        libsais16x64_prefetchw(&SA[i + prefetch_distance]);

        p0 = SA[i + 0]; SA[i + 0] = p0 & (p3 | SAINT_MAX); p0 = (p0 == 0) ? p3 : p0;
        p1 = SA[i + 1]; SA[i + 1] = p1 & (p0 | SAINT_MAX); p1 = (p1 == 0) ? p0 : p1;
        p2 = SA[i + 2]; SA[i + 2] = p2 & (p1 | SAINT_MAX); p2 = (p2 == 0) ? p1 : p2;
        p3 = SA[i + 3]; SA[i + 3] = p3 & (p2 | SAINT_MAX); p3 = (p3 == 0) ? p2 : p3;
    }

    for (j += 3; i < j; i += 1)
    {
        p2 = p3; p3 = SA[i]; SA[i] = p3 & (p2 | SAINT_MAX); p3 = (p3 == 0) ? p2 : p3;
    }
}

static void libsais16x64_clamp_lms_suffixes_length_32s(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3; i < j; i += 4)
    {
        libsais16x64_prefetchw(&SAm[i + prefetch_distance]);

        SAm[i + 0] = (SAm[i + 0] < 0 ? SAm[i + 0] : 0) & SAINT_MAX;
        SAm[i + 1] = (SAm[i + 1] < 0 ? SAm[i + 1] : 0) & SAINT_MAX;
        SAm[i + 2] = (SAm[i + 2] < 0 ? SAm[i + 2] : 0) & SAINT_MAX;
        SAm[i + 3] = (SAm[i + 3] < 0 ? SAm[i + 3] : 0) & SAINT_MAX;
    }

    for (j += 3; i < j; i += 1)
    {
        SAm[i] = (SAm[i] < 0 ? SAm[i] : 0) & SAINT_MAX;
    }
}

static sa_sint_t libsais16x64_renumber_distinct_lms_suffixes_32s_4k_omp(sa_sint_t * RESTRICT SA, sa_sint_t m)
{
    sa_sint_t name = 0;

    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = m - omp_block_start;

    name = libsais16x64_renumber_distinct_lms_suffixes_32s_4k(SA, m, 1, omp_block_start, omp_block_size);

    return name - 1;
}

static void libsais16x64_mark_distinct_lms_suffixes_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = (fast_sint_t)n >> 1;

    libsais16x64_mark_distinct_lms_suffixes_32s(SA, m, omp_block_start, omp_block_size);
}

static void libsais16x64_clamp_lms_suffixes_length_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = (fast_sint_t)n >> 1;

    libsais16x64_clamp_lms_suffixes_length_32s(SA, m, omp_block_start, omp_block_size);
}

static sa_sint_t libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    memset(&SA[m], 0, ((size_t)n >> 1) * sizeof(sa_sint_t));

    sa_sint_t name = libsais16x64_renumber_distinct_lms_suffixes_32s_4k_omp(SA, m);
    if (name < m)
    {
        libsais16x64_mark_distinct_lms_suffixes_32s_omp(SA, n, m);
    }

    return name;
}

static sa_sint_t libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT SAm = &SA[m];

    {
        libsais16x64_gather_lms_suffixes_32s(T, SA, n);

        memset(&SA[m], 0, ((size_t)n - (size_t)m - (size_t)m) * sizeof(sa_sint_t));

        fast_sint_t i, j;
        for (i = (fast_sint_t)n - (fast_sint_t)m, j = (fast_sint_t)n - 1 - prefetch_distance - 3; i < j; i += 4)
        {
            libsais16x64_prefetchr(&SA[i + 2 * prefetch_distance]);

            libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 0]) >> 1]);
            libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 1]) >> 1]);
            libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 2]) >> 1]);
            libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 3]) >> 1]);

            SAm[((sa_uint_t)SA[i + 0]) >> 1] = SA[i + 1] - SA[i + 0] + 1 + SAINT_MIN;
            SAm[((sa_uint_t)SA[i + 1]) >> 1] = SA[i + 2] - SA[i + 1] + 1 + SAINT_MIN;
            SAm[((sa_uint_t)SA[i + 2]) >> 1] = SA[i + 3] - SA[i + 2] + 1 + SAINT_MIN;
            SAm[((sa_uint_t)SA[i + 3]) >> 1] = SA[i + 4] - SA[i + 3] + 1 + SAINT_MIN;
        }

        for (j += prefetch_distance + 3; i < j; i += 1)
        {
            SAm[((sa_uint_t)SA[i]) >> 1] = SA[i + 1] - SA[i] + 1 + SAINT_MIN;
        }

        SAm[((sa_uint_t)SA[n - 1]) >> 1] = 1 + SAINT_MIN;
    }

    {
        libsais16x64_clamp_lms_suffixes_length_32s_omp(SA, n, m);
    }

    sa_sint_t name = 1;

    {
        fast_sint_t i, j, p = SA[0], plen = SAm[p >> 1]; sa_sint_t pdiff = SAINT_MIN;
        for (i = 1, j = m - prefetch_distance - 1; i < j; i += 2)
        {
            libsais16x64_prefetchr(&SA[i + 2 * prefetch_distance]);
            
            libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 0]) >> 1]); libsais16x64_prefetchr(&T[((sa_uint_t)SA[i + prefetch_distance + 0])]);
            libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 1]) >> 1]); libsais16x64_prefetchr(&T[((sa_uint_t)SA[i + prefetch_distance + 1])]);

            fast_sint_t q = SA[i + 0], qlen = SAm[q >> 1]; sa_sint_t qdiff = SAINT_MIN;
            if (plen == qlen) { fast_sint_t l = 0; do { if (T[p + l] != T[q + l]) { break; } } while (++l < qlen); qdiff = (sa_sint_t)(l - qlen) & SAINT_MIN; }
            SAm[p >> 1] = name | (pdiff & qdiff); name += (qdiff < 0);

            p = SA[i + 1]; plen = SAm[p >> 1]; pdiff = SAINT_MIN;
            if (qlen == plen) { fast_sint_t l = 0; do { if (T[q + l] != T[p + l]) { break; } } while (++l < plen); pdiff = (sa_sint_t)(l - plen) & SAINT_MIN; }
            SAm[q >> 1] = name | (qdiff & pdiff); name += (pdiff < 0);
        }

        for (j += prefetch_distance + 1; i < j; i += 1)
        {
            fast_sint_t q = SA[i], qlen = SAm[q >> 1]; sa_sint_t qdiff = SAINT_MIN;
            if (plen == qlen) { fast_sint_t l = 0; do { if (T[p + l] != T[q + l]) { break; } } while (++l < plen); qdiff = (sa_sint_t)(l - plen) & SAINT_MIN; }
            SAm[p >> 1] = name | (pdiff & qdiff); name += (qdiff < 0);

            p = q; plen = qlen; pdiff = qdiff;
        }

        SAm[p >> 1] = name | pdiff; name++;
    }

    if (name <= m)
    {
        libsais16x64_mark_distinct_lms_suffixes_32s_omp(SA, n, m);
    }

    return name - 1;
}

static void libsais16x64_reconstruct_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    const sa_sint_t * RESTRICT SAnm = &SA[n - m];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais16x64_prefetchw(&SA[i + 2 * prefetch_distance]);

        libsais16x64_prefetchr(&SAnm[SA[i + prefetch_distance + 0]]);
        libsais16x64_prefetchr(&SAnm[SA[i + prefetch_distance + 1]]);
        libsais16x64_prefetchr(&SAnm[SA[i + prefetch_distance + 2]]);
        libsais16x64_prefetchr(&SAnm[SA[i + prefetch_distance + 3]]);

        SA[i + 0] = SAnm[SA[i + 0]];
        SA[i + 1] = SAnm[SA[i + 1]];
        SA[i + 2] = SAnm[SA[i + 2]];
        SA[i + 3] = SAnm[SA[i + 3]];
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[i] = SAnm[SA[i]];
    }
}

static void libsais16x64_reconstruct_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = m;

    libsais16x64_reconstruct_lms_suffixes(SA, n, m, omp_block_start, omp_block_size);
}

static void libsais16x64_place_lms_suffixes_interval_16u(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    const sa_sint_t * RESTRICT bucket_end = &buckets[7 * ALPHABET_SIZE];

    fast_sint_t c, j = n;
    for (c = ALPHABET_SIZE - 2; c >= 0; --c)
    {
        fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(c, 1) + BUCKETS_INDEX2(1, 0)] - (fast_sint_t)buckets[BUCKETS_INDEX2(c, 1)];
        if (l > 0)
        {
            fast_sint_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(sa_sint_t));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (sa_sint_t)l], (size_t)l * sizeof(sa_sint_t));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais16x64_place_lms_suffixes_interval_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t m, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t c = k - 1; fast_sint_t i, l = buckets[c];
    for (i = (fast_sint_t)m - 1; i >= prefetch_distance + 3; i -= 4)
    {
        libsais16x64_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 0]]);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 1]]);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 2]]);
        libsais16x64_prefetchr(&T[SA[i - prefetch_distance - 3]]);

        sa_sint_t p0 = SA[i - 0]; if (T[p0] != c) { c = T[p0]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(sa_sint_t)); l = buckets[c]; } SA[--l] = p0;
        sa_sint_t p1 = SA[i - 1]; if (T[p1] != c) { c = T[p1]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(sa_sint_t)); l = buckets[c]; } SA[--l] = p1;
        sa_sint_t p2 = SA[i - 2]; if (T[p2] != c) { c = T[p2]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(sa_sint_t)); l = buckets[c]; } SA[--l] = p2;
        sa_sint_t p3 = SA[i - 3]; if (T[p3] != c) { c = T[p3]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(sa_sint_t)); l = buckets[c]; } SA[--l] = p3;
    }

    for (; i >= 0; i -= 1)
    {
        sa_sint_t p = SA[i]; if (T[p] != c) { c = T[p]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(sa_sint_t)); l = buckets[c]; } SA[--l] = p;
    }

    memset(&SA[0], 0, (size_t)l * sizeof(sa_sint_t));
}

static void libsais16x64_place_lms_suffixes_histogram_32s_6k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    const sa_sint_t * RESTRICT bucket_end = &buckets[5 * (fast_sint_t)k];

    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; --c)
    {
        fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX4(c, 1)];
        if (l > 0)
        {
            fast_sint_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(sa_sint_t));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (sa_sint_t)l], (size_t)l * sizeof(sa_sint_t));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais16x64_place_lms_suffixes_histogram_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    const sa_sint_t * RESTRICT bucket_end = &buckets[3 * (fast_sint_t)k];

    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; --c)
    {
        fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(c, 1)];
        if (l > 0)
        {
            fast_sint_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(sa_sint_t));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (sa_sint_t)l], (size_t)l * sizeof(sa_sint_t));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais16x64_final_bwt_scan_left_to_right_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint16_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL); Ts0--; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint16_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL); Ts1--; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 & SAINT_MAX; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | SAINT_MIN; SA[induction_bucket[T[p0]]++] = p0 | ((sa_sint_t)(T[p0 - (p0 > 0)] < T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 & SAINT_MAX; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | SAINT_MIN; SA[induction_bucket[T[p1]]++] = p1 | ((sa_sint_t)(T[p1 - (p1 > 0)] < T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p & SAINT_MAX; if (p > 0) { p--; SA[i] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | ((sa_sint_t)(T[p - (p > 0)] < T[p]) << (SAINT_BIT - 1)); }
    }
}

static void libsais16x64_final_bwt_aux_scan_left_to_right_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint16_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL); Ts0--; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint16_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL); Ts1--; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 & SAINT_MAX; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | SAINT_MIN; SA[induction_bucket[T[p0]]++] = p0 | ((sa_sint_t)(T[p0 - (p0 > 0)] < T[p0]) << (SAINT_BIT - 1)); if ((p0 & rm) == 0) { I[p0 / (rm + 1)] = induction_bucket[T[p0]]; }}
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 & SAINT_MAX; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | SAINT_MIN; SA[induction_bucket[T[p1]]++] = p1 | ((sa_sint_t)(T[p1 - (p1 > 0)] < T[p1]) << (SAINT_BIT - 1)); if ((p1 & rm) == 0) { I[p1 / (rm + 1)] = induction_bucket[T[p1]]; }}
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p & SAINT_MAX; if (p > 0) { p--; SA[i] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | ((sa_sint_t)(T[p - (p > 0)] < T[p]) << (SAINT_BIT - 1)); if ((p & rm) == 0) { I[p / (rm + 1)] = induction_bucket[T[p]]; } }
    }
}

static void libsais16x64_final_sorting_scan_left_to_right_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint16_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL); Ts0--; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint16_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL); Ts1--; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 ^ SAINT_MIN; if (p0 > 0) { p0--; SA[induction_bucket[T[p0]]++] = p0 | ((sa_sint_t)(T[p0 - (p0 > 0)] < T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 ^ SAINT_MIN; if (p1 > 0) { p1--; SA[induction_bucket[T[p1]]++] = p1 | ((sa_sint_t)(T[p1 - (p1 > 0)] < T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p ^ SAINT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | ((sa_sint_t)(T[p - (p > 0)] < T[p]) << (SAINT_BIT - 1)); }
    }
}

static void libsais16x64_final_sorting_scan_left_to_right_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais16x64_prefetchw(&SA[i + 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i + 2 * prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i + 2 * prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);
        sa_sint_t s2 = SA[i + 1 * prefetch_distance + 0]; if (s2 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s2 - 1]]); libsais16x64_prefetchr(&T[s2] - 2); }
        sa_sint_t s3 = SA[i + 1 * prefetch_distance + 1]; if (s3 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s3 - 1]]); libsais16x64_prefetchr(&T[s3] - 2); }

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 ^ SAINT_MIN; if (p0 > 0) { p0--; SA[induction_bucket[T[p0]]++] = p0 | ((sa_sint_t)(T[p0 - (p0 > 0)] < T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 ^ SAINT_MIN; if (p1 > 0) { p1--; SA[induction_bucket[T[p1]]++] = p1 | ((sa_sint_t)(T[p1 - (p1 > 0)] < T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p ^ SAINT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | ((sa_sint_t)(T[p - (p > 0)] < T[p]) << (SAINT_BIT - 1)); }
    }
}

static void libsais16x64_final_bwt_scan_left_to_right_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket)
{
    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | ((sa_sint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    libsais16x64_final_bwt_scan_left_to_right_16u(T, SA, induction_bucket, 0, n);
}

static void libsais16x64_final_bwt_aux_scan_left_to_right_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket)
{
    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | ((sa_sint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    if ((((sa_sint_t)n - 1) & rm) == 0) { I[((sa_sint_t)n - 1) / (rm + 1)] = induction_bucket[T[(sa_sint_t)n - 1]]; }

    libsais16x64_final_bwt_aux_scan_left_to_right_16u(T, SA, rm, I, induction_bucket, 0, n);
}

static void libsais16x64_final_sorting_scan_left_to_right_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket)
{
    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | ((sa_sint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    libsais16x64_final_sorting_scan_left_to_right_16u(T, SA, induction_bucket, 0, n);
}

static void libsais16x64_final_sorting_scan_left_to_right_32s_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT induction_bucket)
{
    SA[induction_bucket[T[n - 1]]++] = (n - 1) | ((sa_sint_t)(T[n - 2] < T[n - 1]) << (SAINT_BIT - 1));

    libsais16x64_final_sorting_scan_left_to_right_32s(T, SA, induction_bucket, 0, n);
}

static sa_sint_t libsais16x64_final_bwt_scan_right_to_left_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j; sa_sint_t index = -1;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint16_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL); Ts0--; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint16_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL); Ts1--; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);

        sa_sint_t p0 = SA[i - 0]; index = (p0 == 0) ? (sa_sint_t)(i - 0) : index;
        SA[i - 0] = p0 & SAINT_MAX; if (p0 > 0) { p0--; uint16_t c0 = T[p0 - (p0 > 0)], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p0 : t; }

        sa_sint_t p1 = SA[i - 1]; index = (p1 == 0) ? (sa_sint_t)(i - 1) : index;
        SA[i - 1] = p1 & SAINT_MAX; if (p1 > 0) { p1--; uint16_t c0 = T[p1 - (p1 > 0)], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p1 : t; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; index = (p == 0) ? (sa_sint_t)i : index;
        SA[i] = p & SAINT_MAX; if (p > 0) { p--; uint16_t c0 = T[p - (p > 0)], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; }
    }

    return index;
}

static void libsais16x64_final_bwt_aux_scan_right_to_left_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint16_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL); Ts0--; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint16_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL); Ts1--; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);

        sa_sint_t p0 = SA[i - 0];
        SA[i - 0] = p0 & SAINT_MAX; if (p0 > 0) { p0--; uint16_t c0 = T[p0 - (p0 > 0)], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p0 : t; if ((p0 & rm) == 0) { I[p0 / (rm + 1)] = induction_bucket[T[p0]] + 1; } }

        sa_sint_t p1 = SA[i - 1];
        SA[i - 1] = p1 & SAINT_MAX; if (p1 > 0) { p1--; uint16_t c0 = T[p1 - (p1 > 0)], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p1 : t; if ((p1 & rm) == 0) { I[p1 / (rm + 1)] = induction_bucket[T[p1]] + 1; } }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i];
        SA[i] = p & SAINT_MAX; if (p > 0) { p--; uint16_t c0 = T[p - (p > 0)], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; if ((p & rm) == 0) { I[p / (rm + 1)] = induction_bucket[T[p]] + 1; } }
    }
}

static void libsais16x64_final_sorting_scan_right_to_left_16u(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint16_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL); Ts0--; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint16_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL); Ts1--; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);

        sa_sint_t p0 = SA[i - 0]; SA[i - 0] = p0 & SAINT_MAX; if (p0 > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | ((sa_sint_t)(T[p0 - (p0 > 0)] > T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i - 1]; SA[i - 1] = p1 & SAINT_MAX; if (p1 > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | ((sa_sint_t)(T[p1 - (p1 > 0)] > T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p & SAINT_MAX; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | ((sa_sint_t)(T[p - (p > 0)] > T[p]) << (SAINT_BIT - 1)); }
    }
}

static void libsais16x64_final_sorting_scan_right_to_left_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais16x64_prefetchw(&SA[i - 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i - 2 * prefetch_distance - 0]; const sa_sint_t * Ts0 = &T[s0] - 1; libsais16x64_prefetchr(s0 > 0 ? Ts0 : NULL);
        sa_sint_t s1 = SA[i - 2 * prefetch_distance - 1]; const sa_sint_t * Ts1 = &T[s1] - 1; libsais16x64_prefetchr(s1 > 0 ? Ts1 : NULL);
        sa_sint_t s2 = SA[i - 1 * prefetch_distance - 0]; if (s2 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s2 - 1]]); libsais16x64_prefetchr(&T[s2] - 2); }
        sa_sint_t s3 = SA[i - 1 * prefetch_distance - 1]; if (s3 > 0) { libsais16x64_prefetchw(&induction_bucket[T[s3 - 1]]); libsais16x64_prefetchr(&T[s3] - 2); }

        sa_sint_t p0 = SA[i - 0]; SA[i - 0] = p0 & SAINT_MAX; if (p0 > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | ((sa_sint_t)(T[p0 - (p0 > 0)] > T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i - 1]; SA[i - 1] = p1 & SAINT_MAX; if (p1 > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | ((sa_sint_t)(T[p1 - (p1 > 0)] > T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p & SAINT_MAX; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | ((sa_sint_t)(T[p - (p > 0)] > T[p]) << (SAINT_BIT - 1)); }
    }
}

static sa_sint_t libsais16x64_final_bwt_scan_right_to_left_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket)
{
    sa_sint_t index = -1;

    index = libsais16x64_final_bwt_scan_right_to_left_16u(T, SA, induction_bucket, 0, n);

    return index;
}

static void libsais16x64_final_bwt_aux_scan_right_to_left_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket)
{
    libsais16x64_final_bwt_aux_scan_right_to_left_16u(T, SA, rm, I, induction_bucket, 0, n);
}

static void libsais16x64_final_sorting_scan_right_to_left_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket)
{
    libsais16x64_final_sorting_scan_right_to_left_16u(T, SA, induction_bucket, 0, n);
}

static void libsais16x64_final_sorting_scan_right_to_left_32s_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT induction_bucket)
{
    libsais16x64_final_sorting_scan_right_to_left_32s(T, SA, induction_bucket, 0, n);
}

static void libsais16x64_clear_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT bucket_start, sa_sint_t * RESTRICT bucket_end)
{
    fast_sint_t c;

    for (c = 0; c < k; ++c)
    {
        if (bucket_end[c] > bucket_start[c])
        {
            memset(&SA[bucket_start[c]], 0, ((size_t)bucket_end[c] - (size_t)bucket_start[c]) * sizeof(sa_sint_t));
        }
    }
}

static sa_sint_t libsais16x64_induce_final_order_16u_omp(const uint16_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t bwt, sa_sint_t r, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT buckets)
{
    if (!bwt)
    {
        libsais16x64_final_sorting_scan_left_to_right_16u_omp(T, SA, n, k, &buckets[6 * ALPHABET_SIZE]);
        libsais16x64_final_sorting_scan_right_to_left_16u_omp(T, SA, n, k, &buckets[7 * ALPHABET_SIZE]);
        return 0;
    }
    else if (I != NULL)
    {
        libsais16x64_final_bwt_aux_scan_left_to_right_16u_omp(T, SA, n, k, r - 1, I, &buckets[6 * ALPHABET_SIZE]);
        libsais16x64_final_bwt_aux_scan_right_to_left_16u_omp(T, SA, n, k, r - 1, I, &buckets[7 * ALPHABET_SIZE]);
        return 0;
    }
    else
    {
        libsais16x64_final_bwt_scan_left_to_right_16u_omp(T, SA, n, k, &buckets[6 * ALPHABET_SIZE]);
        return libsais16x64_final_bwt_scan_right_to_left_16u_omp(T, SA, n, k, &buckets[7 * ALPHABET_SIZE]);
    }
}

static void libsais16x64_induce_final_order_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(T, SA, n, &buckets[4 * (fast_sint_t)k]);
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(T, SA, n, &buckets[5 * (fast_sint_t)k]);
}

static void libsais16x64_induce_final_order_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(T, SA, n, &buckets[2 * (fast_sint_t)k]);
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(T, SA, n, &buckets[3 * (fast_sint_t)k]);
}

static void libsais16x64_induce_final_order_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_start_32s_1k(k, buckets);
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(T, SA, n, buckets);

    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_end_32s_1k(k, buckets);
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(T, SA, n, buckets);
}

static sa_sint_t libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t f, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT SAm = &SA[m];

    sa_sint_t i, j;
    for (i = (sa_sint_t)omp_block_start, j = (sa_sint_t)omp_block_start + (sa_sint_t)omp_block_size - 2 * (sa_sint_t)prefetch_distance - 3; i < j; i += 4)
    {
        libsais16x64_prefetchr(&SA[i + 3 * prefetch_distance]);

        libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 0]) >> 1]);
        libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 1]) >> 1]);
        libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 2]) >> 1]);
        libsais16x64_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 3]) >> 1]);

        sa_uint_t q0 = (sa_uint_t)SA[i + prefetch_distance + 0]; const sa_sint_t * Tq0 = &T[q0]; libsais16x64_prefetchw(SAm[q0 >> 1] < 0 ? Tq0 : NULL);
        sa_uint_t q1 = (sa_uint_t)SA[i + prefetch_distance + 1]; const sa_sint_t * Tq1 = &T[q1]; libsais16x64_prefetchw(SAm[q1 >> 1] < 0 ? Tq1 : NULL);
        sa_uint_t q2 = (sa_uint_t)SA[i + prefetch_distance + 2]; const sa_sint_t * Tq2 = &T[q2]; libsais16x64_prefetchw(SAm[q2 >> 1] < 0 ? Tq2 : NULL);
        sa_uint_t q3 = (sa_uint_t)SA[i + prefetch_distance + 3]; const sa_sint_t * Tq3 = &T[q3]; libsais16x64_prefetchw(SAm[q3 >> 1] < 0 ? Tq3 : NULL);

        sa_uint_t p0 = (sa_uint_t)SA[i + 0]; sa_sint_t s0 = SAm[p0 >> 1]; if (s0 < 0) { T[p0] |= SAINT_MIN; f++; s0 = i + 0 + SAINT_MIN + f; } SAm[p0 >> 1] = s0 - f;
        sa_uint_t p1 = (sa_uint_t)SA[i + 1]; sa_sint_t s1 = SAm[p1 >> 1]; if (s1 < 0) { T[p1] |= SAINT_MIN; f++; s1 = i + 1 + SAINT_MIN + f; } SAm[p1 >> 1] = s1 - f;
        sa_uint_t p2 = (sa_uint_t)SA[i + 2]; sa_sint_t s2 = SAm[p2 >> 1]; if (s2 < 0) { T[p2] |= SAINT_MIN; f++; s2 = i + 2 + SAINT_MIN + f; } SAm[p2 >> 1] = s2 - f;
        sa_uint_t p3 = (sa_uint_t)SA[i + 3]; sa_sint_t s3 = SAm[p3 >> 1]; if (s3 < 0) { T[p3] |= SAINT_MIN; f++; s3 = i + 3 + SAINT_MIN + f; } SAm[p3 >> 1] = s3 - f;
    }

    for (j += 2 * (sa_sint_t)prefetch_distance + 3; i < j; i += 1)
    {
        sa_uint_t p = (sa_uint_t)SA[i]; sa_sint_t s = SAm[p >> 1]; if (s < 0) { T[p] |= SAINT_MIN; f++; s = i + SAINT_MIN + f; } SAm[p >> 1] = s - f;
    }

    return f;
}

static void libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t * pl, fast_sint_t * pr, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    sa_sint_t * RESTRICT SAl = &SA[0];
    sa_sint_t * RESTRICT SAr = &SA[0];

    fast_sint_t i, j, l = *pl - 1, r = *pr - 1;
    for (i = (fast_sint_t)m + omp_block_start + omp_block_size - 1, j = (fast_sint_t)m + omp_block_start + 3; i >= j; i -= 4)
    {
        libsais16x64_prefetchr(&SA[i - prefetch_distance]);

        sa_sint_t p0 = SA[i - 0]; SAl[l] = p0 & SAINT_MAX; l -= p0 < 0; SAr[r] = p0 - 1; r -= p0 > 0;
        sa_sint_t p1 = SA[i - 1]; SAl[l] = p1 & SAINT_MAX; l -= p1 < 0; SAr[r] = p1 - 1; r -= p1 > 0;
        sa_sint_t p2 = SA[i - 2]; SAl[l] = p2 & SAINT_MAX; l -= p2 < 0; SAr[r] = p2 - 1; r -= p2 > 0;
        sa_sint_t p3 = SA[i - 3]; SAl[l] = p3 & SAINT_MAX; l -= p3 < 0; SAr[r] = p3 - 1; r -= p3 > 0;
    }

    for (j -= 3; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SAl[l] = p & SAINT_MAX; l -= p < 0; SAr[r] = p - 1; r -= p > 0;
    }
    
    *pl = l + 1; *pr = r + 1;
}

static sa_sint_t libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t m)
{
    sa_sint_t f = 0;

    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = m - omp_block_start;

    f = libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s(T, SA, m, 0, omp_block_start, omp_block_size);

    return f;
}

static void libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t f)
{

    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = ((fast_sint_t)n >> 1) - omp_block_start;

    fast_sint_t l = m, r = (fast_sint_t)n + (fast_sint_t)fs;
    libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s(SA, m, &l, &r, omp_block_start, omp_block_size);

    memcpy(&SA[(fast_sint_t)n + (fast_sint_t)fs - (fast_sint_t)m], &SA[(fast_sint_t)m - (fast_sint_t)f], (size_t)f * sizeof(sa_sint_t));
}

static sa_sint_t libsais16x64_compact_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs)
{
    sa_sint_t f = libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s_omp(T, SA, m);
    libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s_omp(SA, n, m, fs, f);

    return f;
}

static void libsais16x64_merge_unique_lms_suffixes_32s(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, fast_sint_t l, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    const sa_sint_t * RESTRICT SAnm = &SA[(fast_sint_t)n - (fast_sint_t)m - 1 + l];

    sa_sint_t i, j; fast_sint_t tmp = *SAnm++;
    for (i = (sa_sint_t)omp_block_start, j = (sa_sint_t)omp_block_start + (sa_sint_t)omp_block_size - 6; i < j; i += 4)
    {
        libsais16x64_prefetchr(&T[i + prefetch_distance]);

        sa_sint_t c0 = T[i + 0]; if (c0 < 0) { T[i + 0] = c0 & SAINT_MAX; SA[tmp] = i + 0; i++; tmp = *SAnm++; }
        sa_sint_t c1 = T[i + 1]; if (c1 < 0) { T[i + 1] = c1 & SAINT_MAX; SA[tmp] = i + 1; i++; tmp = *SAnm++; }
        sa_sint_t c2 = T[i + 2]; if (c2 < 0) { T[i + 2] = c2 & SAINT_MAX; SA[tmp] = i + 2; i++; tmp = *SAnm++; }
        sa_sint_t c3 = T[i + 3]; if (c3 < 0) { T[i + 3] = c3 & SAINT_MAX; SA[tmp] = i + 3; i++; tmp = *SAnm++; }
    }

    for (j += 6; i < j; i += 1)
    {
        sa_sint_t c = T[i]; if (c < 0) { T[i] = c & SAINT_MAX; SA[tmp] = i; i++; tmp = *SAnm++; }
    }
}

static void libsais16x64_merge_nonunique_lms_suffixes_32s(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, fast_sint_t l, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 32;

    const sa_sint_t * RESTRICT SAnm = &SA[(fast_sint_t)n - (fast_sint_t)m - 1 + l];

    fast_sint_t i, j; sa_sint_t tmp = *SAnm++;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3; i < j; i += 4)
    {
        libsais16x64_prefetchr(&SA[i + prefetch_distance]);

        if (SA[i + 0] == 0) { SA[i + 0] = tmp; tmp = *SAnm++; }
        if (SA[i + 1] == 0) { SA[i + 1] = tmp; tmp = *SAnm++; }
        if (SA[i + 2] == 0) { SA[i + 2] = tmp; tmp = *SAnm++; }
        if (SA[i + 3] == 0) { SA[i + 3] = tmp; tmp = *SAnm++; }
    }

    for (j += 3; i < j; i += 1)
    {
        if (SA[i] == 0) { SA[i] = tmp; tmp = *SAnm++; }
    }
}

static void libsais16x64_merge_unique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = n - omp_block_start;
    
    libsais16x64_merge_unique_lms_suffixes_32s(T, SA, n, m, 0, omp_block_start, omp_block_size);
}

static void libsais16x64_merge_nonunique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t f)
{
    fast_sint_t omp_block_start   = 0;
    fast_sint_t omp_block_size    = m - omp_block_start;

    libsais16x64_merge_nonunique_lms_suffixes_32s(SA, n, m, f, omp_block_start, omp_block_size);
}

static void libsais16x64_merge_compacted_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t f)
{
    libsais16x64_merge_unique_lms_suffixes_32s_omp(T, SA, n, m);
    libsais16x64_merge_nonunique_lms_suffixes_32s_omp(SA, n, m, f);
}

static void libsais16x64_reconstruct_compacted_lms_suffixes_32s_2k_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, sa_sint_t fs, sa_sint_t f, sa_sint_t * RESTRICT buckets)
{
    if (f > 0)
    {
        memmove(&SA[n - m - 1], &SA[n + fs - m], (size_t)f * sizeof(sa_sint_t));

        libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets);
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m - f);

        memcpy(&SA[n - m - 1 + f], &SA[0], ((size_t)m - (size_t)f) * sizeof(sa_sint_t));
        memset(&SA[0], 0, (size_t)m * sizeof(sa_sint_t));

        libsais16x64_merge_compacted_lms_suffixes_32s_omp(T, SA, n, m, f);
    }
    else
    {
        libsais16x64_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets, 0, n);
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m);
    }
}

static void libsais16x64_reconstruct_compacted_lms_suffixes_32s_1k_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t f)
{
    if (f > 0)
    {
        memmove(&SA[n - m - 1], &SA[n + fs - m], (size_t)f * sizeof(sa_sint_t));

        libsais16x64_gather_compacted_lms_suffixes_32s(T, SA, n);
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m - f);

        memcpy(&SA[n - m - 1 + f], &SA[0], ((size_t)m - (size_t)f) * sizeof(sa_sint_t));
        memset(&SA[0], 0, (size_t)m * sizeof(sa_sint_t));

        libsais16x64_merge_compacted_lms_suffixes_32s_omp(T, SA, n, m, f);
    }
    else
    {
        libsais16x64_gather_lms_suffixes_32s(T, SA, n);
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m);
    }
}

static sa_sint_t libsais16x64_main_32s_recursion(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t fs, sa_sint_t * RESTRICT local_buffer)
{
    fs = fs < (SAINT_MAX - n) ? fs : (SAINT_MAX - n);

    if (k > 0 && ((fs / k >= 6) || (LIBSAIS_LOCAL_BUFFER_SIZE / k >= 6)))
    {
        sa_sint_t alignment = (fs - 1024) / k >= 6 ? (sa_sint_t)1024 : (sa_sint_t)16;
        sa_sint_t * RESTRICT buckets = (fs - alignment) / k >= 6 ? (sa_sint_t *)libsais16x64_align_up(&SA[n + fs - 6 * (fast_sint_t)k - alignment], (size_t)alignment * sizeof(sa_sint_t)) : &SA[n + fs - 6 * (fast_sint_t)k];
        buckets = (LIBSAIS_LOCAL_BUFFER_SIZE / k >= 6) ? local_buffer : buckets;

        sa_sint_t m = libsais16x64_count_and_gather_lms_suffixes_32s_4k_omp(T, SA, n, k, buckets);
        if (m > 1)
        {
            memset(SA, 0, ((size_t)n - (size_t)m) * sizeof(sa_sint_t));

            sa_sint_t first_lms_suffix    = SA[n - m];
            sa_sint_t left_suffixes_count = libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(T, k, buckets, first_lms_suffix);

            libsais16x64_radix_sort_lms_suffixes_32s_6k_omp(T, SA, n, m, &buckets[4 * (fast_sint_t)k]);

            if ((n / 8192) < k) { libsais16x64_radix_sort_set_markers_32s_6k_omp(SA, k, &buckets[4 * (fast_sint_t)k]); }

            libsais16x64_initialize_buckets_for_partial_sorting_32s_6k(T, k, buckets, first_lms_suffix, left_suffixes_count);
            libsais16x64_induce_partial_order_32s_6k_omp(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count);

            sa_sint_t names = (n / 8192) < k
                ? libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(SA, n, m)
                : libsais16x64_renumber_and_gather_lms_suffixes_omp(SA, n, m, fs);

            if (names < m)
            {
                sa_sint_t f = (n / 8192) < k
                    ? libsais16x64_compact_lms_suffixes_32s_omp(T, SA, n, m, fs)
                    : 0;

                if (libsais16x64_main_32s_recursion(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f, local_buffer) != 0)
                {
                    return -2;
                }

                libsais16x64_reconstruct_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, m, fs, f, buckets);
            }
            else
            {
                libsais16x64_count_lms_suffixes_32s_2k(T, n, k, buckets);
            }

            libsais16x64_initialize_buckets_start_and_end_32s_4k(k, buckets);
            libsais16x64_place_lms_suffixes_histogram_32s_4k(SA, n, k, m, buckets);
            libsais16x64_induce_final_order_32s_4k(T, SA, n, k, buckets);
        }
        else
        {
            SA[0] = SA[n - 1];

            libsais16x64_initialize_buckets_start_and_end_32s_6k(k, buckets);
            libsais16x64_place_lms_suffixes_histogram_32s_6k(SA, n, k, m, buckets);
            libsais16x64_induce_final_order_32s_6k(T, SA, n, k, buckets);
        }

        return 0;
    } 
    else 
    {
        sa_sint_t * buffer = fs < k ? (sa_sint_t *)libsais16x64_alloc_aligned((size_t)k * sizeof(sa_sint_t), 4096) : (sa_sint_t *)NULL;

        sa_sint_t alignment = fs - 1024 >= k ? (sa_sint_t)1024 : (sa_sint_t)16;
        sa_sint_t * RESTRICT buckets = fs - alignment >= k ? (sa_sint_t *)libsais16x64_align_up(&SA[n + fs - k - alignment], (size_t)alignment * sizeof(sa_sint_t)) : fs >= k ? &SA[n + fs - k] : buffer;

        if (buckets == NULL) { return -2; }

        memset(SA, 0, (size_t)n * sizeof(sa_sint_t));

        libsais16x64_count_suffixes_32s(T, n, k, buckets); 
        libsais16x64_initialize_buckets_end_32s_1k(k, buckets);

        sa_sint_t m = libsais16x64_radix_sort_lms_suffixes_32s_1k(T, SA, n, buckets);
        if (m > 1)
        {
            libsais16x64_induce_partial_order_32s_1k_omp(T, SA, n, k, buckets);

            sa_sint_t names = libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(T, SA, n, m);
            if (names < m)
            {
                if (buffer != NULL) { libsais16x64_free_aligned(buffer); buckets = NULL; }

                sa_sint_t f = libsais16x64_compact_lms_suffixes_32s_omp(T, SA, n, m, fs);

                if (libsais16x64_main_32s_recursion(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f, local_buffer) != 0)
                {
                    return -2;
                }

                libsais16x64_reconstruct_compacted_lms_suffixes_32s_1k_omp(T, SA, n, m, fs, f);

                if (buckets == NULL) { buckets = buffer = (sa_sint_t *)libsais16x64_alloc_aligned((size_t)k * sizeof(sa_sint_t), 4096); }
                if (buckets == NULL) { return -2; }
            }
            
            libsais16x64_count_suffixes_32s(T, n, k, buckets);
            libsais16x64_initialize_buckets_end_32s_1k(k, buckets);
            libsais16x64_place_lms_suffixes_interval_32s_1k(T, SA, k, m, buckets);
        }

        libsais16x64_induce_final_order_32s_1k(T, SA, n, k, buckets);
        libsais16x64_free_aligned(buffer);

        return 0;
    }
}

static sa_sint_t libsais16x64_main_32s_entry(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t fs)
{
    sa_sint_t local_buffer[LIBSAIS_LOCAL_BUFFER_SIZE];

    return libsais16x64_main_32s_recursion(T, SA, n, k, fs, local_buffer);
}

static sa_sint_t libsais16x64_main_16u(const uint16_t * T, sa_sint_t * SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t bwt, sa_sint_t r, sa_sint_t * RESTRICT I, sa_sint_t fs, sa_sint_t * freq)
{
    fs = fs < (SAINT_MAX - n) ? fs : (SAINT_MAX - n);

    sa_sint_t m = libsais16x64_count_and_gather_lms_suffixes_16u_omp(T, SA, n, buckets);
    sa_sint_t k = libsais16x64_initialize_buckets_start_and_end_16u(buckets, freq);

    if (m > 0)
    {
        sa_sint_t first_lms_suffix    = SA[n - m];
        sa_sint_t left_suffixes_count = libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_16u(T, buckets, first_lms_suffix);

        libsais16x64_radix_sort_lms_suffixes_16u_omp(T, SA, n, m, buckets);

        libsais16x64_initialize_buckets_for_partial_sorting_16u(T, buckets, first_lms_suffix, left_suffixes_count);
        libsais16x64_induce_partial_order_16u_omp(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count);

        sa_sint_t names = libsais16x64_renumber_and_gather_lms_suffixes_omp(SA, n, m, fs);
        if (names < m)
        {
            if (libsais16x64_main_32s_entry(SA + n + fs - m, SA, m, names, fs + n - 2 * m) != 0)
            {
                return -2;
            }

            libsais16x64_gather_lms_suffixes_16u_omp(T, SA, n);
            libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m);
        }

        libsais16x64_place_lms_suffixes_interval_16u(SA, n, m, buckets);
    }
    else
    {
        memset(SA, 0, (size_t)n * sizeof(sa_sint_t));
    }

    return libsais16x64_induce_final_order_16u_omp(T, SA, n, k, bwt, r, I, buckets);
}

static sa_sint_t libsais16x64_main(const uint16_t * T, sa_sint_t * SA, sa_sint_t n, sa_sint_t bwt, sa_sint_t r, sa_sint_t * I, sa_sint_t fs, sa_sint_t * freq)
{
    sa_sint_t *             RESTRICT buckets        = (sa_sint_t *)libsais16x64_alloc_aligned((size_t)8 * ALPHABET_SIZE * sizeof(sa_sint_t), 4096);

    sa_sint_t index = buckets != NULL
        ? libsais16x64_main_16u(T, SA, n, buckets, bwt, r, I, fs, freq)
        : -2;

    libsais16x64_free_aligned(buckets);

    return index;
}

int64_t libsais16x64(const uint16_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (fs < 0))
    {
        return -1;
    }
    else if (n < 2)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int64_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    return libsais16x64_main(T, SA, n, 0, 0, NULL, fs, freq);
}
