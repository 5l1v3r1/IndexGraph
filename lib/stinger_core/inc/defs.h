/*
 * Copyright 2018-now Yilong Liu
 * @Author: Yilong Liu
 * @Last Modified by: Yilong Liu
 */


#ifndef SRC_COMMON_DEFS_H_
#define SRC_COMMON_DEFS_H_

#ifndef PAGE_SHIFT
#define PAGE_SHIFT 12
#endif

#ifndef PAGE_SIZE
#ifdef __ASSEMBLY__
#define PAGE_SIZE       (1 << PAGE_SHIFT)
#else
#ifdef __x86_64__
#define PAGE_SIZE       (1UL << PAGE_SHIFT)
#else
#define PAGE_SIZE       (1ULL << PAGE_SHIFT)
#endif
#endif
#endif

#ifndef SMALL_PAGE
#define SMALL_PAGE PAGE_SIZE
#endif

#ifndef HUGE_PAGE
#define HUGE_PAGE (PAGE_SIZE << 9)
#endif

#ifndef STACK_SIZE_PAGE_ORDER
#define STACK_SIZE_PAGE_ORDER  4
#endif

#ifndef STACK_SIZE
#define STACK_SIZE             (PAGE_SIZE * (1 << STACK_SIZE_PAGE_ORDER))
#endif

#define EDGES_PER_ROOM 8
#define EMPTY 0
#define FULL  1

#define DEFAULT_CACHE_SIZE (1024ULL * 1024ULL)

#define THRESHOLD EDGES_PER_ROOM

#define STINGER_EDGE_DIRECTION_MASK (0x6000000000000000L)
#define STINGER_EDGE_DIRECTION_OUT (0x4000000000000000L)
#define STINGER_EDGE_DIRECTION_IN (0x2000000000000000L)

#endif  // SRC_COMMON_DEFS_H_
