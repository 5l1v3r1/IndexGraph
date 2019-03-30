/*
 * Copyright 2018-now Zhihui Zhao
 * @Author: Zhihui Zhao
 * @Last Modified by: Yilong Liu
 */

#ifndef SRC_CORE_VERTEX_ARRAY_H_
#define SRC_CORE_VERTEX_ARRAY_H_

#include <cstdint>
#include <mutex>
#include <atomic>
#include <boost/thread/locks.hpp>
#include <boost/thread/shared_mutex.hpp>
#include "stinger_core/stinger_error.h"
#include "stinger_core/stinger_atomics.h"
#include "stinger_core/indirect_index.h"

typedef struct Vertex Vertex;
typedef struct Edge Edge;
typedef struct VertexArray VertexArray;

// A vertex
struct Vertex {
    uint64_t in_degree;
    uint64_t out_degree;
    int64_t out_buffer_index;  // 出边存在shared_buffer的第几个room
    int64_t in_buffer_index;

    // @ADD: in_data, out_data to pointer indirect_index(may exist)
    void *in_data;
    void *out_data;
};

// A single edge in Shared Buffer and Indirect_index
struct Edge {        // stinger 是4个8字节
    int64_t src;
    int64_t dst;
    int64_t weight;
    // uint64_t time_First;
    int64_t time_Recent;
};

struct VertexArray {
    Vertex* vertices;
    int64_t num_vertices;
    IndirectIndex * index;
};

// VertexArray functions
VertexArray *new_vertex_array(int64_t num_vertices);
void free_vertex_array(VertexArray *vertex_array);


inline uint64_t get_outdegree(const VertexArray *va, const int64_t src) { 
    return va->vertices[src].out_degree;
}

inline uint64_t increment_outdegree(const VertexArray *va, const int64_t src) {
    return __sync_add_and_fetch(&va->vertices[src].out_degree, 1);
}

inline uint64_t get_indegree(const VertexArray *va, const int64_t src) { 
    return va->vertices[src].in_degree;
}

inline uint64_t increment_indegree(const VertexArray *va, const int64_t src) {
    return __sync_add_and_fetch(&va->vertices[src].in_degree, 1);
}

inline uint64_t get_bufferindex(volatile uint64_t * i) {
    uint64_t index = *i;
    return index;
}

inline void set_inbuffer(const VertexArray *va, const int64_t src, const int id) {
    if (va->vertices[src].in_buffer_index == -1)
        va->vertices[src].in_buffer_index = id;
}

inline void set_outbuffer(const VertexArray *va, const int64_t src, const int id) {
    __sync_val_compare_and_swap(&va->vertices[src].out_buffer_index, -1, id);
}

inline IndirectIndex* get_indata(const VertexArray * va, const int64_t src) {
    return (IndirectIndex*)va->vertices[src].in_data;
}
inline IndirectIndex* get_outdata(const VertexArray * va, const int64_t src) {
    return (IndirectIndex*)va->vertices[src].out_data;
}

inline void set_indata(const VertexArray *va, const int64_t src,  IndirectIndex * rtn) {
    __sync_val_compare_and_swap(&va->vertices[src].in_data, nullptr, rtn);
}
inline void set_outdata(const VertexArray *va, const int64_t src,  IndirectIndex * rtn) {
    __sync_val_compare_and_swap(&va->vertices[src].out_data, nullptr, rtn);
}

#endif  // SRC_CORE_VERTEX_ARRAY_H_
