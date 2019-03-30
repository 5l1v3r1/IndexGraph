/*
 * Copyright 2018-now Zhihui Zhao
 * @Author: Zhihui Zhao
 * @Last Modified by: Yilong Liu
 */

#ifndef SRC_CORE_SHARED_BUFFER_H_
#define SRC_CORE_SHARED_BUFFER_H_

#include <cstdio>
#include <cstdint>
#include <mutex>

#include "stinger_core/bitmap.hpp"

#include "stinger_core/vertex_array.h"

// 冷点的存储区
class SharedBuffer{
public:
    int64_t room_num;       // 有几个房间
    int64_t room_allocated_limit; // lowest bound index of free room
//    uint32_t room_used;      // occupied room nums
//    uint8_t room_size;      // 一个room的size大小，如假设能存16条边

    int64_t buffer_size;      // 整个buffer区的大小，如4G
    Bitmap* bit_buffer;   // 记录各个房间是否被占用

    std::mutex buffer_lock;
    std::mutex * room_lock;   // lock a room

    //VertexArray *va;

    Edge* buffer;      // buffer区
    int64_t room_size;      // 一个room的size大小，如假设能存16条边

    SharedBuffer();
    ~SharedBuffer();

    //int allocate_room(const Edge &e);
    int64_t allocate_room(const uint64_t src);
    int push_room(int64_t id, const Edge &e,  VertexArray * va, int64_t direction);

    void free_room(int64_t id);
    void set_invalid(Edge *e);

    bool find_in_room(int64_t id, int64_t dest);
    void move_to_index(Edge *e, int64_t id);

    Edge* get_buffer() {return buffer;}

    // printf shared context
    void show();

};

#endif  // SRC_CORE_SHARED_BUFFER_H_

