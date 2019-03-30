/*
 * Copyright 2018-now Zhihui Zhao
 * @Author: Zhihui Zhao
 */

#include <cstdint>
#include <cstdlib>
#include <cstdio>

// @example: read write lock in boost library
// #include <boost/thread/locks.hpp>
// #include <boost/thread/shared_mutex.hpp>

// typedef boost::shared_mutex Lock;
// typedef boost::unique_lock< Lock > WriteLock;
// typedef boost::shared_lock< Lock > ReadLock;

// static Lock myLock;

// void ReadFunction()
// {
//     ReadLock r_lock(myLock);
//     //Do reader stuff
// }

// void WriteFunction()
// {
//      WriteLock w_lock(myLock);
//      //Do writer stuff
// }

#include "stinger_core/defs.h"
#include "stinger_core/stinger_atomics.h"
#include "stinger_core/stinger_error.h"

#include "stinger_core/bitmap.hpp"

#include "stinger_core/vertex_array.h"
#include "stinger_core/shared_buffer.h"

// SharedBuffer 相关的方法，但是这里只考虑了源节点的出边情况

// 申请空间，初始化 SharedBuffer 这个管理结构
SharedBuffer::SharedBuffer() {
    LOG_D("Entering Shared_Buffer Constructor");

    // 或者后面根据顶点数来确定大小
    this->buffer_size = 28ULL * 1024 * 1024 * 1024;

    this->buffer = new Edge[this->buffer_size/sizeof(Edge)];
    for (int64_t i = 0; i < this->buffer_size/sizeof(Edge); ++i) {
        this->buffer[i].src = -1;
        this->buffer[i].dst = -1;
    }

    this->room_size = EDGES_PER_ROOM;     // 每个room能存多少条边
    this->room_num = this->buffer_size/(this->room_size*sizeof(Edge));
    this->room_allocated_limit = 0;

    this->bit_buffer = new Bitmap(size_t(this->room_num));   // 全置为0，表示有空间

    this->room_lock = new std::mutex[room_num];
    //this->room_lock[room_num];

    //this->va = va;
}

// 找到第一个可用的room,返回room_id号
int64_t SharedBuffer::allocate_room(const uint64_t src) {
    // when data size was 100w, 500000 rooms were allocated
    if (room_allocated_limit == room_num-1)
        room_allocated_limit = 0;
    for (int64_t i = room_allocated_limit; i < room_num; ++i) {
        if (bit_buffer->get_bit(i) == EMPTY) {
            std::unique_lock<std::mutex> lock(room_lock[i]);
            if (bit_buffer->get_bit(i) == EMPTY) {
                room_allocated_limit = i;
                bit_buffer->set_bit(i);
                return i;
            } else {
                // occupy by anther thread
                lock.unlock();
                continue;
            }
        }
    }

    LOG_I("No free shared_buffer room");
    return -1;
}

// 根据找到的room_id号，找到房间的空位置来放边
int SharedBuffer::push_room(int64_t id, const Edge &e, VertexArray * va, int64_t direction) {
    int ret = 0;
    if (id == -1) printf("wrong push into shared\n");
 //   std::unique_lock<std::mutex> lock(room_lock[id]);

    for (int64_t i = id * room_size;
        i < ((id + 1) * room_size); ++i) {
        if (buffer[i].src == e.src && buffer[i].dst == e.dst) {
        //    printf("this edge is already in this room\n");
            break;
        }
        if (buffer[i].src == -1 && buffer[i].dst == -1) {
        //    std::unique_lock<std::mutex> lock(room_lock[id]);
            buffer[i].src = e.src;     // 找到房间的一个空隙
            buffer[i].dst = e.dst;
            ret = 1;
            break;
        }
    }

    // 按理说前面会判断 degree大小，所以应该不会出现这种情况
//    printf("there is no empty slot in this room\n");
    return ret;
}

void SharedBuffer::free_room(int64_t id) {
    int64_t room_id = id;    // 得到src所在的room

    if (bit_buffer->get_bit(room_id) != EMPTY) {
        if (room_allocated_limit > room_id) room_allocated_limit = room_id;
        bit_buffer->clear_bit(room_id);
    }

    for (int64_t i = room_id * room_size;
        i < int64_t((room_id + 1) * room_size); ++i) {
        set_invalid(&buffer[i]);   // 将里面的边置为无效
    }
}

void SharedBuffer::move_to_index(Edge * e, int64_t id) {
//    printf("Starting to move all Shared_Buffer edges in room %d\n", id);
 //   std::unique_lock<std::mutex> lock(room_lock[id]);
    int cnt = 0;
    for (int64_t i = id * room_size;
        i < ((id + 1) * room_size); ++i) {
        e[cnt].src = buffer[i].src;
        e[cnt].dst = buffer[i].dst;
        ++cnt;
        //__sync_add_and_fetch(&cnt, 1);

        set_invalid(&buffer[i]);
    }

    // not called free fun(), so push it here
    if (bit_buffer->get_bit(id) != EMPTY) {
        if (room_allocated_limit > id) room_allocated_limit = id;
        bit_buffer->clear_bit(id);
    }
    //bit_buffer->clear_bit(id);

  //  lock.unlock();
}

void SharedBuffer::set_invalid(Edge *e) {    // 将边e设置为无效
    e->src = -1;
    e->dst = -1;
    e->weight = 0;
    e->time_Recent = 0;
}

SharedBuffer::~SharedBuffer() {
    LOG_D("Entering Shared_Buffer Descontructor");
  //  delete []buffer;
    delete this->bit_buffer;
    this->buffer = nullptr;
    delete this->room_lock;
    this->room_lock = nullptr;
    delete this->buffer;
    this->buffer = nullptr;
}

bool SharedBuffer::find_in_room(int64_t id, int64_t dest) {
    for (int64_t i = id * room_size;
        i < ((id + 1) * room_size); ++i) {
        if (buffer[i].dst == dest) {
//            printf("this edge is already in this room\n");
            return true;
        }
    }
    return false;
}

void SharedBuffer::show() {
    printf("show the context of shared buffer\n");
    for (int i = 0; i < room_num; ++i) {
        if (bit_buffer->get_bit(i) == EMPTY)
            continue;
        printf("room %d:\n",i);
        for (int j = i*room_size; j < (i+1)*room_size; ++j) {
            if (buffer[j].src == -1 || buffer[j].dst == -1)
                break;
            printf("<%ld,%ld>,  ", buffer[j].src, buffer[j].dst);
        }
        printf("\n");
    }
}
