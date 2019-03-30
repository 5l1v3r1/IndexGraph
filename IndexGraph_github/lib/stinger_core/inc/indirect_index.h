/*
 * Copyright 2018-now Zhihui Zhao
 * @Author: Zhihui Zhao
 * @Last Modified by: Yilong Liu
 */

#ifndef SRC_CORE_INDIRECT_INDEX_H_
#define SRC_CORE_INDIRECT_INDEX_H_

#include <cstdint>
#include <mutex>
#include <thread>

#include "stinger_core/bitmap.hpp"
//#include "stinger_core/vertex_array.h"
#include "stinger_core/stinger_error.h"
#include "stinger_core/defs.h"

#define INDEX_SIZE ((PAGE_SIZE-16)/sizeof(uint64_t*))
typedef struct Edge Edge;

class IndirectIndex {
public:
    //int16_t index_size;
    uint32_t index_used;

    uint32_t arr_size;

    // add for beap
    uint32_t heap_height;
    uint32_t heap_size;    // total 16 Byte above

    Edge *first_index[INDEX_SIZE];  // 一级索引
    IndirectIndex * next;

    IndirectIndex();
    IndirectIndex(int size);
    ~IndirectIndex();

    void set_array(Edge *arr);

  //  void update_indirect(int64_t const &, int index_id, int arr_id);
    void push_indirect(int64_t const &);
    void push_nobeep(int64_t const &);
    bool find_in_indirect(int64_t const &, int &index_id, int &arr_id);
    bool find_nobeep(int64_t const &, int &index_id, int &arr_id);

    void move_into_indirect(Edge e[], int num);
    //Edge* get_first_index()  { return first_index; }

    void show();
    void do_show();

    void push_append(int64_t const &obj);
    void do_push_append(int64_t const &obj);
    void make_beap();
    void quick_sort(int l, int r);

    bool find(int64_t const &);
    void push(int64_t const &);
    bool do_find(int64_t const &);
    void do_push(int64_t const &);
    void left_push(int64_t const&, int, int);
    void right_push(int64_t const&, int, int);
    static int first(int);
    static int last(int);  
};

#endif  // SRC_CORE_INDIRECT_INDEX_H_

