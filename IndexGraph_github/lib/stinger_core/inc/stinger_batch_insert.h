/*
 * stinger_batch_insert.h
 * Author: Eric Hein <ehein6@gatech.edu>
 * Date: 10/24/2016
 * Purpose:
 *   Provides optimized edge insert/update routines for stinger.
 *   Multiple updates for the same source vertex are dispatched to each thread,
 *   reducing the number of edge-list traversals that need to be done.
 *
 * Accepts a pair of random-access iterators to the range of edge updates to perform.
 * Caller must also provide an adapter template argument: a struct of static functions
 * to access the source, destination, weight, time, and result code of the update.
 * This implementation handles duplicates and sets the return code correctly.
 */

#ifndef STINGER_BATCH_INSERT_H_
#define STINGER_BATCH_INSERT_H_

#include "stinger.h"
#include "stinger_internal.h"
#include "stinger_atomics.h"
#include "x86_full_empty.h"
#define LOG_AT_I
#include "stinger_error.h"
#undef LOG_AT_I

#include "stinger_core/defs.h"
#include "stinger_core/indirect_index.h"
#include "stinger_core/shared_buffer.h"
#include "stinger_core/vertex_array.h"
#include "stinger_net/stinger_server_state.h"
#include "stinger_utils/timer.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <thread>
#include <sys/time.h>
#include <unordered_map>
#include <list>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include<unistd.h>
#define NUM_THREADS 28

using namespace gt::stinger;
using namespace std;
using std::vector;
int core_num = 28;
int mode = 2;  // default mode:hybrid
extern int vertex_num;
extern std::atomic_long rcon_cnt;
extern std::atomic_long rcon_time;
extern double mode_time;
extern double binary_time, beap_time;
extern int who;

// *** Public interface (definitions at end of file) ***
template<typename adapter, typename iterator>
void stinger_batch_incr_edges(SharedBuffer * shared, VertexArray * va, iterator begin, iterator end);
template<typename adapter, typename iterator>
void stinger_batch_insert_edges(stinger_t * G, iterator begin, iterator end);

// re-write
template<typename adapter, typename iterator>
void stinger_batch_incr_edge_pairs(SharedBuffer * shared, VertexArray * va, iterator begin, iterator end);

template<typename adapter, typename iterator>
void stinger_batch_insert_edge_pairs(stinger_t * G, iterator begin, iterator end);

// *** Implementation ***
namespace gt { namespace stinger {

/*
 * Rather than add these template arguments to every function in the file, just wrap everything in a template class
 *
 * adapter - provides static methods for accessing fields of an update
 * iterator - iterator for a collection of updates
 */
template<typename adapter, typename iterator>
class BatchInserter
{
protected:
    // Everything in this class is protected and static, these friend functions are the only public interface
    BatchInserter() {}
    friend void stinger_batch_incr_edges<adapter, iterator>(SharedBuffer * shared,  VertexArray * va, iterator begin, iterator end);
    friend void stinger_batch_insert_edges<adapter, iterator>(stinger_t * G,iterator begin, iterator end);

    // @re-write
    friend void stinger_batch_incr_edge_pairs<adapter, iterator>(SharedBuffer * shared, VertexArray * va, iterator begin, iterator end);
    friend void stinger_batch_insert_edge_pairs<adapter, iterator>(stinger_t * G, iterator begin, iterator end);

    // what 'iterator' points to
    typedef typename std::iterator_traits<iterator>::value_type update;

    // Result codes
    // Caller initializes result code to 0, and expects 0 if edge is already present.
    // But we need to use this field to track if the update has been processed yet.
    // We use a fixed offset, then subtract it before returning to match behavior of functions like stinger_incr_edge.
    static const int64_t result_code_offset = 10;
    enum result_codes {
        PENDING =           0,
        EDGE_ADDED =        1 + result_code_offset,
        EDGE_UPDATED =      0 + result_code_offset,
        EDGE_NOT_ADDED =   -1 + result_code_offset
    };

    // Set of functions to sort/filter a collection of updates by source vertex
    struct source_funcs
    {
        static int64_t
        get(const update &x){
            return adapter::get_source(x);
        }
        static void
        set(update &x, int64_t neighbor){
            adapter::set_source(x, neighbor);
        }
        static bool
        equals(const update &a, const update &b){
            return adapter::get_source(a) == adapter::get_source(b);
        }
        static bool
        compare(const update &a, const update &b){
            return adapter::get_source(a) < adapter::get_source(b);
        }
        static bool
        sort(const update &a, const update &b){
            if (adapter::get_type(a) != adapter::get_type(b))
                return adapter::get_type(a) < adapter::get_type(b);
            if (adapter::get_source(a) != adapter::get_source(b))
                return adapter::get_source(a) < adapter::get_source(b);
            if (adapter::get_dest(a) != adapter::get_dest(b))
                return adapter::get_dest(a) < adapter::get_dest(b);
            if (adapter::get_time(a) != adapter::get_time(b))
                return adapter::get_time(a) < adapter::get_time(b);
            return false;
        }
    };

    // Set of functions to sort/filter a collection of updates by destination vertex
    struct dest_funcs
    {
        static int64_t
        get(const update &x){
            return adapter::get_dest(x);
        }
        static void
        set(update &x, int64_t neighbor){
            adapter::set_dest(x, neighbor);
        }
        static bool
        equals(const update &a, const update &b){
            return adapter::get_dest(a) == adapter::get_dest(b);
        }
        static bool
        compare(const update &a, const update &b){
            return adapter::get_dest(a) < adapter::get_dest(b);
        }
        static bool
        sort(const update &a, const update &b){
            if (adapter::get_type(a) != adapter::get_type(b))
                return adapter::get_type(a) < adapter::get_type(b);
            if (adapter::get_dest(a) != adapter::get_dest(b))
                return adapter::get_dest(a) < adapter::get_dest(b);
            if (adapter::get_source(a) != adapter::get_source(b))
                return adapter::get_source(a) < adapter::get_source(b);
            if (adapter::get_time(a) != adapter::get_time(b))
                return adapter::get_time(a) < adapter::get_time(b);
            return false;
        }
    };

    // Finds an element in a sorted range using binary search
    // http://stackoverflow.com/a/446327/1877086
    template<class Iter, class T, class Compare>
    static Iter
    binary_find(Iter begin, Iter end, T val, Compare comp)
    {
        // Finds the lower bound in at most log(last - first) + 1 comparisons
        Iter i = std::lower_bound(begin, end, val, comp);

        if (i != end && !(comp(val, *i)))
            return i; // found
        else
            return end; // not found
    }

    // Find the first pending update with destination 'dest'
    template<class use_dest>
    static iterator
    find_updates(iterator begin, iterator end, int64_t neighbor)
    {
        // Create a dummy update object with the neighbor we are looking for
        update key;
        use_dest::set(key, neighbor);
        // Find the updates for this neighbor
        iterator pos = binary_find(begin, end, key, use_dest::compare);
        // If the first update is not pending, we have already done all the updates for this neighbor
        if (pos != end && adapter::get_result(*pos) == PENDING) return pos;
        else return end;
    }

    // Keeps track of the next pending update to insert into the graph
    class next_update_tracker
    {
    protected:
        iterator pos;
        iterator end;
    public:
        next_update_tracker(iterator begin, iterator end)
        : pos(begin), end(end) {}

         iterator
        operator () ()
        {
            while (pos != end && adapter::get_result(*pos) != PENDING) { ++pos; }
            return pos;
        }    
    };

    class last_update_tracker {
        public:
        iterator begin;
        iterator pos;
        last_update_tracker(iterator begin, iterator end): pos(end), begin(begin) {}
        iterator operator () () {
            while (pos != begin-1 && adapter::get_result(*pos) != PENDING) {--pos;}
            return pos;
        }
    };

    /*
     * Arguments:
     * result - Result code to set on first update for a neighbor. Result of duplicate updates will always be EDGE_UPDATED.
     * create - Are we inserting into an empty slot, or just updating an existing slot?
     * pos - pointer to first update
     * updates_end - pointer to end of update list
     * G - pointer to STINGER
     * eb - pointer to edge block
     * e - edge index within block
     * operation - EDGE_WEIGHT_SET or EDGE_WEIGHT_INCR
     */
    template<int64_t direction, class use_dest>
    static void
    do_edge_updates(
        int64_t result, bool create, iterator pos, iterator updates_end,
        stinger_t * G, stinger_eb *eb, size_t e, int64_t operation)
    {
        // 'pos' points to the first update for this edge slot; there may be more than one, or none if it equals 'updates_end'
        // Keep incrementing the iterator until we reach an update with a different neighbor
        for (iterator u = pos; u != updates_end && use_dest::get(*u) == use_dest::get(*pos); ++u) {
            int64_t neighbor = use_dest::get(*pos);
            int64_t weight = adapter::get_weight(*u);
            int64_t time = adapter::get_time(*u);
            if (u == pos) {
                adapter::set_result(*u, result);
                update_edge_data_and_direction (G, eb, e, neighbor, weight, time, direction, create ? EDGE_WEIGHT_SET : operation);
            } else {
                adapter::set_result(*u, EDGE_UPDATED);
                update_edge_data_and_direction (G, eb, e, neighbor, weight, time, direction, operation);
            }
        }
    }

    template<int64_t direction, class use_dest>
    static void
    do_edge_updates(
        int64_t result, iterator pos, iterator updates_end) {
            for (iterator u = pos; u != updates_end && use_dest::get(*u)==use_dest::get(*pos); ++u)
                adapter::set_result(*u, result);
        }
    template<int64_t direction, class use_dest>
    static void
    do_edge_updates(
        int64_t result, iterator pos, iterator updates_begin, bool reverse) {
            for (iterator u = pos; u != updates_begin-1 && use_dest::get(*u)==use_dest::get(*pos); --u)
                adapter::set_result(*u, result);
        }

    /*
     * The core algorithm for updating edges in one direction.
     * Similar to stinger_update_directed_edge(), but optimized to perform several updates for the same vertex.
     * direction - are we updating out-edges or in-edges?
     * use_dest - Either source_funcs or dest_funcs (allows us to swap source/destination easily)
     */
    template<int64_t direction, class use_dest>
    static inline void
    update_directed_edges_for_vertex_binary(  // shared: all binary, index: all binary
            SharedBuffer *shared, VertexArray * va, int64_t src, int64_t type,
            iterator updates_begin,
            iterator updates_end,
            int64_t operation)
    {
        assert(direction == STINGER_EDGE_DIRECTION_OUT || direction == STINGER_EDGE_DIRECTION_IN);
        // Track the next edge that should be inserted
        next_update_tracker next_update(updates_begin, updates_end);
        Edge e;  e.src = src;

        if (direction == STINGER_EDGE_DIRECTION_OUT) {
            // 1.first step: find, all binary search
            if (get_outdegree(va, src) <= THRESHOLD) {
                int64_t outid = va->vertices[src].out_buffer_index;
                if (outid == -1)  {
                    outid = shared->allocate_room(src);
                    va->vertices[src].out_buffer_index = outid;
                }
                for (int64_t i = outid * THRESHOLD; i < (outid + 1) * THRESHOLD; ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }
            } else if (get_outdegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_outdata(va, src);
                for (int i = 0; i < rect->index_used; ++i) {
                    for (int j = 0; j < rect->arr_size; ++j) {
                        int64_t dest = rect->first_index[i][j].dst;
                        if (dest == -1)  break;
                        iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                        do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                    }
                }

                if (rect->next != nullptr) {
                    for (int i = 0; i < rect->next->index_used; ++i) {
                        for (int j = 0; j < rect->next->arr_size; ++j) {
                            int64_t dest = rect->next->first_index[i][j].dst;
                            if (dest == -1)  break;
                            iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                        }
                    }
                }
            }

            // 2. second step: insert
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t outid = va->vertices[src].out_buffer_index;
                if (get_outdegree(va, src) < THRESHOLD) {
                    shared->push_room(outid, e, va, direction);
                    increment_outdegree(va, src);
                } else if (get_outdegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].out_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_outdata(va, src, indirect);
                    }
                    Edge tmp_e;
                    IndirectIndex * rect = get_outdata(va, src);
                    for (int64_t i = outid * THRESHOLD;
                        i < ((outid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                    //    rect->do_push(tmp_e.dst);
                        rect->push_append(tmp_e.dst);
                        shared->set_invalid(&shared->buffer[i]);
                    }
                //    rect->do_push(e.dst);
                    rect->push_append(e.dst);
                    increment_outdegree(va, src);
                    va->vertices[src].out_buffer_index = -1;
                    // clear-room-bit
                    shared->bit_buffer->clear_bit(outid);
                } else {
                    IndirectIndex * rect = get_outdata(va, src);
                    rect->push_nobeep(e.dst);
                    increment_outdegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }
        }
        
        if (direction == STINGER_EDGE_DIRECTION_IN) {
            if (get_indegree(va, src) <= THRESHOLD) {
                // find
                int64_t inid = va->vertices[src].in_buffer_index;
                if (inid == -1)  {
                    inid = shared->allocate_room(src);
                    va->vertices[src].in_buffer_index = inid;
                }
                for (int64_t i = inid * THRESHOLD; i < (inid + 1) * THRESHOLD; ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }
            } else if (get_indegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_indata(va, src);
                for (int i = 0; i < rect->index_used; ++i) {
                    for (int j = 0; j < rect->arr_size; ++j) {
                        int64_t dest = rect->first_index[i][j].dst;
                        if (dest == -1)  break;
                        iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                        do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                    }
                }

                if (rect->next != nullptr) {
                    for (int i = 0; i < rect->next->index_used; ++i) {
                        for (int j = 0; j < rect->next->arr_size; ++j) {
                            int64_t dest = rect->next->first_index[i][j].dst;
                            if (dest == -1)  break;
                            iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                        }
                    }
                }
            }

            // 2.second step: insert
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t inid = va->vertices[src].in_buffer_index;
                if (get_indegree(va, src) < THRESHOLD) {
                    shared->push_room(inid, e, va, direction);
                    increment_indegree(va, src);
                } else if (get_indegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].in_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_indata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_indata(va, src);
                    for (int64_t i = inid * THRESHOLD;
                        i < ((inid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                        //rect->do_push(tmp_e.dst);
                        rect->push_append(tmp_e.dst);
                        shared->set_invalid(&shared->buffer[i]);
                    }
                    //rect->do_push(e.dst);
                    rect->push_append(e.dst);
                    increment_indegree(va, src);
                    va->vertices[src].in_buffer_index = -1;
                    // clear-room-bit
                    shared->bit_buffer->clear_bit(inid);
                } else {
                    IndirectIndex * rect = get_indata(va, src);
                    rect->push_nobeep(e.dst);
                    increment_indegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }
        }
    }

// our-cen : change init code to centralized processing search and insert
template<int64_t direction, class use_dest>
    static inline void
    update_directed_edges_for_vertex_beep(   // shared: all binary, index: all beep
            SharedBuffer *shared, VertexArray * va, int64_t src, int64_t type,
            iterator updates_begin,
            iterator updates_end,
            int64_t operation)
    {
        assert(direction == STINGER_EDGE_DIRECTION_OUT || direction == STINGER_EDGE_DIRECTION_IN);
        // Track the next edge that should be inserted
        next_update_tracker next_update(updates_begin, updates_end);
        Edge e;  e.src = src;
        iterator pos = updates_begin;
        int64_t distance = std::distance(updates_begin, updates_end);

        if (direction == STINGER_EDGE_DIRECTION_OUT) {
            // 1.first step: find
            if (get_outdegree(va, src) <= THRESHOLD) {
                int64_t outid = va->vertices[src].out_buffer_index;
                if (outid == -1)  {
                    outid = shared->allocate_room(src);
                    va->vertices[src].out_buffer_index = outid;
                }
                //binary search
                for (int64_t i = (outid * THRESHOLD); i < ((outid + 1) * THRESHOLD); ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }
            } else if (get_outdegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_outdata(va, src);
                // beep search
                while (pos != updates_end) {
                    if (rect->do_find(use_dest::get(*pos)))
                        do_edge_updates<direction, use_dest>(EDGE_UPDATED,  pos, updates_end);
                    ++pos;
                }
            }

            // 2. second step: insert
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t outid = va->vertices[src].out_buffer_index;
                if (get_outdegree(va, src) < THRESHOLD) {
                    shared->push_room(outid, e, va, direction);
                    increment_outdegree(va, src);
                } else if (get_outdegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].out_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_outdata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_outdata(va, src);
                    for (int64_t i = (outid * THRESHOLD);
                        i < ((outid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                        rect->do_push(tmp_e.dst);
                        shared->set_invalid(&shared->buffer[i]);
                    }
                // room_clear_bit
                    shared->bit_buffer->clear_bit(outid);
                    rect->do_push(e.dst);
                    increment_outdegree(va, src);  
                    va->vertices[src].out_buffer_index = -1;                               
                } else {
                    IndirectIndex * rect = get_outdata(va, src);
                    rect->do_push(e.dst);
                    increment_outdegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }
        }
        
        if (direction == STINGER_EDGE_DIRECTION_IN) {                
            // 1. find
            if (get_indegree(va, src) <= THRESHOLD) {
                int64_t inid = va->vertices[src].in_buffer_index;
                if (inid == -1)  {
                    inid = shared->allocate_room(src);
                    va->vertices[src].in_buffer_index = inid;
                }
                // binary search
                for (int64_t i = (inid * THRESHOLD); i < ((inid + 1) * THRESHOLD); ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }
            } else if (get_indegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_indata(va, src);
                // beep search
                while (pos != updates_end) {
                    if (rect->do_find(use_dest::get(*pos)))
                        do_edge_updates<direction, use_dest>(EDGE_UPDATED,  pos, updates_end);
                    ++pos;
                }
            }

            // 2.second step: insert
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t inid = va->vertices[src].in_buffer_index;
                if (get_indegree(va, src) < THRESHOLD) {
                    shared->push_room(inid, e, va, direction);
                    increment_indegree(va, src);
                } else if (get_indegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].in_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_indata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_indata(va, src);
                    for (int64_t i = (inid * THRESHOLD);
                        i < ((inid + 1) * THRESHOLD); ++i) {                                                 
                        tmp_e.dst = shared->buffer[i].dst;
                        rect->do_push(tmp_e.dst);
                        shared->set_invalid(&shared->buffer[i]);
                    }
                // room_clear_bit
                    shared->bit_buffer->clear_bit(inid);
                    rect->do_push(e.dst);
                    increment_indegree(va, src);
                    va->vertices[src].in_buffer_index = -1;
                } else {
                    IndirectIndex * rect = get_indata(va, src);
                    rect->do_push(e.dst);
                    increment_indegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }  
        }
    }

// our-optimization : dynamic choose search mode
template<int64_t direction, class use_dest>
    static inline void
    update_directed_edges_for_vertex_hybrid(   // index: binary/beep
            SharedBuffer *shared, VertexArray * va, int64_t src, int64_t type,
            iterator updates_begin,
            iterator updates_end,
            int64_t operation)
    {
        assert(direction == STINGER_EDGE_DIRECTION_OUT || direction == STINGER_EDGE_DIRECTION_IN);
        // Track the next edge that should be inserted
        next_update_tracker next_update(updates_begin, updates_end);
        Edge e;  e.src = src;
        iterator pos = updates_begin;
        int64_t distance = std::distance(updates_begin, updates_end);
        int64_t benefit = 0;

        if (direction == STINGER_EDGE_DIRECTION_OUT) {
            // 1.first step: find
            if (get_outdegree(va, src) <= THRESHOLD) {
                int64_t outid = va->vertices[src].out_buffer_index;
                if (outid == -1)  {
                    outid = shared->allocate_room(src);
                    va->vertices[src].out_buffer_index = outid;
                }
                // binary search
                for (int64_t i = outid * THRESHOLD; i < (outid + 1) * THRESHOLD; ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }

            } else if (get_outdegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_outdata(va, src);
                int64_t degree = get_outdegree(va,src);
                benefit = (degree*log2(distance)) - (2*distance*sqrt(degree));
                if (benefit > 0) {
                    // beep search
                    if (rect->heap_height == 0) {  // the first time to find in beep, need make_heap
                        // call make_beap;
                        rect->make_beap();
                    }
                    while (pos != updates_end) {
                        if (rect->do_find(use_dest::get(*pos)))
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  pos, updates_end);
                        ++pos;
                    }
                } else {
                    // binary search
                    for (int i = 0; i < rect->index_used; ++i) {
                        for (int j = 0; j < rect->arr_size; ++j) {
                            int64_t dest = rect->first_index[i][j].dst;
                            if (dest == -1)  break;
                            iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                        }
                    }
                    if (rect->next != nullptr) {
                        for (int i = 0; i < rect->next->index_used; ++i) {
                            for (int j = 0; j < rect->next->arr_size; ++j) {
                                int64_t dest = rect->next->first_index[i][j].dst;
                                if (dest == -1)  break;
                                iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                            }
                        }
                    }
                }   
            }

            // 2. second step: insert
        //    LOG_I_A("start to insert");
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t outid = va->vertices[src].out_buffer_index;
                if (get_outdegree(va, src) < THRESHOLD) {
                    shared->push_room(outid, e, va, direction);
                    increment_outdegree(va, src);
                } else if (get_outdegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].out_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_outdata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_outdata(va, src);
                    for (int64_t i = outid * THRESHOLD;
                        i < ((outid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                        rect->push_append(tmp_e.dst);  
                        shared->set_invalid(&shared->buffer[i]);
                    }
                    // room_clear_bit
                    shared->bit_buffer->clear_bit(outid);
                    rect->push_append(e.dst);
                    increment_outdegree(va, src);   
                    va->vertices[src].out_buffer_index = -1;                              
                } else {
                    IndirectIndex * rect = get_outdata(va, src);
                    if (rect->heap_height != 0) {
                        rect->do_push(e.dst);
                    } else {
                        rect->push_append(e.dst);
                    }        
                    increment_outdegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }
        }
        
        if (direction == STINGER_EDGE_DIRECTION_IN) {
            if (get_indegree(va, src) <= THRESHOLD) {
                // find
                int64_t inid = va->vertices[src].in_buffer_index;
                if (inid == -1)  {
                    inid = shared->allocate_room(src);
                    va->vertices[src].in_buffer_index = inid;
                }
                // binary search
                for (int64_t i = inid * THRESHOLD; i < (inid + 1) * THRESHOLD; ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }
                
            } else if (get_indegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_indata(va, src);
                int64_t degree = get_indegree(va,src);
                benefit = (degree*log2(distance)) - (2*distance*sqrt(degree));
                if (benefit > 0) {
                    // beep search
                    if (rect->heap_height == 0) {  // the first time to find in beep, need make_heap
                        // call make_beap;
                        rect->make_beap();
                    }
                    while (pos != updates_end) {
                        if (rect->do_find(use_dest::get(*pos)))
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  pos, updates_end);
                        ++pos;
                    }
                } else {
                    // binary search
                    for (int i = 0; i < rect->index_used; ++i) {
                        for (int j = 0; j < rect->arr_size; ++j) {
                            int64_t dest = rect->first_index[i][j].dst;
                            if (dest == -1)  break;
                            iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                        }
                    }
                    if (rect->next != nullptr) {
                        for (int i = 0; i < rect->next->index_used; ++i) {
                            for (int j = 0; j < rect->next->arr_size; ++j) {
                                int64_t dest = rect->next->first_index[i][j].dst;
                                if (dest == -1)  break;
                                iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                            }
                        }
                    } 
                }   
            }

            // 2.second step: insert
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t inid = va->vertices[src].in_buffer_index;
                if (get_indegree(va, src) < THRESHOLD) {
                    shared->push_room(inid, e, va, direction);
                    increment_indegree(va, src);
                } else if (get_indegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].in_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_indata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_indata(va, src);
                    for (int64_t i = inid * THRESHOLD;
                        i < ((inid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                        rect->push_append(tmp_e.dst);
                        shared->set_invalid(&shared->buffer[i]);
                    }
                    // room_clear_bit
                    shared->bit_buffer->clear_bit(inid);
                    rect->push_append(e.dst);
                    increment_indegree(va, src);
                    va->vertices[src].in_buffer_index = -1;
                } else {
                    IndirectIndex * rect = get_indata(va, src);
                    if (rect->heap_height != 0) {
                        rect->do_push(e.dst);
                    } else {
                        rect->push_append(e.dst);
                    }                       
                    increment_indegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }  
        }
    }

    template<int64_t direction, class use_dest>
    static inline void
    update_directed_edges_for_vertex_hybrid2(   // shared:binary, index: binary/beep
            SharedBuffer *shared, VertexArray * va, int64_t src, int64_t type,
            iterator updates_begin,
            iterator updates_end,
            int64_t operation)
    {
        assert(direction == STINGER_EDGE_DIRECTION_OUT || direction == STINGER_EDGE_DIRECTION_IN);
        // Track the next edge that should be inserted
        next_update_tracker next_update(updates_begin, updates_end);
        Edge e;  e.src = src;
        iterator pos = updates_begin;
        int64_t distance = std::distance(updates_begin, updates_end);

        if (direction == STINGER_EDGE_DIRECTION_OUT) {
            // 1.first step: find
            if (get_outdegree(va, src) <= THRESHOLD) {
                int64_t outid = va->vertices[src].out_buffer_index;
                if (outid == -1)  {
                    outid = shared->allocate_room(src);
                    va->vertices[src].out_buffer_index = outid;
                }
                // binary search
                for (int64_t i = outid * THRESHOLD; i < (outid + 1) * THRESHOLD; ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }

            } else if (get_outdegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_outdata(va, src);
                if (rect == nullptr)  exit(-1);
                int64_t degree = get_outdegree(va,src);
            //    if (degree*log2(distance) - distance*sqrt(degree) > 0) {
                if (degree*log2(distance)+distance - 2*distance*sqrt(degree) > 0) {
            //    if (degree*log2(distance)+1 - 2*distance*sqrt(degree) > 0) {
                    // beep search
                    while (pos != updates_end) {
                        if (rect->do_find(use_dest::get(*pos)))
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  pos, updates_end);
                        ++pos;
                    }
                } else {
                    // binary search
                    for (int i = 0; i < rect->index_used; ++i) {
                        for (int j = 0; j < rect->arr_size; ++j) {
                            int64_t dest = rect->first_index[i][j].dst;
                            if (dest == -1)  break;
                            iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                        }
                    }
                }   
            }

            // 2. second step: insert
        //    LOG_I_A("start to insert");
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t outid = va->vertices[src].out_buffer_index;
                if (get_outdegree(va, src) < THRESHOLD) {
                    shared->push_room(outid, e, va, direction);
                    increment_outdegree(va, src);
                } else if (get_outdegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].out_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_outdata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_outdata(va, src);
                    for (int64_t i = outid * THRESHOLD;
                        i < ((outid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                        rect->do_push(tmp_e.dst);  // must change to append
                        shared->set_invalid(&shared->buffer[i]);
                    }
                    // room_clear_bit
                    shared->bit_buffer->clear_bit(outid);
                    rect->do_push(e.dst);
                    increment_outdegree(va, src);   
                    va->vertices[src].out_buffer_index = -1;                              
                } else {
                    IndirectIndex * rect = get_outdata(va, src);
                    rect->do_push(e.dst);
                    increment_outdegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }
        }
        
        if (direction == STINGER_EDGE_DIRECTION_IN) {
            if (get_indegree(va, src) <= THRESHOLD) {
                // find
                int64_t inid = va->vertices[src].in_buffer_index;
                if (inid == -1)  {
                    inid = shared->allocate_room(src);
                    va->vertices[src].in_buffer_index = inid;
                }
                // binary search
                for (int64_t i = inid * THRESHOLD; i < (inid + 1) * THRESHOLD; ++i) {
                    int64_t dest = shared->buffer[i].dst;
                    if (dest == -1)  break;
                    iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                    do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                }

            } else if (get_indegree(va, src) > THRESHOLD) {
                IndirectIndex * rect = get_indata(va, src);
                int64_t degree = get_indegree(va,src);
            //    if (degree*log2(distance) - distance*sqrt(degree) > 0) {
                if (degree*log2(distance)+distance - 2*distance*sqrt(degree) > 0) {
            //    if (degree*log2(distance)+1 - 2*distance*sqrt(degree) > 0) {
                    // beep search
                    while (pos != updates_end) {
                        if (rect->do_find(use_dest::get(*pos)))
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  pos, updates_end);
                        ++pos;
                    }
                } else {
                    // binary search
                    for (int i = 0; i < rect->index_used; ++i) {
                        for (int j = 0; j < rect->arr_size; ++j) {
                            int64_t dest = rect->first_index[i][j].dst;
                            if (dest == -1)  break;
                            iterator u = find_updates<use_dest>(next_update(), updates_end, dest);
                            do_edge_updates<direction, use_dest>(EDGE_UPDATED,  u, updates_end);
                        }
                    }
                }   
            }

            // 2.second step: insert
            while (next_update() != updates_end) {
                e.dst = use_dest::get(*next_update());
                int64_t inid = va->vertices[src].in_buffer_index;
                if (get_indegree(va, src) < THRESHOLD) {
                    shared->push_room(inid, e, va, direction);
                    increment_indegree(va, src);
                } else if (get_indegree(va, src) == THRESHOLD) {
                    if (va->vertices[src].in_data == nullptr) {
                        IndirectIndex * indirect = new IndirectIndex();
                        set_indata(va, src, indirect);
                    }

                    Edge tmp_e;
                    IndirectIndex * rect = get_indata(va, src);
                    for (int64_t i = inid * THRESHOLD;
                        i < ((inid + 1) * THRESHOLD); ++i) {                           
                        tmp_e.dst = shared->buffer[i].dst;
                        rect->do_push(tmp_e.dst);
                        shared->set_invalid(&shared->buffer[i]);
                    }
                    // room_clear_bit
                    shared->bit_buffer->clear_bit(inid);
                    rect->do_push(e.dst);
                    increment_indegree(va, src);
                    va->vertices[src].in_buffer_index = -1;
                } else {
                    IndirectIndex * rect = get_indata(va, src);
                    rect->do_push(e.dst);
                    increment_indegree(va, src);
                }
                do_edge_updates<direction, use_dest>(EDGE_UPDATED,  next_update(), updates_end);
            }  
        }
    }

    template<int64_t direction, class use_dest>
    static inline void
    update_directed_edges_for_vertex(
            SharedBuffer *shared, VertexArray * va, int64_t src, int64_t type,
            iterator updates_begin,
            iterator updates_end,
            int64_t operation)
    {
        assert(direction == STINGER_EDGE_DIRECTION_OUT || direction == STINGER_EDGE_DIRECTION_IN);
        // Track the next edge that should be inserted
        next_update_tracker next_update(updates_begin, updates_end);

        if (va == nullptr)
            LOG_D("this static vertex_array is wrong");

        Edge e;
        e.src = src;

        int in_id = 0, in_degree = 0;
        int out_id = 0, out_degree = 0;

        if (direction == STINGER_EDGE_DIRECTION_OUT) {
            out_id = va->vertices[src].out_buffer_index;
            iterator pos = updates_begin;

             while (pos != updates_end) {
                e.dst = use_dest::get(*pos);
            //    if (get_outdegree(va, src) < THRESHOLD && get_outdata(va, src) == nullptr) {
                if (get_outdegree(va, src) < THRESHOLD) {
                    if (out_id == -1) {
                    //    out_id = va->vertices[src].out_buffer_index;
                        if (out_id == -1) {
                            out_id = shared->allocate_room(src);
                        }
                        va->vertices[src].out_buffer_index = out_id;
                        //__sync_val_compare_and_swap(&va->vertices[src].out_buffer_index, -1, out_id);
                    }
                    
                    int ret = shared->push_room(out_id, e, va, direction);                   
                    if (ret == 1)
                        increment_outdegree(va, src);

            //    } else if (get_outdegree(va, src) == THRESHOLD && out_id != -1) {
                } else if (get_outdegree(va, src) == THRESHOLD) {
                    if (!(shared->find_in_room(out_id, e.dst))) {
                     //   Edge move[THRESHOLD+1];
                     //   shared->move_to_index(move, out_id);

                     //   move[THRESHOLD].src = src;
                     //   move[THRESHOLD].dst = e.dst;

                        if (va->vertices[src].out_data == nullptr) {
                            IndirectIndex * indirect = new IndirectIndex();
                            set_outdata(va, src, indirect); 
                            // for (int i = 0; i <= THRESHOLD; ++i)
                            //     indirect->do_push(move[i].dst);
                        } else {  // if another thread have done move operation
                            IndirectIndex * rect = get_outdata(va, src);
                            //IndirectIndex * rect = (IndirectIndex*)va->vertices[src].out_data;
                          //  rect->do_push(e.dst);
                        }

                        Edge tmp_e;
                        IndirectIndex * rect = get_outdata(va, src);
                        for (size_t i = out_id * THRESHOLD;
                            i < size_t((out_id + 1) * THRESHOLD); ++i) {
                            //tmp_e.src = shared->buffer[i].src;
                            tmp_e.dst = shared->buffer[i].dst;
                            rect->do_push(tmp_e.dst);
                            shared->set_invalid(&shared->buffer[i]);
                        }
                        rect->do_push(e.dst);

                        // if (va->vertices[src].out_data == nullptr)
                        //     printf("move here, out data is null, null out data\n");
                        increment_outdegree(va, src);
                        va->vertices[src].out_buffer_index = -1;
                    }

                } else if (get_outdegree(va, src) > THRESHOLD) {  // && out_id == -1) {
                    //IndirectIndex * rect = (IndirectIndex*)va->vertices[src].out_data;
                    IndirectIndex * rect = get_outdata(va, src);
                    
                    // if (rect == nullptr)
                    //     printf("out data is null, null out data\n");
                   
                    // int index_id = 0, arr_id = 0;
                    // if (!rect->find_nobeep(e.dst, index_id, arr_id)) {
                    //     rect->push_nobeep(e.dst, index_id, arr_id);
                    //     increment_outdegree(va, src);
                    // }
                    if (!(rect->do_find(e.dst))) {
                        rect->do_push(e.dst);
                        increment_outdegree(va, src);
                    }
                }
               ++pos;
            }

        // -----------------------change the direction-----------------------------------
        } else if (direction == STINGER_EDGE_DIRECTION_IN) {
            in_id = va->vertices[src].in_buffer_index;
            iterator pos = updates_begin;

            while (pos != updates_end) {
                e.dst = use_dest::get(*pos);

                if (get_indegree(va, src) < THRESHOLD) {  // && get_indata(va, src) == nullptr) {
                    if (in_id == -1) {
                        in_id = va->vertices[src].in_buffer_index;
                        if (in_id == -1) {
                            in_id = shared->allocate_room(src);
                        }
                            
                        if (in_id == -1)
                            LOG_D("there is no empty room in this whole buffer");
                        //__sync_val_compare_and_swap(&va->vertices[src].in_buffer_index, -1, in_id);
                        va->vertices[src].in_buffer_index = in_id;
                    }
                    int ret = shared->push_room(in_id, e, va, direction);
                    
                    if (ret == 1)
                        increment_indegree(va, src);
                    //    __sync_add_and_fetch(&va->vertices[src].in_degree, 1);

                } else if (get_indegree(va, src) == THRESHOLD) {  // && in_id != -1) {
                    if (!(shared->find_in_room(in_id, e.dst))) {
                        // really a new edge, need remove and then insert this new edge
                        Edge move[THRESHOLD+1];
                        shared->move_to_index(move, in_id);

                        move[THRESHOLD].src = src;
                        move[THRESHOLD].dst = e.dst;

                        if (va->vertices[src].in_data == nullptr) {
                            IndirectIndex * indirect = new IndirectIndex();
                            set_indata(va, src, indirect);
                        //    va->vertices[src].in_data = indirect;
                            for (int i = 0; i <= THRESHOLD; ++i)
                                indirect->do_push(move[i].dst);
                        } else {
                            IndirectIndex * rect = get_indata(va, src);
                            rect->do_push(e.dst);
                        }
                        if (va->vertices[src].in_data == nullptr)
                            LOG_D("move here, out data is null, null out data");
                    //    indirect->move_into_indirect(move, THRESHOLD+1);

                        increment_indegree(va, src);
                        va->vertices[src].in_buffer_index = -1;
                    }

                } else if (get_indegree(va, src) > THRESHOLD) {  // && in_id == -1) {
                    IndirectIndex * rect = get_indata(va, src);
                    // if (rect == nullptr)
                    //     LOG_D("out data is null, null out data");
                    // int index_id = 0, arr_id = 0;
                    // if (!rect->find_nobeep(e.dst, index_id, arr_id)) {
                    //     rect->push_nobeep(e.dst, index_id, arr_id);
                    //     increment_indegree(va, src);
                    // }
                    if (!(rect->do_find(e.dst))) {
                        rect->do_push(e.dst);
                        increment_indegree(va, src);
                    }
                }
                ++pos;
            }
        }
    }

    template <class use_source>
    static bool same_source_and_type(const iterator &a, const iterator &b)
    {
        return adapter::get_type(*a) == adapter::get_type(*b)
            && use_source::equals(*a, *b);
    }

    /*
     * Splits a range of updates into chunks that all update the same source,
     * then calls update_directed_edges_for_vertex() in parallel on each range
     *
     * Template arguments:
     * direction - are we updating out-edges or in-edges?
     * use_source - set of functions to use for working with the source vertex (source_funcs or dest_funcs)
     * use_dest - set of functions to use for working with the destination vertex (source_funcs or dest_funcs)
     */

    // @re-write
    template<int64_t direction, class use_source, class use_dest>
    static void
    do_batch_update( SharedBuffer * shared, VertexArray * va, iterator updates_begin, iterator updates_end, int64_t operation)
    {
        typedef typename std::vector<iterator>::iterator iterator_ptr;
        typedef typename std::pair<iterator, iterator> range;
        typedef typename std::vector<range>::iterator range_iterator;

        // Sort by type, src, dst, time ascending
 //       LOG_V("Sorting...");
        std::sort(updates_begin, updates_end, use_source::sort);

        // Get a list of pointers to each element
        int64_t num_updates = std::distance(updates_begin, updates_end);
        std::vector<iterator> pointers(num_updates);
     //   int size_av = num_updates/core_num;
        OMP("omp parallel for")
        for (int64_t i = 0; i < num_updates; ++i) { pointers[i] = updates_begin + i; }

        std::vector<iterator> unique_sources(num_updates);   // 拷贝指定范围的唯一化
        iterator_ptr last_unique_source = std::unique_copy(
            pointers.begin(), pointers.end(), unique_sources.begin(), same_source_and_type<use_source>);
        unique_sources.erase(last_unique_source, unique_sources.end());
        // Now each consecutive pair of elements represents a range of updates for the same type and source ID
        unique_sources.push_back(updates_end);

    //     // Split up long ranges of updates for the same source vertex
    // //    LOG_V("Splitting ranges...");
        std::vector<range> update_ranges;
        iterator_ptr b_ptr = unique_sources.begin();
        iterator_ptr e_ptr = unique_sources.end();
        vector<vector<range>> local_ranges(core_num);
        OMP("omp parallel for")
        for (iterator_ptr ptr = b_ptr; ptr < e_ptr-1; ++ptr)
        {
            iterator begin = *ptr;
            iterator end = *(ptr + 1);

            // Calculate number of updates for each range
            size_t num_updates = std::distance(begin, end);
            size_t num_ranges = core_num;
            size_t updates_per_range = std::floor((double)num_updates / num_ranges);   // floor向下取整

        //    std::vector<range> local_ranges; 
            int k = omp_get_thread_num();
            local_ranges[k].push_back(make_pair(begin, end));
        }
        
    //    omp_set_num_threads(1);
        for (int i = 0; i < core_num; ++i) {
            update_ranges.insert(update_ranges.end(), local_ranges[i].begin(), local_ranges[i].end());
        }

        // LOG_I_A("Entering parallel update loop: %ld updates for %ld vertices.",
        //     std::distance(updates_begin, updates_end), unique_sources.size()-1);

        range_iterator range_begin = update_ranges.begin();
        range_iterator range_end = update_ranges.end();
        size_t range_size = update_ranges.size();
    //    LOG_I_A("the range size is %ld \n", range_size);
        
      if (mode == 0) {
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < range_size; ++i)
        {
            // Get this thread's range of updates
            iterator begin = update_ranges[i].first;
            iterator end = update_ranges[i].second;
            // Each range guaranteed to have same edge type and source vertex
            int64_t type = adapter::get_type(*begin);
            int64_t source = use_source::get(*begin);
        //    LOG_I_A("Thread %d processing %ld updates (%ld -> *)", omp_get_thread_num(), std::distance(begin, end), source);
            update_directed_edges_for_vertex_binary<direction, use_dest>(shared, va, source, type, begin, end, operation);
        }
      } else if (mode == 1) {
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < range_size; ++i)
        {
            iterator begin = update_ranges[i].first;
            iterator end = update_ranges[i].second;
            int64_t type = adapter::get_type(*begin);
            int64_t source = use_source::get(*begin);
            update_directed_edges_for_vertex_beep<direction, use_dest>(shared, va, source, type, begin, end, operation);
        }
      } else if (mode == 2) {
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < range_size; ++i)
        {
            iterator begin = update_ranges[i].first;
            iterator end = update_ranges[i].second;
            int64_t type = adapter::get_type(*begin);
            int64_t source = use_source::get(*begin);
            update_directed_edges_for_vertex_hybrid<direction, use_dest>(shared, va, source, type, begin, end, operation);
        }
      } else {
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < range_size; ++i)
        {
            iterator begin = update_ranges[i].first;
            iterator end = update_ranges[i].second;
            int64_t type = adapter::get_type(*begin);
            int64_t source = use_source::get(*begin);
            update_directed_edges_for_vertex_hybrid2<direction, use_dest>(shared, va, source, type, begin, end, operation);
        }
      }       
    }

    static bool
    source_less_than_destination(const update &x){
        return adapter::get_source(x) < adapter::get_dest(x);
    }

    // We use the result field to keep track of which updates have been performed
    // Between IN and OUT iterations, we clear it out to make sure all the updates are performed again
    static void

    clear_results(iterator begin, iterator end)
    {
        OMP("omp parallel for")
        for (iterator u = begin; u < end; ++u) {
            int64_t result = adapter::get_result(*u);
            switch (result){
                case EDGE_ADDED:
                case EDGE_UPDATED:
                case EDGE_NOT_ADDED:
                    // Reset return code so we process this update next iteration
                    adapter::set_result(*u, PENDING);
                    break;
                case EDGE_ADDED-result_code_offset:
                case EDGE_NOT_ADDED-result_code_offset:
                    // Caller must have set the return code, leave it be
                    break;
                case PENDING:
                default:
                    // We didn't finish all the updates, or invalid codes were set
                    assert(0);

            }
        }
    }

    static void
    remap_results(iterator begin, iterator end)
    {
        OMP("omp parallel for")
        for (iterator u = begin; u < end; ++u) {
            int64_t result = adapter::get_result(*u);
            switch (result){
                case EDGE_ADDED:
                case EDGE_UPDATED:
                case EDGE_NOT_ADDED:
                    // Subtract 10 to get back to return code caller expects
                    adapter::set_result(*u, result - result_code_offset);
                    break;
                case EDGE_ADDED-result_code_offset:
                case EDGE_NOT_ADDED-result_code_offset:
                    // Caller must have set the return code, leave it be
                    break;
                case PENDING:
                default:
                    // We didn't finish all the updates, or invalid codes were set
                    assert(0);

            }
        }
    }

    // @re-write
    static void
    batch_update_dispatch(SharedBuffer * shared, VertexArray * va, iterator updates_begin, iterator updates_end, int64_t operation, bool directed)
    {
        LOG_V("Partitioning batch to obey ordering constraints...");
        const iterator pos = std::partition(updates_begin, updates_end, source_less_than_destination);

        const int64_t OUT = STINGER_EDGE_DIRECTION_OUT;
        const int64_t IN = STINGER_EDGE_DIRECTION_IN;

        // All elements between begin and pos have src < dest. Update the out-edge slot first
        LOG_V("Beginning OUT updates for first half of batch...");
        do_batch_update<OUT, source_funcs, dest_funcs>(shared, va, updates_begin, pos, operation);
        clear_results(updates_begin, pos);

        LOG_V("Beginning IN updates for first half of batch...");
        do_batch_update<IN, dest_funcs, source_funcs>(shared, va, updates_begin, pos, operation);

        // All elements between pos and end have src > dest. Update the in-edge slot first
        LOG_V("Beginning IN updates for second half of batch...");
        do_batch_update<IN, dest_funcs, source_funcs>(shared, va, pos, updates_end, operation);
        clear_results(pos, updates_end);
        LOG_V("Beginning OUT updates for second half of batch...");
        do_batch_update<OUT, source_funcs, dest_funcs>(shared, va, pos, updates_end, operation);
    }
}; // end of class batch insert

}}; // end namespace gt::stinger

template<typename adapter, typename iterator>
// re-write
void
stinger_batch_incr_edges(SharedBuffer * shared, VertexArray * va, iterator begin, iterator end)
{
    gt::stinger::BatchInserter<adapter, iterator>::batch_update_dispatch(shared,  va, begin, end, EDGE_WEIGHT_INCR, true);
}
template<typename adapter, typename iterator>
void
stinger_batch_insert_edges(stinger_t * G, iterator begin, iterator end)
{
    gt::stinger::BatchInserter<adapter, iterator>::batch_update_dispatch(G, begin, end, EDGE_WEIGHT_SET, true);
}

// @re-write
template<typename adapter, typename iterator>
void
stinger_batch_incr_edge_pairs(SharedBuffer * shared, VertexArray * va, iterator begin, iterator end)
{
    gt::stinger::BatchInserter<adapter, iterator>::batch_update_dispatch(shared, va, begin, end, EDGE_WEIGHT_INCR, false);
}

template<typename adapter, typename iterator>
void
stinger_batch_insert_edge_pairs(stinger_t * G, iterator begin, iterator end)
{
    gt::stinger::BatchInserter<adapter, iterator>::batch_update_dispatch(G, begin, end, EDGE_WEIGHT_SET, false);
}

#endif //STINGER_BATCH_INSERT_H_
