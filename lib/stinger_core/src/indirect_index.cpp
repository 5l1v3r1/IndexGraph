/*
 * Copyright 2018-now Zhihui Zhao
 * @Author: Zhihui Zhao
 * @Last Modified by: Yilong Liu
 */

#include <cstdio>
#include <cstdint>

#include "stinger_core/defs.h"

#include "stinger_core/bitmap.hpp"

#include "stinger_core/shared_buffer.h"
#include "stinger_core/indirect_index.h"
#include "stinger_core/stinger_atomics.h"

IndirectIndex::IndirectIndex() {
//    LOG_D("Entering IndirectIndex constructor and initial");
    this->index_used = 0;

    //this->index_size = (PAGE_SIZE-16)/sizeof(uint64_t*);
    this->arr_size = 1024/sizeof(Edge);

    //this->first_index = new Edge*[this->index_size];
    for (int i = 0; i < INDEX_SIZE; ++i) {
        this->first_index[i] = nullptr;
    }
    
    // add for beap
    this->heap_height = 0;
    this->heap_size = 0;
//    this->flag = false;

    next = nullptr;
}

IndirectIndex::IndirectIndex(int size) {
    this->index_used = 0;
    //this->index_size = (PAGE_SIZE-16)/sizeof(uint64_t*);
    this->arr_size = HUGE_PAGE/sizeof(Edge);

    //this->first_index = new Edge*[this->index_size];
    for (int i = 0; i < INDEX_SIZE; ++i) {
        this->first_index[i] = nullptr;
    }
    
    this->heap_height = 0;
    this->heap_size = 0;
//    this->flag = true;

    this->next = nullptr;
}

IndirectIndex::~IndirectIndex() {
 //   LOG_D("Entering IndirectIndex destructor");
    if (this->next != nullptr) {
        IndirectIndex * pnext = this->next;
        delete pnext;
    }

    for (int i = 0; i < this->index_used; ++i) {
        delete [] (this->first_index)[i];
    }

    //delete [] (this->first_index);
}

void IndirectIndex::set_array(Edge *arr) {
    for (int i = 0; i < this->arr_size; ++i) {
        arr[i].src = -1;
        arr[i].dst = -1;
    }
}

// void IndirectIndex::update_indirect(int64_t const &obj, int index_id, int arr_id) {
//    // this->first_index[index_id][arr_id].src = e.src;
//     this->first_index[index_id][arr_id].dst = obj;
// }

// -------------------------- no beep , sequence ------------------------------------------
bool IndirectIndex::find_in_indirect(int64_t const &obj, int &index_id, int &arr_id) {
    // only use to find a empty slot
    for (int i = 0; i < this->index_used; ++i) {
        for (int j = 0; j < this->arr_size; ++j) {
            // if (this->first_index[i][j].dst == obj) {
            // //    printf("case 1. this edge is already here in index[%d][%d]\n", i, j);
            //     index_id = i;
            //     arr_id = j;
            //     return true;
            // }

            if (i == this->index_used-1 && this->first_index[i][j].dst == -1) {
            //    printf("2. this edge is not found, "
            //    "but this array has empty slot in index[%d][%d]\n", i, j);
                index_id = i;
                arr_id = j;
                return false;
            }
        }
    }

//    printf("3. this edge is not found, "
//    "and there is no empty slot, need a new array\n");
    index_id = this->index_used;
    arr_id = this->arr_size;
    return false;
}

bool IndirectIndex::find_nobeep(int64_t const &obj, int &index_id, int &arr_id) {
    bool ret = this->find_in_indirect(obj, index_id, arr_id);
    if (this->next != nullptr) {
        ret = ret | this->next->find_in_indirect(obj, index_id, arr_id);
    }
    return ret;
}

//void IndirectIndex::push_indirect(int64_t const &obj, int index_id, int arr_id) {
void IndirectIndex::push_indirect(int64_t const &obj) {
    int index_id = 0, arr_id = 0;
    bool exist = find_in_indirect(obj, index_id, arr_id);

    if (arr_id != this->arr_size) {
        this->first_index[index_id][arr_id].dst = obj;        
    } else {
        if (index_id < INDEX_SIZE-1) {
            // because the [index_size-1] used to 2-level indirect
            this->first_index[index_id] = new Edge[arr_size];
            set_array(this->first_index[index_id]);
            this->first_index[index_id][0].dst = obj;
            ++this->index_used;
        }
    }
}

//void IndirectIndex::push_nobeep(int64_t const &obj, int index_id, int arr_id) {
void IndirectIndex::push_nobeep(int64_t const &obj) {
    //printf("******######*******######this insert use sequence\n");
    if ((this->index_used == INDEX_SIZE - 1) && (this->heap_size % this->arr_size) == 0) {
        // employ the secondnary indexing
        if (this->next == nullptr) {
            IndirectIndex * pnext = new IndirectIndex(INDEX_SIZE);
            next = pnext;
        }
        next->push_indirect(obj);
    } else {
        this->push_indirect(obj);
    }
}

// -----------------------------other function, move to index, show------------------------------
void IndirectIndex::move_into_indirect(Edge * e, int num) {
    int cnt = 0;
    while (num > this->arr_size) {
        this->first_index[this->index_used] = new Edge[this->arr_size];
        set_array(this->first_index[this->index_used]);
        ++this->index_used;
     //   __sync_add_and_fetch(&this->index_used, 1);

        for (int i = 0; i < this->arr_size; ++i) {
            this->first_index[index_used-1][i].src = e[cnt*arr_size+i].src;
            this->first_index[index_used-1][i].dst = e[cnt*arr_size+i].dst;
        }

        num -= this->arr_size;
        ++cnt;
    }

    if (num) {
        this->first_index[index_used] = new Edge[arr_size];
        set_array(this->first_index[this->index_used]);
        ++this->index_used;
    //    __sync_add_and_fetch(&this->index_used, 1);

        for (int i = 0; i < num; ++i) {
            this->first_index[this->index_used-1][i].src = e[cnt*arr_size+i].src;
            this->first_index[this->index_used-1][i].dst = e[cnt*arr_size+i].dst;
        }
    }
}

void IndirectIndex::show() {
    for (int i = 0; i < index_used; ++i) {
        printf("arr %d :", i);
        for (int j = 0; j < arr_size; ++j) {
            printf("<%ld,%ld>  ", first_index[i][j].src, first_index[i][j].dst);
            if (first_index[i][j].dst == -1)
                break;
        }
        printf("\n");
    }     
}

void IndirectIndex::do_show() {
    LOG_D("Showing the context of Indirect_Index");
    this->show();
    if (next != nullptr) {
        LOG_D("showing the next");
        this->next->show();
    }
}

// --------------------------------another push, just append, don't need traverse -----------------
void IndirectIndex::push_append(int64_t const &obj) {
    if ((this->index_used == INDEX_SIZE - 1) && (this->heap_size % this->arr_size) == 0) {
        // employ the secondnary indexing
        if (this->next == nullptr) {
            IndirectIndex * pnext = new IndirectIndex(INDEX_SIZE);
            next = pnext;
        }
    //    LOG_D("pushing into the next");
        next->do_push_append(obj);
    } else {
        this->do_push_append(obj);
    }
}

void IndirectIndex::do_push_append(int64_t const &obj) {
    if (heap_size == 0) {
        first_index[index_used] = new Edge[arr_size];
        set_array(first_index[index_used]);
        ++index_used;

        first_index[0][0].dst = obj;
        heap_size = 1;
        return;
    }

    int posn = heap_size;
    
    if (heap_size % arr_size == 0) {
        first_index[index_used] = new Edge[arr_size];
        set_array(first_index[index_used]);
        ++index_used;
    }
    
    // new object
    ++heap_size;
    
    first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = obj;
}

// ------------------------------  make_beap------------------------------------------------
/*bool IndirectIndex::cmp(int64_t a, int64_t b) {
	int x_a = a/arr_size;
	int y_a = a-x_a*arr_size;
	int x_b = b/arr_size;
	int y_b = b-x_b*arr_size;
	return first_index[x_a][y_a].dst > first_index[x_b][y_b].dst;
}*/

void IndirectIndex::quick_sort(int l, int r) {
    if (l < r) {
        int i = l, j = r, x = first_index[l/arr_size][l-(l/arr_size)*arr_size].dst;
        int x_i, y_i, x_j, y_j;
        while (i < j) {
            x_j = j/arr_size;  y_j = j-x_j*arr_size;
            x_i = i/arr_size;  y_i = i-x_i*arr_size;
            while (i < j && first_index[x_j][y_j].dst >= x) {
                j--;
                x_j = j/arr_size;  y_j = j-x_j*arr_size;
            }   
            if (i < j) {
                first_index[x_i][y_i].dst = first_index[x_j][y_j].dst;
                i++;
                x_i = i/arr_size;  y_i = i-x_i*arr_size;
            }
            while (i < j && first_index[x_i][y_i].dst < x) {
                i++;
                x_i = i/arr_size;  y_i = i-x_i*arr_size;
            }    
            if (i < j) {
                first_index[x_j][y_j].dst = first_index[x_i][y_i].dst;
                j--;
                x_j = j/arr_size;  y_j = j-x_j*arr_size;
            }
        }
        first_index[x_i][y_i].dst = x;
        quick_sort(l, i-1);
        quick_sort(i+1, r);
    }
}

void IndirectIndex::make_beap() {
    quick_sort(0, heap_size-1);
	int current = 1;
	while (current < heap_size) {
		++heap_height;
		current+=(heap_height+1);
	}
}

// --------------------------------use beep, search and find--------------------------------------
bool IndirectIndex::find(int64_t const &obj) {
    if (heap_size == 0)  return false;
    int h = heap_height;
    int posn = h*(h+1)/2;

    while (posn < heap_size && h >= 0) {
        int i = posn/arr_size;
        int j = posn - i*arr_size;
        if (obj > first_index[i][j].dst) {
            if (posn == (heap_height+1) * (heap_height+2)/2-1)
                return false;
            if (posn == heap_height*(heap_height+1)/2-1 && heap_size != (heap_height+1)*(heap_height+2)/2)
                return false;

            if (posn + h + 2 < heap_size) {
                posn += h+2;
                ++h;
            } else {
                if (posn-h >= 0 && posn+1 < heap_size)
                    ++posn;
                else
                    return false;
            }
        } else if (first_index[i][j].dst > obj) {
            if (posn == (heap_height+1) * (heap_height+2)/2-1)
                return false;
            if (posn == heap_height*(heap_height+1)/2-1 && heap_size != (heap_height+1)*(heap_height+2)/2)
                return false;
            else {
                if (posn == (h+1)*(h+2)/2-1)
                    return false;
                posn -= h;
                --h;
            }
        } else if (first_index[i][j].dst == obj)
            return true;
    }
    return false;
}

bool IndirectIndex::do_find(int64_t const &obj) {
    bool ret = this->find(obj);
    if (this->next != nullptr) {
        ret = ret | this->next->find(obj);
    }
    return ret;
}

void IndirectIndex::do_push(int64_t const &obj) {
    if ((this->index_used == INDEX_SIZE - 1) && (this->heap_size % this->arr_size) == 0) {
        // employ the secondnary indexing
        if (this->next == nullptr) {
            IndirectIndex * pnext = new IndirectIndex(INDEX_SIZE);
            next = pnext;
        }
    //    LOG_D("pushing into the next");
        next->push(obj);
    } else {
        this->push(obj);
    }
}

void IndirectIndex::push(int64_t const &obj) {
    if (heap_size == 0) {
        first_index[index_used] = new Edge[this->arr_size];
        set_array(first_index[index_used]);
        ++index_used;
     //   __sync_add_and_fetch(&index_used, 1);

        first_index[0][0].dst = obj;
        heap_size = 1;
        heap_height = 0;
        return;
    }

    int posn = heap_size;
    if (posn == first(heap_height + 1)) {
        ++heap_height;
      //  __sync_add_and_fetch(&heap_height, 1);
    }
    if ((heap_size % this->arr_size) == 0) {
        first_index[index_used] = new Edge[this->arr_size];
        set_array(first_index[index_used]);
        ++index_used;
    }
    
    // new object
    ++heap_size;

    for (int h = heap_height; h >= 2; --h) {
        if (posn == (h*(h+1)/2))
            return left_push(obj, posn, h);
        else if (posn == ((h+1)*(h+2)/2-1))
            return right_push(obj, posn, h);
        else {
            //swap with the larger parent
            int left_parent = posn - h - 1;
            int right_parent = posn - h;
            int left_i = left_parent/arr_size;
            int left_j = left_parent - left_i * arr_size;
            int right_i = right_parent/arr_size;
            int right_j = right_parent - right_i * arr_size;
            int pos_i = posn/arr_size;
            int pos_j = posn - pos_i * arr_size;

            if (obj < first_index[left_i][left_j].dst && first_index[left_i][left_j].dst > first_index[right_i][right_j].dst) {
                first_index[pos_i][pos_j].dst = first_index[left_i][left_j].dst;
                posn = left_parent;
            } else if (obj < first_index[right_i][right_j].dst) {
                first_index[pos_i][pos_j].dst = first_index[right_i][right_j].dst;
                posn = right_parent;
            } else {
                first_index[pos_i][pos_j].dst = obj;
                return;
            }
        }
    }
    if (obj < first_index[0][0].dst) {
        first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = first_index[0][0].dst;
        first_index[0][0].dst = obj;
    } else {
        first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = obj;
    }
}

// the left edge of the beap
void IndirectIndex::left_push(int64_t const &obj, int posn, int height) {
    // swap if the parent is greater, otherwise, store
    for (int h = height; h >= 1; --h) {
        int parent = posn - h;
        int i = parent/arr_size;
        int j = parent - i*arr_size;
        if (obj < first_index[i][j].dst) {
            first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = first_index[i][j].dst;
            posn = parent;
        } else {
            first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = obj;
            return;
        }
    }
    first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = obj;
}

void IndirectIndex::right_push(int64_t const &obj, int posn, int height) {
    for (int h = height; h >= 1; --h) {
        int parent = posn - h -1;
        int i = parent/arr_size;
        int j = parent - i*arr_size;
        if (obj < first_index[i][j].dst) {
            first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = first_index[i][j].dst;
            posn = parent;
        } else {
            first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = obj;
            return;
        }
    }
    first_index[posn/arr_size][posn-(posn/arr_size)*arr_size].dst = obj;
}

// first entry of a row of the beap
int IndirectIndex::first(int h) {
    return h*(h+1)/2;
}

// last entry of a row
int IndirectIndex::last(int h) {
    return (h+1)*(h+2)/2-1;
}
