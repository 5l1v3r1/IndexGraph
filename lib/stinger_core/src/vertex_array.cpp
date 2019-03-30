/*
 * Copyright 2018-now Zhihui Zhao
 * @Author: Zhihui Zhao
 * @Last Modified by: Yilong Liu
 */

#include <cstdio>
#include <cstdlib>

#include "stinger_core/vertex_array.h"
#include "stinger_core/indirect_index.h"

// vertex 相关的方法
Vertex* new_vertices(int64_t num_vertices) {
    printf("new vertex array of size %ld\n", num_vertices);

    Vertex* vertices = new Vertex[num_vertices];

    for (int i = 0; i < num_vertices; ++i) {
        vertices[i].in_degree = 0;
        vertices[i].out_degree = 0;
        vertices[i].out_buffer_index = -1;
        vertices[i].in_buffer_index = -1;

        vertices[i].in_data = nullptr;
        vertices[i].out_data = nullptr;
    }

    return vertices;
}

void free_vertices(Vertex* vertices, int num) {
    int out = 0, in = 0;
    for (int i = 0; i < num; ++i) {
        if (vertices[i].in_data != nullptr) {
            IndirectIndex * rect = (IndirectIndex *)vertices[i].in_data;
            //++in;
            __sync_fetch_and_add(&in, 1);
            delete rect;
        }
        if (vertices[i].out_data != nullptr) {
            IndirectIndex * rect = (IndirectIndex *)vertices[i].out_data;
            //++out;
            __sync_fetch_and_add(&out, 1);
        //    printf("------free a out_Indirect index\n");
            delete rect;
        }
    }

    delete vertices;
  //  free(vertices);
    vertices = nullptr;
}

VertexArray *new_vertex_array(int64_t num_vertices) {
    VertexArray *vertex_array = new VertexArray();
    vertex_array->vertices = new_vertices(num_vertices);
    vertex_array->num_vertices = num_vertices;

    return vertex_array;
}

void free_vertex_array(VertexArray *vertex_array) {
    LOG_D("Entering Vertex_Array destructor");
    free_vertices(vertex_array->vertices, vertex_array->num_vertices);
    delete vertex_array;
    vertex_array = nullptr;
}

