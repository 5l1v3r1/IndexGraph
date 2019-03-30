#ifndef SRC_CORE_TRAVERSAL_H_
#define SRC_CORE_TRAVERSAL_H_

#include <omp.h>
#include "stinger_core/vertex_array.h"
#include "stinger_core/indirect_index.h"
#include "stinger_core/shared_buffer.h"
#include "stinger_core/x86_full_empty.h"

#define PARALLEL_ #pragma omp parallel for

// Generic macro for iterating over all edges of a vertex in shared. Edges are writable.
#define SHARED_GENERIC_FORALL_EDGES_OF_VTX_BEGIN(SHARED_,VA_,VTX_,PARALLEL_)\
  do {                                                              \
    PARALLEL_                                          \
    for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {    \
      struct Edge current_edge__ = SHARED_->buffer[i];     \
      if (current_edge__.dst == -1)  break;

#define SHARED_GENERIC_FORALL_EDGES_OF_VTX_END()         \
    } /* end while not last edge */                       \
  } while (0)

// Generic macro for iterating over all edges of a vertex in Index. Edges are writable.
#define INDEX_GENERIC_FORALL_EDGES_OF_VTX_BEGIN(SHARED_,VA_,VTX_,PARALLEL_)\
  do {                                                              \
    IndirectIndex * rect = get_outdata(VA_, VTX_);    \
    PARALLEL_ {                              \
        for (int64_t i = 0; i < rect->index_used; ++i) {    \
            for (int64_t j = 0; j < rect->arr_size; ++j) {    \
                struct Edge current_edge__ = rect->first_index[i][j];     \
                if (current_edge__.dst == -1)  break;

#define INDEX_GENERIC_FORALL_EDGES_OF_VTX_END()         \
            }                                             \
        } /* end while not last edge */                       \
    }                                                      \
  } while (0)

// #define SHARED_FORALL_OUT_EDGES_OF_VTX_BEGIN(SHARED_,VA_,VTX_) \
//     int64_t room_ = VA_->vertices[VTX_].out_buffer_index; \
//     if (room_ != -1) {  \
//         SHARED_GENERIC_FORALL_EDGES_OF_VTX_BEGIN(SHARED_,VA_,VTX_,) \
//     } else {    \
//         INDEX_GENERIC_FORALL_EDGES_OF_VTX_BEGIN(INDEX_,VA,VTX_,) \
//     }
// #define SHARED_FORALL_OUT_EDGES_OF_VTX_END() \
//     if (room_ ！= -1) {  \
//         SHARED_GENERIC_FORALL_EDGES_OF_VTX_END()  \
//     } else {              \
//         INDEX_GENERIC_FORALL_EDGES_OF_VTX_END()   \
//     }

// For all out-edges of vertex in shared
#define SHARED_FORALL_OUT_EDGES_OF_VTX_BEGIN(SHARED_,VA_,VTX_) \
  SHARED_GENERIC_FORALL_EDGES_OF_VTX_BEGIN(SHARED_,VA_,VTX_,)
#define SHARED_FORALL_OUT_EDGES_OF_VTX_END() \
  SHARED_GENERIC_FORALL_EDGES_OF_VTX_END()

// // For all out-edges of vertex in index
#define INDEX_FORALL_OUT_EDGES_OF_VTX_BEGIN(INDEX_,VA,VTX_) \
  INDEX_GENERIC_FORALL_EDGES_OF_VTX_BEGIN(INDEX_,VA,VTX_,)
#define INDEX_FORALL_OUT_EDGES_OF_VTX_END() \
  INDEX_GENERIC_FORALL_EDGES_OF_VTX_END()

#define EDGE_DEST (current_edge__.dst)

/* --------- diameter ---------------------*/
extern int vertex_num;
typedef struct {
  int64_t vertex;
  int64_t cost;
}weighted_vertex_t;

bool comp(weighted_vertex_t a, weighted_vertex_t b) {
  return a.cost < b.cost;
}
typedef bool(*CompareType) (weighted_vertex_t a, weighted_vertex_t b);

std::vector<int64_t> dijkstra(SharedBuffer * shared, VertexArray * va, int64_t source_vertex) {
//int64_t dijkstra(SharedBuffer * shared, VertexArray * va, int64_t source_vertex) {
//  printf("enter into the dijkstra step\n");
  std::vector<int64_t> cost_so_far(vertex_num);
  for (int64_t v = 0; v < vertex_num; v++) {
    cost_so_far[v] = std::numeric_limits<int64_t>::max();
  }
  cost_so_far[source_vertex] = 0;
  std::priority_queue<weighted_vertex_t, std::vector<weighted_vertex_t>, CompareType> frontier (comp);
  weighted_vertex_t source;
  source.vertex = source_vertex;
  source.cost = cost_so_far[source_vertex];
  frontier.push(source);

  while (!frontier.empty()) {
    weighted_vertex_t current = frontier.top();
    frontier.pop();
    // for all out-edges of vertex
    do {
      int64_t outdegree = get_outdegree(va, current.vertex);
      if (outdegree > 0 && outdegree <= THRESHOLD) {
    //    printf("traverse shared \n");
        int64_t room_ = va->vertices[current.vertex].out_buffer_index;
        if (room_ == -1)  {printf("wrong in room\n");}  // exit(-1);}
        for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {
          Edge current_edge__ = shared->buffer[i];
          if (current_edge__.dst == -1)  break;
          int64_t new_cost;
          new_cost = cost_so_far[current.vertex]+1;
          if (new_cost < cost_so_far[EDGE_DEST]) {
            cost_so_far[EDGE_DEST] = new_cost;
            weighted_vertex_t next;
            next.vertex = EDGE_DEST;
            next.cost = new_cost;
            frontier.push(next);
          }
        }
      } else if (outdegree > THRESHOLD) {
    //    printf("traverse index \n");
          IndirectIndex * rect = get_outdata(va, current.vertex);
          if (rect == nullptr)  {printf("wrong in index\n");}  // exit(-1);}
          for (int64_t i = 0; i < rect->index_used; ++i) {
            for (int64_t j = 0; j < rect->arr_size; ++j) {
              Edge current_edge__ = rect->first_index[i][j];
              if (current_edge__.dst == -1)  break;
              int64_t new_cost;
              new_cost = cost_so_far[current.vertex]+1;
              if (new_cost < cost_so_far[EDGE_DEST]) {
                cost_so_far[EDGE_DEST] = new_cost;
                weighted_vertex_t next;
                next.vertex = EDGE_DEST;
                next.cost = new_cost;
                frontier.push(next);
              }
            }
          }
      }
    } while (0);
  }
//  return 0;
  return cost_so_far;
}

int64_t pseudo_diameter(SharedBuffer * shared, VertexArray * va, int64_t source, int64_t dist) {
  dist = 0;
  int64_t target = source;

  printf("vertex num is %ld\n", vertex_num);
  std::vector<int64_t> paths(vertex_num);
  while (1) {
    int64_t new_source = target;
    paths = dijkstra(shared, va, new_source);
    int64_t max = std::numeric_limits<int64_t>::min();
    int64_t max_index = 0;
    for (int64_t i = 0; i < vertex_num; i++) {
      if (paths[i] > max && paths[i] != std::numeric_limits<int64_t>::max()) {
        max = paths[i];
        max_index = i;
      }
    }
    if (max > dist) {
      target = max_index;
      dist = max;
    } else {
      break;
    }
  }
  return dist;
}

/* --------- pagerank ---------------------*/
#define EPSILON_DEFAULT 1e-8
#define DAMPINGFACTOR_DEFAULT 0.85
#define MAXITER_DEFAULT 20

inline double * set_tmp_pr(double * tmp_pr_in, int64_t nv) {
  double * tmp_pr = NULL;
  if (tmp_pr_in) {
    tmp_pr = tmp_pr_in;
  } else {
    tmp_pr = (double*)xmalloc(sizeof(double)*nv);
  }
  return tmp_pr;
}

inline void unset_tmp_pr(double * tmp_pr, double * tmp_pr_in) {
  if (!tmp_pr_in)
    free(tmp_pr);
}

// only works on undirected graphs
int64_t page_rank(SharedBuffer * shared, VertexArray * va, double *pr, double *tmp_pr_in, 
double epsilon, double dampingfactor, int64_t maxiter) {
  double *tmp_pr = set_tmp_pr(tmp_pr_in, vertex_num);
  int64_t iter = maxiter;
  double delta = 1;
  int64_t iter_count = 0;

  while (delta > epsilon && iter > 0) {
    iter_count++;
    double pr_constant = 0.0;

    OMP("omp parallel for reduction(+:pr_constant)")
    for (uint64_t v = 0; v < vertex_num; v++) {
      tmp_pr[v] = 0;
      if (get_outdegree(va, v) == 0) {
        pr_constant += pr[v];
      } else {
        // for all in-edges of vertex
        do {
          int64_t indegree = get_indegree(va, v);
          if (indegree > 0 && indegree <= THRESHOLD) {
            int64_t room_ = va->vertices[v].in_buffer_index;
            if (room_ == -1)  {printf("wrong in room\n");}  // exit(-1);}
            for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {
              Edge current_edge__ = shared->buffer[i];
              if (current_edge__.dst == -1)  break;
              int64_t outdegree = get_outdegree(va, EDGE_DEST);
              tmp_pr[v] += (((double)pr[EDGE_DEST]) / ((double)(outdegree? outdegree: vertex_num-1)));
            }
          } else if (indegree > THRESHOLD) {
            IndirectIndex * rect = get_indata(va, v);
            if (rect == nullptr)  {printf("wrong in index\n");}  // exit(-1);}
            for (int64_t i = 0; i < rect->index_used; ++i) {
              for (int64_t j = 0; j < rect->arr_size; ++j) {
                Edge current_edge__ = rect->first_index[i][j];
                if (current_edge__.dst == -1)  break;
                int64_t outdegree = get_outdegree(va, EDGE_DEST);
                tmp_pr[v] += (((double)pr[EDGE_DEST]) / ((double)(outdegree? outdegree: vertex_num-1)));
              }
            }
          }
        } while (0);
        // traverse over
      }
    }

    OMP("omp parallel for")
    for (uint64_t v = 0; v < vertex_num; v++) {
      tmp_pr[v] = (tmp_pr[v] + pr_constant/(double)vertex_num) * dampingfactor + (((double)(1-dampingfactor))/((double)vertex_num));
    }

    delta = 0;
    OMP("omp parallel for reduction(+:delta)")
    for (uint64_t v = 0; v < vertex_num; v++) {
      double mydelta = tmp_pr[v] - pr[v];
      if (mydelta < 0)
        mydelta = -mydelta;
      delta += mydelta;
    }

    OMP("omp parallel for")
    for (uint64_t v = 0; v < vertex_num; v++) {
      pr[v] = tmp_pr[v];
    }
    iter--;
  }
  LOG_I_A("PageRank iteration count : %ld", iter_count);
  unset_tmp_pr(tmp_pr,tmp_pr_in);
}

/*---------------------WCC------------------------------*/
#define EDGE_SOURCE source__

int64_t parallel_shiloach_vishkin_components_of_type(SharedBuffer * shared, VertexArray * va, 
int64_t *component_map, int64_t type) {
  int64_t nv = vertex_num;

  OMP ("omp parallel for")
  for (uint64_t i = 0; i < nv; i++) {
    component_map[i] = i;
  }

  while (1) {
    int changed = 0;
    // for all edges of a given type, in parallel
    OMP("omp parallel for")
    for (uint64_t v = 0; v < nv; v++) {
      do {
        int64_t outdegree = get_outdegree(va, v);
        if (outdegree > 0 && outdegree <= THRESHOLD) {
          int64_t room_ = va->vertices[v].out_buffer_index;
          if (room_ == -1)  {printf("wrong in room\n");} // exit(-1);}
          int64_t source__ = v; 

        //  OMP("omp parallel for")      
          for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {
            Edge current_edge__ = shared->buffer[i];
            if (current_edge__.dst == -1)  continue;
            int64_t c_src = component_map[EDGE_SOURCE]; 
            int64_t c_dest = component_map[EDGE_DEST];
            if (c_dest < c_src) {
              component_map[v] = c_dest;
              changed++;
            }
            if (c_src < c_dest) {
              component_map[EDGE_DEST] = c_src;
              changed++;
            }
          }
        } else if (outdegree > THRESHOLD) {
          IndirectIndex * rect = get_outdata(va, v);
          if (rect == nullptr)  {printf("wrong in index\n");}  // exit(-1);}
          int64_t source__ = v; 

        //  OMP("omp parallel for")
          for (int64_t i = 0; i < rect->index_used; ++i) {
            for (int64_t j = 0; j < rect->arr_size; ++j) {
              Edge current_edge__ = rect->first_index[i][j];
              if (current_edge__.dst == -1)  break;
              int64_t c_src = component_map[EDGE_SOURCE];
              int64_t c_dest = component_map[EDGE_DEST];
              if (c_dest < c_src) {
                component_map[v] = c_dest;
                changed++;
              }
              if (c_src < c_dest) {
                component_map[EDGE_DEST] = c_src;
                changed++;
              }
            }
          }
        }
      } while (0);
    }

    if (!changed)  break;

    OMP ("omp parallel for")
    for (uint64_t i = 0; i < nv; i++) {
      while (component_map[i] != component_map[component_map[i]]) {
        component_map[i] = component_map[component_map[i]];
      }
    }
  }
  return 0;
}

int64_t compute_component_sizes (SharedBuffer * shared, VertexArray * va,
int64_t * component_map, int64_t * component_size) {
  int64_t nv = vertex_num;

  OMP ("omp parallel for")
  for (uint64_t i = 0; i < nv; i++) {
    component_size[i] = 0;
  }

  OMP ("omp parallel for")
  for (uint64_t i = 0; i < nv; i++) {
    int64_t c_num = component_map[i];
    stinger_int64_fetch_add(&component_size[c_num], 1);
  }
  return 0;
}

/*----------------------------bfs----------------------------*/
int64_t parallel_breadth_first_search (SharedBuffer * shared, VertexArray * va,
                            int64_t source, int64_t * marks, int64_t * level)
{
  int64_t nv = vertex_num;
  std::vector<int64_t> Qhead(nv);
  std::vector<int64_t> queue(nv);
  for (int64_t i = 0; i < nv; i++) {
    level[i] = -1;
    marks[i] = 0;
  }

  int64_t nQ, Qnext, Qstart, Qend;
  /* initialize */
  queue[0] = source;
  level[source] = 0;
  marks[source] = 1;
  Qnext = 1;    /* next open slot in the queue */
  nQ = 1;         /* level we are currently processing */
  Qhead[0] = 0;    /* beginning of the current frontier */
  Qhead[1] = 1;    /* end of the current frontier */

  Qstart = Qhead[nQ-1];
  Qend = Qhead[nQ];

  while (Qstart != Qend) {
    OMP ("omp parallel for")
    for (int64_t j = Qstart; j < Qend; j++) {
      //For all out-edges of vertex
      do {
        int64_t outdegree = get_outdegree(va, queue[j]);
        if (outdegree > 0 && outdegree <= THRESHOLD) {
          int64_t room_ = va->vertices[queue[j]].out_buffer_index;
          if (room_ == -1)  {printf("wrong in room\n");}  // exit(-1);}
          for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {
            Edge current_edge__ = shared->buffer[i];
            if (current_edge__.dst == -1)  break;
            int64_t d = level[EDGE_DEST];
            if (d < 0) {
              if (stinger_int64_fetch_add (&marks[EDGE_DEST], 1) == 0) {
                level[EDGE_DEST] = nQ;
                int64_t mine = stinger_int64_fetch_add(&Qnext, 1);
                queue[mine] = EDGE_DEST;
              }
            }
          }
        } else if (outdegree > THRESHOLD) {
          IndirectIndex * rect = get_outdata(va, queue[j]);
          if (rect == nullptr)  {printf("wrong in index\n");}  // exit(-1);}
          for (int64_t i = 0; i < rect->index_used; ++i) {
            for (int64_t j = 0; j < rect->arr_size; ++j) {
              Edge current_edge__ = rect->first_index[i][j];
              if (current_edge__.dst == -1)  break;
              int64_t d = level[EDGE_DEST];
              if (d < 0) {
                if (stinger_int64_fetch_add (&marks[EDGE_DEST], 1) == 0) {
                  level[EDGE_DEST] = nQ;
                  int64_t mine = stinger_int64_fetch_add(&Qnext, 1);
                  queue[mine] = EDGE_DEST;
                }
              }
            }
          }
        }
      } while (0);
    }

    Qstart = Qhead[nQ-1];
    Qend = Qnext;
    Qhead[nQ++] = Qend;
  }
  return nQ;
}

// template ___  all in-edges
      // do {
      //     int64_t indegree = get_indegree(va, v);
      //     if (indegree > 0 && indegree <= THRESHOLD)
      //       int64_t room_ = va->vertices[v].in_buffer_index;
      //       if (room_ == -1)  {printf("wrong in room\n"); exit(-1);}
      //       for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {
      //         Edge current_edge__ = shared->buffer[i];
      //         if (current_edge__.dst == -1)  break;

      //       }
      //     } else if (indegree > THRESHOLD) {
      //       IndirectIndex * rect = get_indata(va, v);
      //       if (rect == nullptr)  {printf("wrong in index\n"); exit(-1);}
      //       for (int64_t i = 0; i < rect->index_used; ++i) {
      //         for (int64_t j = 0; j < rect->arr_size; ++j) {
      //           Edge current_edge__ = rect->first_index[i][j];
      //           if (current_edge__.dst == -1)  break;

      //         }
      //       }
      //     }
      //   } while (0);

// template ___ all out-edges
      // do {
      //     int64_t outdegree = get_outdegree(va, v);
      //     if (outdegree > 0 && outdegree <= THRESHOLD)
      //       int64_t room_ = va->vertices[v].out_buffer_index;
      //       if (room_ == -1)  {printf("wrong in room\n"); exit(-1);}
      //       for (int64_t i = room_ * THRESHOLD; i < (room_ + 1) * THRESHOLD; ++i) {
      //         Edge current_edge__ = shared->buffer[i];
      //         if (current_edge__.dst == -1)  break;

      //       }
      //     } else if (outdegree > THRESHOLD) {
      //       IndirectIndex * rect = get_outdata(va, v);
      //       if (rect == nullptr)  {printf("wrong in index\n"); exit(-1);}
      //       for (int64_t i = 0; i < rect->index_used; ++i) {
      //         for (int64_t j = 0; j < rect->arr_size; ++j) {
      //           Edge current_edge__ = rect->first_index[i][j];
      //           if (current_edge__.dst == -1)  break;

      //         }
      //       }
      //     }
      //   } while (0);
#endif
