#include <stdint.h>
#include <stdio.h>
#include "stinger_core/vertex_array.h"
#include "stinger_core/indirect_index.h"
#include "stinger_core/shared_buffer.h"
#include "stinger_core/traversal.h"
#include "stinger_net/proto/stinger-batch.pb.h"
//#include "stinger_net/proto/stinger-connect.pb.h"
//#include "stinger_net/proto/stinger-alg.pb.h"
#include "stinger_net/stinger_server_state.h"

using namespace gt::stinger;

int main(int argc, char const *argv[])
{
    printf("start to exec diameter\n");
    StingerServerState & server_state = StingerServerState::get_server_state();
    SharedBuffer * s = server_state.get_shared();
    VertexArray * va = server_state.get_varray();
    if (s == NULL || va == nullptr) {
        printf("its null\n");
        exit(-1);
    }
    int64_t source_vertex = 0;
    int64_t dist = 0;
    pseudo_diameter(server_state.get_shared(), server_state.get_varray(), source_vertex, dist);
    return 0;
}
