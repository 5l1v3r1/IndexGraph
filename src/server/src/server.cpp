#include <cstdio>
#include <algorithm>
#include <limits>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <netinet/in.h>
#include <unistd.h>
#include <signal.h>
#include <sys/syscall.h>
#include <atomic>

#include "server.h"
#include "stinger_net/stinger_server_state.h"

// @ADD
#include "stinger_core/shared_buffer.h"
#include "stinger_core/vertex_array.h"

extern "C" {
#include "stinger_core/stinger_shared.h"
#include "stinger_core/xmalloc.h"
#include "stinger_utils/stinger_utils.h"
#include "stinger_utils/timer.h"
}

//#define LOG_AT_D
#include "stinger_core/stinger_error.h"

using namespace gt::stinger;

static int start_pipe[2] = {-1, -1};
static char * graph_name = NULL;

pid_t master_tid;
static pthread_t batch_server_tid, alg_server_tid;

static StingerServerState & server_state = StingerServerState::get_server_state();

static void cleanup (void);
extern "C" {
  static void sigterm_cleanup (int);
}

extern int core_num;
extern int mode;
int vertex_num = 5000000;
std::atomic_long rcon_cnt(0);
std::atomic_long rcon_time(0);
double mode_time = 0.0;
double binary_time, beap_time;
int who = 0;

int main(int argc, char *argv[])
{
  /* default global options */
  graph_name = (char*)xmalloc(128*sizeof(char));
  sprintf(graph_name,"/strea---graph");
  int port_streams = 10104;
  int port_algs = 10106;
  int unleash_daemon = 0;

  //int vertex_num = 80000000;
  

  /* parse command line configuration */
  int opt = 0;
  while(-1 != (opt = getopt(argc, argv, "c:v:s:?h"))) {
    switch(opt) {
      case 'c': {
        core_num = atoi(optarg);
      } break;
      case 'v': {
        vertex_num = atoi(optarg);
      } break;
      case 's': {
        mode = atoi(optarg);
      } break;
      // case 'a': {
      //   a = atoi(optarg);
      // } break;
      // case 'b': {
      //   b = atoi(optarg);
      // } break;
      case '?':
      case 'h': {
		  printf("Usage:    %s\n"
			 "   [-c cap number of history files to keep per alg]  \n", argv[0]);
		  printf("Defaults:\n\tport_algs: %d\n\tport_streams: %d\n\tgraph_name:\n", port_algs, port_streams);
		  exit(0);
		} break;

    }
  }

#ifdef __APPLE__
  master_tid = syscall(SYS_thread_selfid);
#else
  master_tid = syscall(SYS_gettid);
#endif

  struct stinger * S = stinger_shared_new_full(&graph_name,NULL);
  size_t graph_sz = S->length + sizeof(struct stinger);

  // @TODO: set vertex array to server state
  VertexArray *va = new_vertex_array(vertex_num);
  // @ADD: new shared-area and init
  SharedBuffer * shared = new SharedBuffer();
  server_state.set_shared(shared);
  server_state.set_varray(va);

  printf("set server_state\n");
  server_state.set_stinger(S);
  server_state.set_stinger_sz(graph_sz);
  server_state.set_port(port_streams, port_algs);

  /* this thread will handle the batch & alg servers */
  /* TODO: bring the thread creation for the alg server to this level */
  pthread_create(&batch_server_tid, NULL, start_batch_server, NULL);
  pthread_create(&alg_server_tid, NULL, start_alg_handling, NULL);

  {
    /* Inform the parent that we're ready for connections. */
    struct sigaction sa;
    sa.sa_flags = 0;
    sigemptyset (&sa.sa_mask);
    sa.sa_handler = sigterm_cleanup;
    /* Ignore the old handlers. */
    sigaction (SIGINT, &sa, NULL);
    sigaction (SIGTERM, &sa, NULL);
    sigaction (SIGHUP, &sa, NULL);
  }

  if(unleash_daemon) {
    int exitcode = EXIT_SUCCESS;
    size_t bytes = write (start_pipe[1], &exitcode, sizeof (exitcode));
    LOG_I_A("Written %lu bytes", bytes);
    close (start_pipe[1]);
    while(1) { sleep(10); }
  } else {
    LOG_I("Press Ctrl-C to shut down the server...");
    while(1) { sleep(10); }
  }

  pthread_join(batch_server_tid, NULL);
  pthread_join(alg_server_tid, NULL);

  return 0;
}

void
cleanup (void)
{
  pid_t tid;
#ifdef __APPLE__
  tid = syscall(SYS_thread_selfid);
#else
  tid = syscall(SYS_gettid);
#endif
  /* Only the main thread executes */
  if (tid == master_tid) {
    LOG_I("Shutting down the batch server..."); fflush(stdout);
    pthread_cancel(batch_server_tid);
    pthread_join(batch_server_tid, NULL);
    LOG_I("done."); fflush(stdout);

    LOG_I("Shutting down the alg server..."); fflush(stdout);
    pthread_cancel(alg_server_tid);
    pthread_join(alg_server_tid, NULL);
    LOG_I("done."); fflush(stdout);

	struct stinger * S = server_state.get_stinger();
	size_t graph_sz = S->length + sizeof(struct stinger);

	/* clean up */
	stinger_shared_free(S, graph_name, graph_sz);
	free(graph_name);

	/* @ADD: free shared_buffer */
    SharedBuffer * shared = server_state.get_shared();
    delete shared;
    // @DONE: free vertex_array
    VertexArray * va = server_state.get_varray();
    free_vertex_array(va);

#ifndef STINGER_USE_TCP
    /* Clean up unix sockets, which were created when the batch/alg server started up */
    char socket_path[128];
    snprintf(socket_path, sizeof(socket_path)-1, "/tmp/stinger.sock.%i", server_state.get_port_algs());
    unlink(socket_path);
    snprintf(socket_path, sizeof(socket_path)-1, "/tmp/stinger.sock.%i", server_state.get_port_streams());
    unlink(socket_path);
#endif
  }
}

void
sigterm_cleanup (int)
{
  cleanup ();
  exit (EXIT_SUCCESS);
}
