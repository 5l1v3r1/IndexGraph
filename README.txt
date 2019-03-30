IndexGraph is a package, modified from Stinger, with indexing mechanism that organizes the adjacency lists of vertices in a dynamic graph with indexed structures, to improve the effciency of dynamic graph processing systems that work in the memory of a single COTS server. 


---------------------------------------  Prerequisites for building  ------------------------------
Linux (Ubuntu 16.04.3 LTS, 4.4.0-87-generic kernel, x86_64 for our case);
basic build tools (e.g., stdlibc, gcc, etc);
g++ (5.4.0 for our case);

-------------------------------------------  Explanations  ----------------------------------------
1) Where is the source code?
	see {Project root}/

2) How to build it?
    IndexGraph is built using CMake, and we provide a script to build in a simple way:
	{Project root}$ make
    More specifically, if you want to build step by step, first create a build directory:
	{Project root}$ mkdir build
    Then call CMake to automatically configure the build and to create a Makefile:
	{Project root}$ ./build/cmake ..
    Finally, call make to build all libraries and executable targets (or call make and the name of an executable or library to build):
	{Project root}$ ./build/make -j3
    Note: the -j flag is a multi-threaded build. Typically you should match the argument to the number of cores on your system.
    All binary targets will be built and placed in build/bin. They are named according to the folder from which they were built (so src/bin/server produces build/bin/stinger_server, src/bin/clients/tools/json_rpc_server produces build/bin/stinger_json_rpc_server, etc.).

3) How to run it?
    To run an example using the server and two terminals:

	term1: {Project root}$ ./build/bin/stream_graph -v 1432693 -s 0
	term2: {Project root}$ cat /space/pokec.txt | ./build/bin/stinger_csv_stream ./template.txt -x 100000

    This will start a stream of pokec edges over 1432,693 vertices in batches of 100,000 edges.
    -v input the the number of vertices of the source file.
    -s determines the searching and inserting method. 0/1/2: BU / Beap / Hybrid scheme.
     '/space/pokec.txt' is the location of the input file. You can replace it as needed.
    './template.txt' is the sample format of the input file. You can edit it according to yours.
    -x determines the size of a batch.

4) Server Configuration
    CMake Parameters: The IndexGraph is configured via several methods when using the IndexGraph Server. There are several CMake Parameters that can be set when configuring the build via CMake. These are listed below.

	STINGER_DEFAULT_VERTICES -> Specifies the number of vertices
	STINGER_DEFAULT_NEB_FACTOR -> Specifies the Number of Edge Blocks allocated for each vertex per edge type
	STINGER_DEFAULT_NUMETYPES -> Specifies the Number of Edge Types
	STINGER_DEFAULT_NUMVTYPES -> Specifies the Number of Vertex Types
	STINGER_EDGEBLOCKSIZE -> Specifies the minimum number of edges of each type that can be attached to a vertex. This should be a multiple of 2
	STINGER_NAME_STR_MAX -> Specifies the maximum length of strings when specifying edge type names and vertex type names
	STINGER_NAME_USE_SQLITE -> (Alpha - Still under development) Replaces the static names structure with a dynamic structure that uses SQLITE
	STINGER_USE_TCP -> Uses TCP instead of Unix Sockets to connect to the IndexGraph server

    Note: IndexGraph, by default, uses one edge type name and vertex type name for the "None" type. This can be disabled using the IndexGraph server config file described in a later section.

    Example: Suppose there is a graph you want to store with a power-law distribution of edges, and the average degree of each vertex is 3. There are 3 named edge types and 4 name vertex types. You wish to store up to 2 million vertices. A good configuration would be

	STINGER_DEFAULT_VERTICES = 1 << 21
	STINGER_EDGEBLOCKSIZE = 4            # On average one edge block is sufficient for most vertices
	STINGER_DEFAULT_NEB_FACTOR = 2       # Depending on the skewness 3 may also be necessary.  This will create 2 * nv edge blocks for each edge types
	STINGER_DEFAULT_NUMETYPES = 4        # There are 3 named types and the default "None" type
	STINGER_DEFAULT_NUMVTYPES = 5        # There are 4 named types and the default "None" type

    Environment Variables: By default, IndexGraph will attempt to allocate 1/2 of the available memory on the system. This can be overwritten using an environment variable, STINGER_MAX_MEMSIZE. For example: 

	STINGER_MAX_MEMSIZE=4g ./bin/stinger_server

    If the IndexGraph parameters provided would create a IndexGraph that is larger than STINGER_MAX_MEMSIZE, then the IndexGraph is shrunk to fit the specified memory size. This is done by repeatedly reducing the number of vertices by 25% until the IndexGraph fits in the specified memory size.

    Config File: The IndexGraph server accepts a config file format in the libconfig style. A config filename can be passed as an argument to the stinger_server with the -C flag. An example config file is provided in stinger.cfg in the root directory of the IndexGraph repository. All of the cmake configurations regarding the size of the IndexGraph can be overridden via the config file except for the STINGER_EDGEBLOCKSIZE. The config file values will override any CMake values. If the STINGER_MAX_MEMSIZE environment variable is set, it will be considered an upper limit on the memory size.
    The following values are allowed in stinger.cfg:
	num_vertices -> A long integer representing the number of vertices. Must have the L at the end of the integer
	edges_per_type -> A long integer representing the number of edges for each vertex type. Must have the L at the end of the integer
	num_edge_types -> Specifies number of edge types.
	num_vertex_types -> Specifies number of vertex types.
	max_memsize -> Specifies maximum memory for IndexGraph. Cannot be bigger than the STINGER_MAX_MEMSIZE environment variable. Uses the same format as the environment variable (e.g. "4G", "1T", "512M")
	edge_type_names -> A list of edge type names to map at IndexGraph startup. If this list is specified the "None" type will not be mapped.
	vertex_type_names -> A list of vertex type names to map at IndexGraph startup. If this list is specified the "None" type will not be mapped.
	map_none_etype -> If set to false, the "None" edge type will not be mapped at startup. Likewise, if true, the "None" edge type will be mapped at startup.
	map_none_vtype -> If set to false, the "None" vertex type will not be mapped at startup. Likewise, if true, the "None" vertex type will be mapped at startup.
	no_resize -> If set to true, the the IndexGraph will not try to resize and fit in the specified memory size. Instead it will throw an error and exit.

5) Example dataset
    We will use pokec as example dataset which will be put in {Project root}/space.
    pokec.txt is a complete graph.
    pokec_base.txt is a graph with 50% edges.
    pokec_update.txt is a graph with the remaining 50% edges.

---------------------------------------------------------------------
The end. enjoy!

