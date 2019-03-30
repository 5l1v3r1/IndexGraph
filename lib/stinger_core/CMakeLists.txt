set(sources
	src/vertex_array.cpp
	src/indirect_index.cpp
	src/shared_buffer.cpp
	src/core_util.c
	src/stinger.c
    src/stinger_deprecated.c
	src/stinger_names.c
	src/stinger_names_sqlite.c
	src/stinger_physmap.c
	src/stinger_shared.c
	src/stinger_vertex.c
	src/xmalloc.c
	src/x86_full_empty.c
)
set(headers
	inc/vertex_array.h
	inc/indirect_index.h
	inc/bitmap.hpp
	inc/defs.h
	inc/shared_buffer.h
	inc/traversal.h
	inc/core_util.h
	inc/stinger.h
	inc/stinger_atomics.h
    inc/stinger_deprecated.h
	inc/stinger_error.h
	inc/stinger_internal.h
	inc/stinger_physmap.h
	inc/stinger_return.h
	inc/stinger_shared.h
	inc/stinger_traversal.h
	inc/stinger_vertex.h
	inc/x86_full_empty.h
	inc/xmalloc.h
	inc/stinger_batch_insert.h
)

set(config
	stinger_names.h
	stinger_defs.h
)

publish_headers(headers "${CMAKE_BINARY_DIR}/include/stinger_core")

# find specific components from boost
find_package(Boost COMPONENTS system thread REQUIRED)

# include headers for boost
if(Boost_THREAD_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories("${CMAKE_BINARY_DIR}/include/stinger_core")
include_directories("${CMAKE_BINARY_DIR}/include/stinger_utils")

set_source_files_properties(${config} PROPERTIES GENERATED TRUE)
add_library(stinger_core ${sources} ${headers} ${config})

if(OPENMP_FOUND)
  target_compile_definitions(stinger_core PUBLIC _GLIBCXX_PARALLEL)
endif()

# link library of boost component to target
if(Boost_THREAD_FOUND)
  target_link_libraries(stinger_core LINK_PUBLIC ${Boost_LIBRARIES})
endif()
