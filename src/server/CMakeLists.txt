set(sources
	src/batch_server.cpp
	src/server.cpp
	src/alg_handling.cpp
	src/batch.cpp
)

set(headers
	inc/server.h
)

include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR}/stinger_net)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
link_libraries(_libconfig++_library)
add_executable(stream_graph ${sources} ${headers})
target_link_libraries(stream_graph stinger_net)
