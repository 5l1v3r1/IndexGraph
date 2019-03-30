set(_diameter_sources
	diameter/main.cpp
)

include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR}/stinger_net)
include_directories(${CMAKE_BINARY_DIR}/stinger_core)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/diameter)
link_libraries(_libconfig++_library)
add_executable(stream_diameter ${_diameter_sources})
target_link_libraries(stream_diameter stinger_net stinger_core)

#============================================================