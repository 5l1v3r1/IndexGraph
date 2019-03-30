set(sources
	src/csv.c
	src/stinger_utils.c
	src/timer.c
)
set(headers
	inc/csv.h
	inc/stinger_utils.h
	inc/timer.h
)

publish_headers(headers "${CMAKE_BINARY_DIR}/include/stinger_utils")

include_directories("${CMAKE_BINARY_DIR}/include/stinger_utils")
include_directories("${CMAKE_BINARY_DIR}/include/stinger_core")

# Need to declare an explicit dependency on rapidjson headers, since there's no associated target
set_source_files_properties("${CMAKE_BINARY_DIR}/include/rapidjson/document.h" GENERATED)
add_library(stinger_utils ${sources} ${headers})
target_link_libraries(stinger_utils stinger_core)
