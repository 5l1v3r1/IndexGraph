cmake_minimum_required (VERSION 2.8.11)
project (STINGER)

include (ExternalProject)
include (CTest)

set(STINGER_VERSION_MAJOR 2017)
set(STINGER_VERSION_MINOR 02)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
#set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "--no-as-needed -ldl")
set(STINGER_NAME_USE_SQLITE FALSE CACHE BOOL "Use sqlite for virtual to physmap")
if(STINGER_NAME_USE_SQLITE)
  add_definitions(-DNAME_USE_SQLITE)
endif()

set(STINGER_USE_TCP FALSE CACHE BOOL "Use TCP instead of Unix Sockets")
if(STINGER_USE_TCP)
  add_definitions(-DSTINGER_USE_TCP)
endif()

find_package( OpenMP )
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

link_libraries(m)

option(BUILD_SHARED_LIBS "Build libraries as SHARED" ON)

IF(APPLE)
  # Apple doesn't support librt, using CoreServices instead
  link_libraries(/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices)
  set(CMAKE_MACOSX_RPATH ON)
ELSEIF(WIN32)
  # Not sure what to do here as far as librt...?
  link_libraries(rt)
ELSE()
  link_libraries(rt)
  link_libraries(dl)
ENDIF()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -Wno-discarded-qualifiers -Wno-int-conversion -Wno-implicit-function-declaration -Wno-unused-result")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -std=c++11")

# Doxygen
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if (BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
  endif()
endif()

if(APPLE)
  set(SHARED_LIB_EXT .dylib)
else()
  set(SHARED_LIB_EXT .so)
endif()

#protobuf
find_package(Protobuf)
# Look for installed version
if (${PROTOBUF_FOUND})
  add_custom_target("protobuf")

  add_library(_protobuf_library STATIC IMPORTED)
  add_dependencies(_protobuf_library protobuf)
  set_target_properties(_protobuf_library PROPERTIES IMPORTED_LOCATION ${PROTOBUF_LIBRARY})

  add_library(_protobuf_lite_library STATIC IMPORTED)
  add_dependencies(_protobuf_lite_library protobuf)
  set_target_properties(_protobuf_lite_library PROPERTIES IMPORTED_LOCATION ${PROTOBUF_LITE_LIBRARY})

  include_directories("${PROTOBUF_INCLUDE_DIR}")
# Download and build it ourselves
else ()
  ExternalProject_Add(protobuf
    PREFIX              ${CMAKE_BINARY_DIR}/external
    URL                 file://${CMAKE_SOURCE_DIR}/external/protobuf-cpp-3.1.0.tar.gz
    UPDATE_COMMAND      ""
    PATCH_COMMAND       ""
    CONFIGURE_COMMAND   CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_BINARY_DIR}/external/src/protobuf/configure
      --prefix=${CMAKE_BINARY_DIR}
      --enable-shared
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
  )

  set(PROTOBUF_PROTOC_EXECUTABLE "./protoc")

  add_library(_protobuf_library STATIC IMPORTED)
  add_dependencies(_protobuf_library protobuf)
  set_target_properties(_protobuf_library PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libprotobuf${SHARED_LIB_EXT})

  add_library(_protobuf_lite_library STATIC IMPORTED)
  add_dependencies(_protobuf_lite_library protobuf)
  set_target_properties(_protobuf_lite_library PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libprotobuf-lite${SHARED_LIB_EXT})

endif()

#libconfig
ExternalProject_Add(libconfig
  PREFIX              ${CMAKE_BINARY_DIR}/external
  URL                 file://${CMAKE_SOURCE_DIR}/external/libconfig-1.5.tar.gz
  URL_MD5             a939c4990d74e6fc1ee62be05716f633
  UPDATE_COMMAND      ""
  PATCH_COMMAND       ""
  CONFIGURE_COMMAND   CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_BINARY_DIR}/external/src/libconfig/configure
    --prefix=${CMAKE_BINARY_DIR}
    --enable-shared
    --disable-examples
  BUILD_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/libconfig/lib/libconfig.h ${CMAKE_BINARY_DIR}/external/src/libconfig-build/lib/ && $(MAKE)
  INSTALL_COMMAND $(MAKE) install
)

add_library(_libconfig_library STATIC IMPORTED)
add_dependencies(_libconfig_library libconfig)
set_target_properties(_libconfig_library PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libconfig${SHARED_LIB_EXT})

add_library(_libconfig++_library STATIC IMPORTED)
add_dependencies(_libconfig++_library libconfig)
set_target_properties(_libconfig++_library PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libconfig++${SHARED_LIB_EXT})

# Memory size configuration bits for stinger_core
set(STINGER_DEFAULT_VERTICES "(1L<<24)" CACHE STRING "Default number of vertices")
set(STINGER_DEFAULT_NUMETYPES "5" CACHE STRING "Default number of edge types")
set(STINGER_DEFAULT_NUMVTYPES "128" CACHE STRING "Default number of vertex types")
set(STINGER_DEFAULT_NEB_FACTOR "4" CACHE STRING "Default number of edge blocks per vertex")
set(STINGER_EDGEBLOCKSIZE "14" CACHE STRING "Number of edges per edge block")
set(STINGER_NAME_STR_MAX "255" CACHE STRING "Max string length in physmap")

configure_file(${CMAKE_SOURCE_DIR}/lib/stinger_core/inc/stinger_defs.h.in ${CMAKE_BINARY_DIR}/include/stinger_core/stinger_defs.h @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/lib/stinger_core/inc/stinger_names.h.in ${CMAKE_BINARY_DIR}/include/stinger_core/stinger_names.h @ONLY)

# Define a custom function for copying headers to the build directory
function(publish_headers header_list destination)
  set(published_headers "")
  foreach(header ${${header_list}})
    get_filename_component(name ${header} NAME)
    set(output "${destination}/${name}")
    list(APPEND published_headers ${output})
    add_custom_command(
      OUTPUT ${output}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${destination}
      COMMAND ${CMAKE_COMMAND} -E copy ${header} ${destination}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS ${header}
    )
  endforeach()
  # Overwrite the list of headers in the caller, so targets depend on the published version
  set(${header_list} ${published_headers} PARENT_SCOPE)
endfunction()

include_directories("${CMAKE_BINARY_DIR}/include")

add_subdirectory(lib)
add_subdirectory(src)
