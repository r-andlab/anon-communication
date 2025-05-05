cmake_minimum_required(VERSION 3.30)

# Use CMake's new Boost package policy (instead of deprecated FindBoost)
cmake_policy(SET CMP0167 NEW)

project(anonymous_communication LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
  -march=native # needed for AES intrinsics
  -Wall
  -Wno-ignored-attributes
  -ggdb
  -O3
)

set(SOURCE_DIR src)

# Dependencies

find_package(Boost REQUIRED system context chrono thread)
find_library(BSD_LIB NAMES bsd)

add_subdirectory(vendor/rpclib)

# Common Sources

set(COMMON_SRCS
  ${SOURCE_DIR}/avl.cpp
  ${SOURCE_DIR}/bst.cpp
  ${SOURCE_DIR}/cdpf.cpp
  ${SOURCE_DIR}/rdpf.cpp
  ${SOURCE_DIR}/cell.cpp
  ${SOURCE_DIR}/duoram.cpp
  ${SOURCE_DIR}/heap.cpp
  ${SOURCE_DIR}/heapsampler.cpp
  ${SOURCE_DIR}/remise.cpp
  ${SOURCE_DIR}/mpcio.cpp
  ${SOURCE_DIR}/mpcops.cpp
)

add_library(anoncomm_common OBJECT ${COMMON_SRCS})

target_link_libraries(anoncomm_common
  ${BSD_LIB}
  Boost::system
  Boost::context
  Boost::chrono
  Boost::thread
  pthread
)

target_include_directories(anoncomm_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
)

# Client

set(CLIENT_SRCS
  ${SOURCE_DIR}/client/main.cpp
)

add_executable(anoncomm_client ${CLIENT_SRCS})

target_link_libraries(anoncomm_client
  anoncomm_common
  ${BSD_LIB}
  Boost::system
  Boost::context
  Boost::chrono
  Boost::thread
  pthread
  rpc
)

target_include_directories(anoncomm_client PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
  ${RPCLIB_INCLUDE_DIR}
)

# Server

option(VERBOSE_COMMS "Enable verbose communication logging" OFF)
if(VERBOSE_COMMS)
  add_compile_definitions(VERBOSE_COMMS)
endif()

set(SERVER_SRCS
  ${SOURCE_DIR}/server/prac.cpp
  ${SOURCE_DIR}/server/preproc.cpp
  ${SOURCE_DIR}/server/online.cpp
  ${SOURCE_DIR}/server/anoncomm.cpp
)

add_executable(anoncomm_server ${SERVER_SRCS})

target_link_libraries(anoncomm_server PRIVATE
  anoncomm_common
  ${BSD_LIB}
  Boost::system
  Boost::context
  Boost::chrono
  Boost::thread
  pthread
  rpc
)

target_include_directories(anoncomm_server PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
  ${RPCLIB_INCLUDE_DIR}
)
