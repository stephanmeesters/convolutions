cmake_minimum_required(VERSION 2.8)

PROJECT(Convolutions)

LIST(APPEND
CMAKE_MODULE_PATH "${CMAKE_SOURCE_PREFIX}")
set(CMAKE_MODULE_PATH
${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} lib)


  include_directories(/usr/local/Cellar/niftilib/2.0.0/include/)
  #link_directories(/usr/local/Cellar/niftilib/2.0.0/lib/libniftiio.a)




find_package(zlib)
include_directories(${ZLIB_INCLUDE_DIR})

file(GLOB CXX_FILES *.cxx)
file(GLOB H_FILES *.h)

if (APPLE)
  add_executable(Convolutions MACOSX_BUNDLE ${CXX_FILES} ${H_FILES})
else()
  add_executable(Convolutions ${CXX_FILES} ${H_FILES})
endif()

target_link_libraries(Convolutions /usr/local/Cellar/niftilib/2.0.0/lib/libniftiio.a /usr/local/Cellar/niftilib/2.0.0/lib/libznz.a /usr/local/Cellar/niftilib/2.0.0/lib/libnifticdf.a /usr/local/Cellar/niftilib/2.0.0/lib/libfslio.a ${ZLIB_LIBRARY})




























