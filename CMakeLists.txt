cmake_minimum_required(VERSION 3.5.0)
project(fastaconvtr VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()


add_compile_options("-DLOG_USE_COLOR")
IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)
add_compile_options("-DDEBUG")

file (GLOB_RECURSE SOURCE_FILE CONFIGURE_DEPENDS  "sources/*.c*")
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(${PROJECT_NAME} ${SOURCE_FILE})


target_link_libraries(${PROJECT_NAME} -lz)
target_link_libraries(${PROJECT_NAME} -lm)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)
