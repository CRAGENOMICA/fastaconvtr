cmake_minimum_required(VERSION 3.5.0)
project(fastaconvtr VERSION 1.0.0 LANGUAGES C)

#include(CTest)
#enable_testing()



# Get the current date and time
string(TIMESTAMP CURRENT_DATE "%Y%m%d")
string(TIMESTAMP CURRENT_TIME "%H%M%S")


# Concatenate date and time to create a build number
set(BUILD_NUMBER "${CURRENT_DATE}${CURRENT_TIME}")

# Define the version number
set(VERSION_NUMBER "${fastaconvtr_VERSION}")

# Define preprocessor macros for date, time, and build number
add_definitions(-DCURRENT_DATE=\"${CURRENT_DATE}\")
add_definitions(-DCURRENT_TIME=\"${CURRENT_TIME}\")
add_definitions(-DBUILD_NUMBER=\"${BUILD_NUMBER}\")
add_definitions(-DVERSION_NUMBER=\"${VERSION_NUMBER}\")

# Write version number to a text file
file(WRITE "${CMAKE_BINARY_DIR}/version.txt" "${VERSION_NUMBER}")

add_compile_options("-DLOG_USE_COLOR")
IF(CMAKE_BUILD_TYPE MATCHES Debug)

# Set compiler options for warnings and AddressSanitizer
# if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     # add_compile_options(-Wall -Wextra -pedantic -fsanitize=address -g)
#     add_compile_options(-Wall -Wextra -pedantic -fsanitize=address -g -O0)
#     add_link_options(-fsanitize=address)
# endif()

add_compile_options("-DDEBUG")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

file (GLOB_RECURSE SOURCE_FILE CONFIGURE_DEPENDS  "sources/*.c*")
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# log library
add_library(logging "lib/log.c")
add_library(tfasta "lib/tfasta.c" "lib/files_util.c")


# Assuming HTSlib was installed in /usr/local
find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS /usr/local/include)
find_library(HTSLIB_LIBRARY NAMES hts PATHS /usr/local/lib)

## find htslib on windows NEED TO BE TESTED
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "C:/Program Files/htslib/include")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "C:/Program Files/htslib/lib")
endif()
## find htslib on macos NEED TO BE TESTED
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "/usr/local/Cellar/htslib")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "/usr/local/Cellar/htslib")
endif()
# /opt/homebrew/opt/htslib
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "/opt/homebrew/opt/htslib/include")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "/opt/homebrew/opt/htslib/lib")
endif()
## if not search in  /opt/homebrew/Cellar/htslib/1.20
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "/opt/homebrew/Cellar/htslib")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "/opt/homebrew/Cellar/htslib") 
endif()

if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  message(FATAL_ERROR "Could not find htslib")
endif()


include_directories(${HTSLIB_INCLUDE_DIR} )

target_link_libraries(tfasta  ${HTSLIB_LIBRARY}   logging )


add_executable(${PROJECT_NAME} ${SOURCE_FILE})
target_include_directories(${PROJECT_NAME} PRIVATE "lib")
# target_link_libraries(${PROJECT_NAME} -lz)
# target_link_libraries(${PROJECT_NAME} -lm)
target_link_libraries(${PROJECT_NAME}  ${HTSLIB_LIBRARY}  -lm  -lz logging tfasta)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})


install(TARGETS ${PROJECT_NAME} DESTINATION bin)

install(FILES README.md DESTINATION .)
include(CPack)
include(CPack)
