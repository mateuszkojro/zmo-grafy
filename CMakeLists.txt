cmake_minimum_required(VERSION 3.20)
project(application_name VERSION 0.0.1)

message("")
message("=== CONFIGURING ENVIROMENT ${CMAKE_PROJECT_NAME} - START  ===")
message("")

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_STATIC_LIBS "Build static libraries" OFF)
option(BUILD_TESTS "Build tests" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_DOCS "Build documentation" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

option(WITH_EIGEN "Enable Eigen support" OFF)
option(WITH_ROOT "Enable ROOT support" OFF)
option(WITH_BOOST "Enable BOOST support" OFF)
option(WITH_SCIPLOT "Enable Scilplot support" OFF)
option(BUILD_TYPE "Build type" "Release")


# Operation system detection
message(CHECK_START "Checking operating system")
if(WIN32)
  set(OS Windows)
  add_compile_definitions(MK_OS_WINDOWS)
  message(CHECK_PASS "Windows detected")
elseif(UNIX)
  set(OS Unix)
  add_compile_definitions(MK_OS_UNIX)
  message(CHECK_PASS "Unix detected")
else()
  message(CHECK_FAIL "Unknown operating system")
endif()

# Build type detection
message(CHECK_START "Checking build type")
if(CMAKE_BUILD_TYPE STREQUAL Release OR CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
  set(RELEASE)
  add_compile_definitions(MK_BUILD_RELEASE)
  message(CHECK_PASS "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(RELEASE)
  add_compile_definitions(MK_BUILD_DEBUG)
  message(CHECK_PASS "Debug")
elseif(CMAKE_BUILD_TYPE STREQUAL "")
  set(DEBUG)
  add_compile_definitions(MK_BUILD_DEBUG)
  message(CHECK_PASS "Unknown build type defaulting to debug")
else()
set(DEBUG)
  add_compile_definitions(MK_BUILD_DEBUG)
  message(CHECK_FAIL "Unsuported build type: \"${CMAKE_BUILD_TYPE}\" defaulting to debug")
endif()

# Compiler detection
message(CHECK_START "Checking compiler")
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(COMPILER Clang)
  message(CHECK_PASS "Clang compiler detected")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER GCC)
  message(CHECK_PASS "GNU compiler detected")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(COMPILER Intel)
  message(CHECK_PASS "Intel compiler detected")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(COMPILER MSVC)
  message(CHECK_PASS "MSVC detected")
else()
  message(CHECK_FAIL "Unknown compiler")
endif()

set(CMAKE_CXX_STANDARD 14)

# Check what libraries are available

set(LIBRARIES_TO_LINK)
set(FILES_TO_COPY_TO_BUILD_DIR)

if(${WITH_EIGEN})
message(CHECK_START "Looking for Eigen")
set(EIGEN_DIR ${CMAKE_SOURCE_DIR}/lib/eigen/)
if(EXISTS ${EIGEN_DIR} AND ${WITH_EIGEN})
    list(APPEND CMAKE_MESSAGE_INDENT "    ")
    include_directories(${EIGEN_DIR})
    add_subdirectory(${EIGEN_DIR})
    add_compile_definitions(MK_WITH_EIGEN)
    list(POP_BACK CMAKE_MESSAGE_INDENT)
    message(CHECK_PASS "found")
else()
    message(CHECK_FAIL "not found")
endif()
endif()

if(${WITH_ROOT})
message(CHECK_START "Looking for ROOT")
# Full list of avaliable ROOT libraries 
# https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
find_package(ROOT CONFIG)
if(${ROOT_FOUND} AND ${WITH_ROOT})
    include(${ROOT_USE_FILE})
    include_directories(${ROOT_INCLUDE_DIRS})
    list(APPEND LIBRARIES_TO_LINK ROOT::Core ROOT::Gui)
    add_compile_definitions(MK_WITH_ROOT)
    message(CHECK_PASS "found")
else()
    message(CHECK_FAIL "not found")
endif()
endif()

if(${WITH_EIGEN})
message(CHECK_START "Looking for Sciplot")
set(SCIPLOT_DIR ${CMAKE_SOURCE_DIR}/lib/sciplot/)
if(EXISTS ${SCIPLOT_DIR})
    include_directories(${SCIPLOT_DIR})
    add_compile_definitions(MK_WITH_SCIPLOT)
    message(CHECK_PASS "found")
else()
    message(CHECK_FAIL "not found")
endif()
endif()

if(${WITH_BOOST})
message(CHECK_START "Looking for BOOST")
if(true)
    list(APPEND CMAKE_MESSAGE_INDENT "    ")
    set(BOOST_URL "https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.bz2" CACHE STRING "Boost download URL")

    include(FetchContent)
    FetchContent_Declare(
      Boost
      URL ${BOOST_URL}
    )
    FetchContent_GetProperties(Boost)

    if(NOT Boost_POPULATED)
      message(STATUS "Fetching Boost")
      FetchContent_Populate(Boost)
      message(STATUS "Fetching Boost - done")
      set(BOOST_SOURCE ${boost_SOURCE_DIR})
    else()
      message(STATUS "Boost already fetched")
    endif()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

    set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

    # Define the header-only Boost target
    add_library(Boost::boost INTERFACE IMPORTED GLOBAL)
    target_include_directories(Boost::boost SYSTEM INTERFACE ${BOOST_SOURCE})

    # Disable autolink
    target_compile_definitions(Boost::boost INTERFACE BOOST_ALL_NO_LIB=1)

    target_compile_definitions(Boost::boost INTERFACE
      $<$<CONFIG:Release>:BOOST_DISABLE_ASSERT>
      BOOST_SYSTEM_NO_DEPRECATED
      BOOST_THREAD_VERSION=4
      BOOST_THREAD_USES_CHRONO
      BOOST_THREAD_PROVIDES_EXECUTORS
    )

    list(APPEND LIBRARIES_TO_LINK Boost::boost)
    add_compile_definitions(MK_WITH_BOOST)
    list(POP_BACK CMAKE_MESSAGE_INDENT)
    message(CHECK_PASS "downloaded to: ${BOOST_SOURCE}")
else()
    message(CHECK_FAIL "not found")
endif()
endif()

# List of all source files to compile and link
set(SOURCES src/main.cpp)

# List of folders contatinig header files
set(INCLUDE_DIRS include)

# Include all the folders contatining header files
include_directories(${INCLUDE_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES})


target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${LIBRARIES_TO_LINK})

# foreach(FILE FILES_TO_COPY_TO_BUILD_DIR)
# add_custom_command(
#         TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         FILE
#         $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>
# )
# endforeach()