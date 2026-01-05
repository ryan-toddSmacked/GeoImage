cmake_minimum_required(VERSION 3.14)

project(GeoImage
    VERSION 0.1.0
    DESCRIPTION "A C++ library for working with GeoTIFF images"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include CPM for package management
include(cmake/CPM.cmake)

# Fetch libtiff from GitLab
CPMAddPackage(
    NAME libtiff
    GIT_REPOSITORY https://gitlab.com/libtiff/libtiff.git
    GIT_TAG v4.6.0
    OPTIONS
        "tiff-tools OFF"
        "tiff-tests OFF"
        "tiff-contrib OFF"
        "tiff-docs OFF"
)

# GeoImage library
add_library(GeoImage)

target_sources(GeoImage
    PRIVATE
        src/GeoImage.cpp
)

target_include_directories(GeoImage
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(GeoImage
    PUBLIC
        tiff
)

# Optional: Build example/test executable
option(GEOIMAGE_BUILD_EXAMPLES "Build example executables" ON)

if(GEOIMAGE_BUILD_EXAMPLES)
    add_executable(geoimage_example examples/main.cpp)
    target_link_libraries(geoimage_example PRIVATE GeoImage)
endif()
