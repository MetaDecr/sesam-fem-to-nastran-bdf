# fem2bdf -- GeniE/Sesam FEM (T-file) to Nastran BDF converter.
#
# Builds the vendored ExtFEMIO library (DNV GL, github.com/dnvgl/ExtFEMIO)
# as a static library and links the fem2bdf converter against it.
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#
# Tested compilers: MSVC 2017+, GCC 7+, Clang 6+, zig c++ (0.16).

cmake_minimum_required(VERSION 3.12)
project(fem2bdf CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ---------------------------------------------------------------------------
# ExtFEMIO static library (vendored under extfemio/)
# ---------------------------------------------------------------------------
file(GLOB EXT_FEM_SRC CONFIGURE_DEPENDS extfemio/src/fem/*.cpp)
file(GLOB EXT_BDF_SRC CONFIGURE_DEPENDS extfemio/src/bdf/*.cpp)
file(GLOB EXT_ROOT_SRC CONFIGURE_DEPENDS extfemio/src/*.cpp)
file(GLOB EXT_SUP_SRC CONFIGURE_DEPENDS extfemio/src/support/*.cpp)

# The library ships its own command-line tools with a main() each; drop them.
list(REMOVE_ITEM EXT_FEM_SRC
     "${CMAKE_CURRENT_SOURCE_DIR}/extfemio/src/fem/fem2fem.cpp")
list(REMOVE_ITEM EXT_BDF_SRC
     "${CMAKE_CURRENT_SOURCE_DIR}/extfemio/src/bdf/bdf2bdf.cpp")

add_library(extfemio STATIC
    ${EXT_FEM_SRC} ${EXT_BDF_SRC} ${EXT_ROOT_SRC} ${EXT_SUP_SRC})

target_include_directories(extfemio PUBLIC extfemio/include)

# The library carries a std::make_unique back-port guarded by _MY_CPP14_H_;
# with a real C++14 compiler the back-port conflicts, so skip it.
target_compile_definitions(extfemio PUBLIC _MY_CPP14_H_)

# The vendored sources predate modern compilers; silence their warnings.
if(MSVC)
    target_compile_options(extfemio PRIVATE /W0)
else()
    target_compile_options(extfemio PRIVATE -w)
endif()

# ---------------------------------------------------------------------------
# fem2bdf converter
# ---------------------------------------------------------------------------
add_executable(fem2bdf fem2bdf/fem2bdf.cpp)
target_link_libraries(fem2bdf PRIVATE extfemio)

install(TARGETS fem2bdf RUNTIME DESTINATION bin)
