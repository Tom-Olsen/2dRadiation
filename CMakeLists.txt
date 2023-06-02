project(Radiation2d)
cmake_minimum_required(VERSION 3.16)

# Set output path:
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Enable Debug Symbols:
set(CMAKE_BUILD_TYPE Debug) # has had no impact on performance yet

# Set cpp version:
set(CMAKE_CXX_STANDARD 20)

# Find packages:
find_package(OpenMP REQUIRED)

# Find all src files in src folder:
file(GLOB srcs src/*.cpp)

# Executables to generate and their dependencies:
add_executable(main.out exe/main.cpp ${srcs})
add_executable(test.out exe/test.cpp ${srcs})

# c++ flags:
# target_compile_options(main PUBLIC)
target_compile_options(main.out PUBLIC -O3 -ffast-math)
target_compile_options(test.out PUBLIC -O3 -ffast-math)

# libraries to link:
target_link_libraries(main.out OpenMP::OpenMP_CXX)
target_link_libraries(test.out OpenMP::OpenMP_CXX)