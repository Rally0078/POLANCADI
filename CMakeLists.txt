cmake_minimum_required(VERSION 3.21)

project(polan
LANGUAGES Fortran)
enable_testing()

if(UNIX)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
    add_compile_options(-g -mtune=native -std=legacy)
    add_link_options(-static)
  endif()
endif()
if(WIN32)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
    add_compile_options(-g -mtune=native -std=legacy)
    add_link_options(-static)
  endif()
endif()

add_executable(polan)
add_test(NAME basic COMMAND $<TARGET_FILE:polan> ${CMAKE_CURRENT_SOURCE_DIR}/examples/in.dat)

add_subdirectory(src)

if(UNIX)
  add_test(NAME CompareOut COMMAND bash -c "diff <(tail -n +2 out.dat) <(tail -n +2 ../examples/out.dat)")
  set_tests_properties(CompareOut PROPERTIES FIXTURES_SETUP basic)
endif()
