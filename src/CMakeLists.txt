add_library(
    matrix
    matrix_mkl.cpp
)

set(MKL_INCLUDE_DIRS $ENV{MKLROOT}/include)
set(MKL_LIBRARY_DIR $ENV{MKLROOT}/lib/intel64)

target_include_directories(
    matrix PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MKL_INCLUDE_DIRS}
)

target_link_libraries(
    matrix PUBLIC
    ${MKL_LIBRARY_DIR}/libmkl_rt.so
)



add_executable(
    opt 
    main.cpp
    vec3f64.cpp
    hessian.cpp
    geomopt.cpp
    cell.cpp
    rfo.cpp
    bfgs.cpp
    gdiis.cpp
    gediis.cpp
    line_search.cpp
)

target_link_libraries(
    opt PUBLIC
    matrix
)