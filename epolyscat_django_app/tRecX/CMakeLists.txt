cmake_minimum_required (VERSION 3.10)

project (tRecX)
set(PROJECT_VERSION "1.0.1")

# Ctype is a minimal build, avoid extra libraries
if(NOT CMAKE_BUILD_TYPE STREQUAL "Ctypes")
    find_package(Boost COMPONENTS system filesystem)    
    if (Boost_FOUND)
        include_directories( ${Boost_INCLUDE_DIR} $ENV{GSL_INCDIR} )
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_BOOST_")
        list(APPEND tRecX_LIBRARIES ${Boost_LIBRARIES} )
        message("\nBoost found -- full functionality")
    else()
        message("\nBoost not found, breaks haCC at present. Specify -DBOOST_ROOT=[path-to-boost] to find")
    endif()
#	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_UI_")
endif()

find_package(MPI)
if (MPI_CXX_FOUND)
    if(_MPI STREQUAL "false")
    message("MPI found on system, but MPI compiliation suppressed by -D_MPI=false")
    else()
    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
    message("MPI found on system, building with ${CMAKE_CXX_COMPILER} , for suppress specify cmake . -D_MPI=false")
    endif()

endif()



set (CMAKE_CXX_STANDARD 11)



set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# set flags for compilation choices
if (CMAKE_BUILD_TYPE STREQUAL "Develop")
    if(NOT ${MPI_CXX_COMPILER} STREQUAL ${CMAKE_CXX_COMPILER})
    message(FATAL_ERROR "not an MPI compiler ${CMAKE_CXX_COMPILER}, cannot build  ${CMAKE_BUILD_TYPE} cmake . -UCMAKE_BUILD_TYPE for choices")
    endif()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg -D_DEVELOP_ -pthread")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    
elseif (CMAKE_BUILD_TYPE STREQUAL "Ctypes")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg -D_DEVELOP_  -pthread -D_NOMPI_")

elseif (CMAKE_BUILD_TYPE STREQUAL "DevelopSimple")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg -D_DEVELOP_  -pthread -D_LOG_ -D_NOMPI_")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")

elseif (CMAKE_BUILD_TYPE STREQUAL "OpenMP")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -pthread -D_NOMPI_")
     
elseif (CMAKE_BUILD_TYPE STREQUAL "Opt0")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -pg -D_DEVELOP_  -pthread ")
     
elseif (CMAKE_BUILD_TYPE STREQUAL "Opt1")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g -pg -D_DEVELOP_  -pthread ")
     
elseif (CMAKE_BUILD_TYPE STREQUAL "Opt2")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -pg -D_DEVELOP_  -pthread ")
     
elseif (CMAKE_BUILD_TYPE STREQUAL "Parallel")
    if(NOT ${MPI_CXX_COMPILER} STREQUAL ${CMAKE_CXX_COMPILER})
    message(FATAL_ERROR "\nNot an MPI compiler ${CMAKE_CXX_COMPILER}, cannot build  ${CMAKE_BUILD_TYPE}\ncmake . -UCMAKE_BUILD_TYPE for choices")
    endif()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -pthread ")
     
elseif (CMAKE_BUILD_TYPE STREQUAL "Scalar")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_NOMPI_  -pthread ")
     
elseif(NOT CMAKE_BUILD_TYPE)
    if(NOT ${MPI_CXX_COMPILER} STREQUAL ${CMAKE_CXX_COMPILER})
    if(MPI_CXX_FOUND)     
    message("MPI found, but not an MPI compiler: ${CMAKE_CXX_COMPILER}")
    message("Build sequential, for Parallel manually define -DCMAKE_CXX_COMPILER=[mpi-compile-script] -DCMAKE_BUILD_TYPE=Parallel")
    endif()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_NOMPI_  -pthread ")

    endif()

else()
     message("\nUnknown CMAKE_BUILD_TYPE: \"" ${CMAKE_BUILD_TYPE} "\"")
     message("Set by -DCMAKE_BUILD_TYPE=[value], allowed values: Parallel|Scalar|Develop|DevelopSimple|OptN,N=0,1,2")
     message("Reset default by -UCMAKE_BUILD_TYPE\n"  )
     message(FATAL_ERROR "Build type unknown")
endif()

# compiler-specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_EXE_LINKER_FLAGS "-mkl")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -diag-disable 858 ")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()

include("${CMAKE_SOURCE_DIR}/LIBRARIES/CMake.include")

if(NOT CMAKE_BUILD_TYPE)
	message("\nDefault build, other options: -DCMAKE_BUILD_TYPE=Parallel|Scalar|Develop|DevelopSimple, return to default by -UCMAKE_BUILD_TYPE\n")
else()
       message("\nBuild \"" ${CMAKE_BUILD_TYPE} "\" with C++ flags " ${CMAKE_CXX_FLAGS} 
       		    ", return to default by -UCMAKE_BUILD_TYPE\n")
endif()

include_directories("${CMAKE_SOURCE_DIR}/LIBRARIES/Eigen")
include_directories("${CMAKE_SOURCE_DIR}/LIBRARIES/alglib")
include_directories("${CMAKE_SOURCE_DIR}/TOOLS")
#include_directories("${CMAKE_SOURCE_DIR}/TOOLS/zslapack")
include_directories("${CMAKE_SOURCE_DIR}/SRCS/core")
include_directories("${CMAKE_SOURCE_DIR}/SRCS/user")


add_subdirectory(LIBRARIES)
add_subdirectory(SRCS/core)
add_subdirectory(SRCS/user)
add_subdirectory(TOOLS)
#add_subdirectory(TOOLS/zslapack)

message("            Compiler: " ${CMAKE_CXX_COMPILER} ", change by -DCMAKE_CXX_COMPILER=[path_to_complier]\n")                     
message(STATUS "Linking to libraries ${tRecX_LIBRARIES}")

add_executable(tRecX SRCS/main/main_trecx.cpp SRCS/main/run_trecx.cpp SRCS/main/tRecX.cpp)
include_directories(${MKL_INCLUDE_DIRS})
target_link_libraries(tRecX LINK_PUBLIC  tRecX_core tRecX_user tRecX_tools ${tRecX_LIBRARIES} ${MKL_LIBRARIES})

add_executable(Harmonics SRCS/main/main_harmonics.cpp SRCS/main/tRecX.cpp)
target_link_libraries(Harmonics LINK_PUBLIC tRecX_core tRecX_tools tRecX_user ${tRecX_LIBRARIES})


if(CMAKE_CXX_FLAGS MATCHES "-O0")
add_custom_command(TARGET tRecX POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red "--- build w/o optimization - tRecX will be slow ---"
                   COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --green "consider, e.g. GNU version 7, compiles with -O3")
endif()

