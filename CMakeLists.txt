cmake_minimum_required(VERSION 2.8)

project(RKDG-1D)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SRC_LIST src/RKDG-1D.cpp)

FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h
  )

list(REMOVE_ITEM ALLFILES "${CMAKE_SOURCE_DIR}/src/Flux/FluxHLLE.cpp")
list(REMOVE_ITEM ALLFILES "${CMAKE_SOURCE_DIR}/src/Flux/FluxKIR.cpp")
list(REMOVE_ITEM ALLFILES "${CMAKE_SOURCE_DIR}/src/Limiter/LimiterWENOH.cpp")  
  
add_executable(${PROJECT_NAME} ${ALLFILES})

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)

#MESSAGE( STATUS "variable is: " ${CMAKE_CONFIGURATION_TYPE} )

include_directories(${dirlist})

#file(GLOB_RECURSE REMOVE_CMAKE "build/*")
#list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/build/CMakeFiles/3.9.1/CompilerIdCXX/CMakeCXXCompilerId.cpp")




