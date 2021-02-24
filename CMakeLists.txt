cmake_minimum_required(VERSION 2.8)
project(HeatMethod_igl)

#Add the project cmake scripts to the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#Download external dependencies NOTE: If the user has one of these libs it shouldn't be downloaded again.
include(${CMAKE_MODULE_PATH}/DownloadProject.cmake)
if (CMAKE_VERSION VERSION_LESS 3.2)
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
else()
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
endif()

##vcglib devel branch
download_project(PROJ		    vcglib_devel
                 GIT_REPOSITORY     https://github.com/cnr-isti-vclab/vcglib.git
                 GIT_TAG	    devel
                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
)

#libigl
download_project(PROJ               LIBIGL
                 GIT_REPOSITORY     https://github.com/libigl/libigl.git
                 GIT_TAG	    master
                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
)
set(LIBIGL_INCLUDE_DIR ${LIBIGL_SOURCE_DIR}/include)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
find_package(LIBIGL REQUIRED QUIET)

#Polyscope
#download_project(PROJ               POLYSCOPE
#                 GIT_REPOSITORY     https://github.com/nmwsharp/polyscope.git
#                 GIT_TAG	    master
#                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
#)
#add_subdirectory(${POLYSCOPE_SOURCE_DIR})



##Eigen 3 NOTE: Eigen is required on the system the code is ran
find_package(Eigen3 3.3 REQUIRED)

#Add the project sources
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} )

set(LIBIGL_INCLUDE_DIR "/home/iason/Coding/Libraries/libigl/include")

target_include_directories(${PROJECT_NAME}
        PRIVATE ${vcglib_devel_SOURCE_DIR}
        PRIVATE ${LIBIGL_SOURCE_DIR}
#        PRIVATE ${POLYSCOPE_SOURCE_DIR}
        )

#Use C++17
if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++20")
else(MSVC)
        target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
endif(MSVC)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui
)


