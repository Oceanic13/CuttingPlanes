cmake_minimum_required(VERSION 3.5)

project(CuttingPlanes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/libs")
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

include_directories("libs/")
file(GLOB_RECURSE LIBS_CPP "${LIBS_DIR}/*.cpp")
file(GLOB_RECURSE LIBS_H "${LIBS_DIR}/*.h")

add_executable(CuttingPlanes
    main.cpp

    Utils.h
    SoPlexSolver.h
    CuttingPlanes.h
    AlgLibSimplexSolver.h
    ToblexSolver.h
    LinearProgram.h




    #${LIBS_CPP}
    #${LIBS_H}
)

include(GNUInstallDirs)
install(TARGETS CuttingPlanes
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)

if (NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 3.4
        SOURCE_DIR "${EXTERNAL_DIR}/Eigen"
        )
    FetchContent_MakeAvailable(eigen)
endif()

if (NOT TARGET libsoplex)
    FetchContent_Declare(soplex
        GIT_REPOSITORY https://github.com/scipopt/soplex
        SOURCE_DIR "${EXTERNAL_DIR}/SoPlex"
    )
    FetchContent_MakeAvailable(soplex)
endif()


if (NOT TARGET nlohmann::json)
    FetchContent_Declare(json
        URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
        #SOURCE_DIR "${EXTERNAL_DIR}/Json"
    )
    FetchContent_MakeAvailable(json)
endif()

find_package(Boost REQUIRED COMPONENTS system filesystem)
#include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries (CuttingPlanes PUBLIC Boost::boost libsoplex nlohmann_json::nlohmann_json Eigen3::Eigen)


