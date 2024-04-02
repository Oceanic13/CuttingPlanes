cmake_minimum_required(VERSION 3.5)

project(CuttingPlanes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/libs")
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

include_directories("libs/")
file(GLOB_RECURSE LIBS_CPP "${LIBS_DIR}/*.cpp")
file(GLOB_RECURSE LIBS_H "${LIBS_DIR}/*.h")

add_executable(${PROJECT_NAME}
    main.cpp

    Utils.h
    SoPlexSolver.h
    CuttingPlanes.h
    AlgLibSimplexSolver.h
    ToblexSolver.h
    MixedIntegerLinearProgram.h

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

target_link_libraries (${PROJECT_NAME} PUBLIC
    Boost::boost
    libsoplex
    nlohmann_json::nlohmann_json
    Eigen3::Eigen
)

###########
# testing #
###########
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()
enable_testing()
include(GoogleTest)
find_package(Boost REQUIRED COMPONENTS system filesystem)
add_executable(${PROJECT_NAME}-test
    tests.cpp
)
target_link_libraries(${PROJECT_NAME}-test
    gtest gtest_main
    Boost::boost
    libsoplex
    nlohmann_json::nlohmann_json
    Eigen3::Eigen
)
gtest_add_tests(TARGET ${PROJECT_NAME}-test
    EXTRA_ARGS "--gtest_color=yes"
)
set_target_properties(${PROJECT_NAME}-test PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
)
