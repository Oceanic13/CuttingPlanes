cmake_minimum_required(VERSION 3.5)

project(CuttingPlanes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(CuttingPlanes main.cpp
    SimplexWrapper.h
    Typedefs.h
    CuttingPlanes.h
    CuttingPlanes.cpp)

include(GNUInstallDirs)
install(TARGETS CuttingPlanes
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

if (NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 3.4
        SOURCE_DIR "${EXTERNAL_DIR}/Eigen"
        )
    FetchContent_MakeAvailable(eigen)
endif()

if (NOT TARGET SoPlex::SoPlex)
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


set(SOURCE_FILES
    main.cpp)

target_link_libraries (CuttingPlanes PUBLIC nlohmann_json::nlohmann_json Eigen3::Eigen)


