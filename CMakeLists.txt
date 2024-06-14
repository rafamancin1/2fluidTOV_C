cmake_minimum_required(VERSION 3.22)
project(2FLUIDTOV_C)


add_executable(TWO_FLUID_TOV src/main.cpp)

#if(NOT TARGET spdlog)
    # Stand-alone build
#find_package(spdlog REQUIRED)
#endif()
#target_link_libraries(TWO_FLUID_TOV PRIVATE spdlog::spdlog)
#add_library( MODULE)