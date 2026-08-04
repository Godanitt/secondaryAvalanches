cmake_minimum_required(VERSION 3.16)

if(APPLE)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "" FORCE)
endif()

project(secondaryAvalanches LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist MathCore)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

add_executable(secondaryAvalanches
  secondaryAvalanches.cxx
  src/aux/PhotonModel.cxx
  src/aux/QuantumEfficiency.cxx
)

target_include_directories(secondaryAvalanches PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/aux
)

target_link_libraries(secondaryAvalanches PRIVATE
  Garfield::Garfield
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::MathCore
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(secondaryAvalanches PRIVATE
    -O3 -Wall -Wextra -Wpedantic
    # Garfield++ implements photon transport in AvalancheMicroscopic but keeps
    # TransportPhoton as an internal method. We intentionally use that method
    # in this target without modifying the installed Garfield++ headers.
  )
else()
  message(FATAL_ERROR
    "secondaryAvalanches requires GCC or Clang because Garfield++ keeps "
    "AvalancheMicroscopic::TransportPhoton private.")
endif()
