include_directories(${ROOT_INCLUDE_DIRS})
include_directories(../include)

add_library(Particle ParticleType.cpp ResonanceType.cpp Particle.cpp)

add_executable(Generation main.cpp)
target_link_libraries(Generation Particle)
target_link_libraries(Generation ${ROOT_LIBRARIES})

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Histograms.root
    COMMAND Generation
    DEPENDS Generation
)

add_custom_target(runGeneration DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Histograms.root)