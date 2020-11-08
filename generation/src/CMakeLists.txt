include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(Particle ParticleType.cpp ResonanceType.cpp Particle.cpp)

add_executable(Generation main.cpp)
set_target_properties(Generation PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_NAME Generation.out)
target_link_libraries(Generation Particle)
target_link_libraries(Generation ${ROOT_LIBRARIES})

