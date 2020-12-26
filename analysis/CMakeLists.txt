include_directories(${ROOT_INCLUDE_DIRS})

add_executable(Analysis analysis.cpp)

target_link_libraries(Analysis ${ROOT_LIBRARIES})

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Generation.pdf ${CMAKE_CURRENT_BINARY_DIR}/InvMass.pdf
    COMMAND Analysis ${CMAKE_BINARY_DIR}/generation/Histograms.root
    DEPENDS runGeneration Analysis ${CMAKE_BINARY_DIR}/generation/Histograms.root
)

add_custom_target(runAnalysis DEPENDS 
runGeneration ${CMAKE_BINARY_DIR}/generation/Histograms.root 
${CMAKE_CURRENT_BINARY_DIR}/Generation.pdf ${CMAKE_CURRENT_BINARY_DIR}/InvMass.pdf)
