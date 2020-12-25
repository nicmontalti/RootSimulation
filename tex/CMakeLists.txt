include(UseLATEX.cmake)

ADD_LATEX_DOCUMENT(
    Relazione.tex 
    TARGET_NAME Relazione
    DEPENDS Analysis Generation runAnalysis ${CMAKE_BINARY_DIR}/analysis/Generation.pdf ${CMAKE_BINARY_DIR}/analysis/InvMass.pdf)
