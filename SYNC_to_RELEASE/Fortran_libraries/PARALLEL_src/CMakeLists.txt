# Get sources from each subdirectory
FILE(GLOB COMM_SRC_FILES COMM_src/*.[fF]90)
FILE(GLOB SPH_SHELL_SRC_FILES SPH_SHELL_src/*.[fF]90)

SET(CALYPSO_PARALLEL_SRC_FILES ${COMM_SRC_FILES} ${SPH_SHELL_SRC_FILES})

