PROJECT( ACME )

# CMake >= 2.6 MUST be installed
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

IF( COMMAND CMAKE_POLICY )
  CMAKE_POLICY( VERSION 2.4 )
  CMAKE_POLICY( SET CMP0005 OLD )
  CMAKE_POLICY( SET CMP0003 NEW )
ENDIF(COMMAND CMAKE_POLICY)

#---------------------------------------------------------------------------
# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


#---------------------------------------------------------------------------
# Define

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY
 ${ACME_BINARY_DIR}/bin
 CACHE PATH "Single output directory for building all executables."
)

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
 ${ACME_BINARY_DIR}/lib
 CACHE PATH "Single output directory for building all libraries."
)

ADD_LIBRARY( Random ${ACME_SOURCE_DIR}/ConfocalImageGenerator/itkThreadSafeMersenneTwisterRandomVariateGenerator.cxx )

SET( DIRECTORIES
  CellPreprocess
  PlanarityFilter
  TensorVoting
  ImageProcessingFilters
  Examples
)

SUBDIRS( ${DIRECTORIES} )


