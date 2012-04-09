PROJECT( MembraneReconstruction )

# CMake >= 2.6 MUST be installed
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

IF( COMMAND CMAKE_POLICY )
  CMAKE_POLICY( VERSION 2.4 )
  CMAKE_POLICY( SET CMP0005 OLD )
  CMAKE_POLICY( SET CMP0003 NEW )
ENDIF(COMMAND CMAKE_POLICY)

# On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security
# warnings. Copied from ITK CMakeLists.
IF( WIN32 )
  IF( NOT BORLAND )
    IF( NOT CYGWIN )
      IF( NOT MINGW )
          ADD_DEFINITIONS(
            -D_CRT_FAR_MAPPINGS_NO_DEPRECATE
            -D_CRT_IS_WCTYPE_NO_DEPRECATE
            -D_CRT_MANAGED_FP_NO_DEPRECATE
            -D_CRT_NONSTDC_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE_GLOBALS
            -D_CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
            -D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
            -D_CRT_VCCLRIT_NO_DEPRECATE
            -D_SCL_SECURE_NO_DEPRECATE
            )
      ENDIF( NOT MINGW )
    ENDIF( NOT CYGWIN )
  ENDIF( NOT BORLAND )
ENDIF( WIN32 )


#---------------------------------------------------------------------------
# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  IF( ITK_MINOR_VERSION GREATER 14 )
    ADD_DEFINITIONS( -DITKLEVELSETCODE )
  ELSE( ITK_MINOR_VERSION GREATER 14 )
    REMOVE_DEFINITIONS( -DITKLEVELSETCODE )
  ENDIF( ITK_MINOR_VERSION GREATER 14 )
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


#---------------------------------------------------------------------------
# Define

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY
 ${MembraneReconstruction_BINARY_DIR}/bin
 CACHE PATH "Single output directory for building all executables."
)

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
 ${MembraneReconstruction_BINARY_DIR}/lib
 CACHE PATH "Single output directory for building all libraries."
)

ADD_LIBRARY( Random ${MembraneReconstruction_SOURCE_DIR}/ConfocalImageGenerator/itkThreadSafeMersenneTwisterRandomVariateGenerator.cxx )

SET( DIRECTORIES
  CellPreprocess
  PlanarityFilter
  TensorVoting
  WatershedImageFilter
  Examples
)

SUBDIRS( ${DIRECTORIES} )


