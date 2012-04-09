INCLUDE_DIRECTORIES( BEFORE
  ${MembraneReconstruction_SOURCE_DIR}/CellPreprocess
  ${MembraneReconstruction_BINARY_DIR}/CellPreprocess
  ${MembraneReconstruction_SOURCE_DIR}/PlanarityFilter
  ${MembraneReconstruction_BINARY_DIR}/PlanarityFilter
  ${MembraneReconstruction_SOURCE_DIR}/TensorVoting
  ${MembraneReconstruction_BINARY_DIR}/TensorVoting
  ${MembraneReconstruction_SOURCE_DIR}/ConfocalImageGenerator
  ${MembraneReconstruction_BINARY_DIR}/ConfocalImageGenerator
  ${MembraneReconstruction_SOURCE_DIR}/WatershedImageFilter
  ${MembraneReconstruction_BINARY_DIR}/WatershedImageFilter
)

ADD_EXECUTABLE(cellPreprocess CellPreprocess.cxx)
TARGET_LINK_LIBRARIES(cellPreprocess ${ITK_LIBRARIES})

ADD_EXECUTABLE(resample Resample.cxx)
TARGET_LINK_LIBRARIES(resample ${ITK_LIBRARIES})

ADD_EXECUTABLE(multiscalePlateMeasureImageFilter MultiscalePlateMeasureImageFilter.cxx)
TARGET_LINK_LIBRARIES(multiscalePlateMeasureImageFilter ${ITK_LIBRARIES})

ADD_EXECUTABLE(ballVotingIntersectingCirclesExample2D BallVotingIntersectingCirclesExample2D.cxx)
TARGET_LINK_LIBRARIES(ballVotingIntersectingCirclesExample2D ${ITK_LIBRARIES})

ADD_EXECUTABLE(membraneVotingField3D MembraneVotingField3D.cxx)
TARGET_LINK_LIBRARIES(membraneVotingField3D ${ITK_LIBRARIES})

ADD_EXECUTABLE( membraneSegmentation MembraneSegmentationImageFilter.cxx)
TARGET_LINK_LIBRARIES( membraneSegmentation ${ITK_LIBRARIES})

ADD_EXECUTABLE(cvtGenerator CVTGenerator.cxx)
TARGET_LINK_LIBRARIES(cvtGenerator ${ITK_LIBRARIES} Random )

