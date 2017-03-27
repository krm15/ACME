INCLUDE_DIRECTORIES( BEFORE
  ${ACME_SOURCE_DIR}/CellPreprocess
  ${ACME_BINARY_DIR}/CellPreprocess
  ${ACME_SOURCE_DIR}/PlanarityFilter
  ${ACME_BINARY_DIR}/PlanarityFilter
  ${ACME_SOURCE_DIR}/TensorVoting
  ${ACME_BINARY_DIR}/TensorVoting
  ${ACME_SOURCE_DIR}/ConfocalImageGenerator
  ${ACME_BINARY_DIR}/ConfocalImageGenerator
  ${ACME_SOURCE_DIR}/ImageProcessingFilters
  ${ACME_BINARY_DIR}/ImageProcessingFilters
)

SUBDIRS(Tools)
SUBDIRS(Segmentation)

ADD_EXECUTABLE(cellPreprocess CellPreprocess.cxx)
TARGET_LINK_LIBRARIES(cellPreprocess ${ITK_LIBRARIES})

ADD_EXECUTABLE(membraneOrientationBasedHessianCalculator MembraneOrientationBasedHessianCalculator.cxx)
TARGET_LINK_LIBRARIES(membraneOrientationBasedHessianCalculator ${ITK_LIBRARIES})

ADD_EXECUTABLE(multiscalePlateMeasureImageFilter MultiscalePlateMeasureImageFilter.cxx)
TARGET_LINK_LIBRARIES(multiscalePlateMeasureImageFilter ${ITK_LIBRARIES})

ADD_EXECUTABLE(ballVotingIntersectingCirclesExample2D BallVotingIntersectingCirclesExample2D.cxx)
TARGET_LINK_LIBRARIES(ballVotingIntersectingCirclesExample2D ${ITK_LIBRARIES} Random)

ADD_EXECUTABLE(membraneVotingField2D MembraneVotingField2D.cxx)
TARGET_LINK_LIBRARIES(membraneVotingField2D ${ITK_LIBRARIES} Random)

ADD_EXECUTABLE(membraneVotingField3D MembraneVotingField3D.cxx)
TARGET_LINK_LIBRARIES(membraneVotingField3D ${ITK_LIBRARIES} Random)

ADD_EXECUTABLE(multiscalePlanarityAndVoting3D MultiscalePlanarityAndVoting3D.cxx)
TARGET_LINK_LIBRARIES(multiscalePlanarityAndVoting3D ${ITK_LIBRARIES} Random)

# For computing the tables in the paper
ADD_EXECUTABLE(labelMapOverlapMeasures LabelMapOverlapMeasures.cxx)
TARGET_LINK_LIBRARIES(labelMapOverlapMeasures ${ITK_LIBRARIES} )

ADD_EXECUTABLE(cellSegmentationStatistics CellSegmentationStatistics.cxx )
TARGET_LINK_LIBRARIES(cellSegmentationStatistics ${ITK_LIBRARIES})



