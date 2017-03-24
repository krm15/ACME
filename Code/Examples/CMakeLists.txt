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

ADD_EXECUTABLE(cellPreprocess CellPreprocess.cxx)
TARGET_LINK_LIBRARIES(cellPreprocess ${ITK_LIBRARIES})

ADD_EXECUTABLE(intensityThreshold IntensityThreshold.cxx )
TARGET_LINK_LIBRARIES(intensityThreshold ${ITK_LIBRARIES})

ADD_EXECUTABLE(resample Resample.cxx)
TARGET_LINK_LIBRARIES(resample ${ITK_LIBRARIES})

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

ADD_EXECUTABLE( membraneSegmentation MembraneSegmentationImageFilter.cxx)
TARGET_LINK_LIBRARIES( membraneSegmentation ${ITK_LIBRARIES})

ADD_EXECUTABLE( membraneSegmentationWithMarkers MembraneSegmentationWithMarkersImageFilter.cxx)
TARGET_LINK_LIBRARIES( membraneSegmentationWithMarkers ${ITK_LIBRARIES})

ADD_EXECUTABLE(cvtGenerator CVTGenerator.cxx)
TARGET_LINK_LIBRARIES(cvtGenerator ${ITK_LIBRARIES} Random )

# For comparing the three segmentation versions
ADD_EXECUTABLE(membraneSegmentationEvaluation MembraneSegmentationEvaluation.cxx)
TARGET_LINK_LIBRARIES(membraneSegmentationEvaluation ${ITK_LIBRARIES} Random)

# For computing the tables in the paper
ADD_EXECUTABLE(labelMapOverlapMeasures LabelMapOverlapMeasures.cxx)
TARGET_LINK_LIBRARIES(labelMapOverlapMeasures ${ITK_LIBRARIES} )

ADD_EXECUTABLE( MorphologicalErosionOnLabelImageFilter MorphologicalErosionOnLabelImageFilter.cxx )
TARGET_LINK_LIBRARIES( MorphologicalErosionOnLabelImageFilter ${ITK_LIBRARIES} )

ADD_EXECUTABLE(cellSegmentationStatistics CellSegmentationStatistics.cxx )
TARGET_LINK_LIBRARIES(cellSegmentationStatistics ${ITK_LIBRARIES})

ADD_EXECUTABLE( DistanceFromMask DistanceFromMask.cxx )
TARGET_LINK_LIBRARIES( DistanceFromMask ${ITK_LIBRARIES} )

ADD_EXECUTABLE(sizeThresh SizeThreshold.cxx )
TARGET_LINK_LIBRARIES(sizeThresh ${ITK_LIBRARIES} )

