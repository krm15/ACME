ADD_EXECUTABLE(cvtGenerator CVTGenerator.cxx)
TARGET_LINK_LIBRARIES(cvtGenerator ${ITK_LIBRARIES} Random )

# For comparing the three segmentation versions
ADD_EXECUTABLE(membraneSegmentationEvaluation MembraneSegmentationEvaluation.cxx)
TARGET_LINK_LIBRARIES(membraneSegmentationEvaluation ${ITK_LIBRARIES} Random)
