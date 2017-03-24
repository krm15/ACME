ADD_EXECUTABLE( DistanceFromMask DistanceFromMask.cxx )
TARGET_LINK_LIBRARIES( DistanceFromMask ${ITK_LIBRARIES} )

ADD_EXECUTABLE(sizeThresh SizeThreshold.cxx )
TARGET_LINK_LIBRARIES(sizeThresh ${ITK_LIBRARIES} )

ADD_EXECUTABLE(intensityThreshold IntensityThreshold.cxx )
TARGET_LINK_LIBRARIES(intensityThreshold ${ITK_LIBRARIES})

ADD_EXECUTABLE(resample Resample.cxx)
TARGET_LINK_LIBRARIES(resample ${ITK_LIBRARIES})

ADD_EXECUTABLE( MorphologicalErosionOnLabelImageFilter MorphologicalErosionOnLabelImageFilter.cxx )
TARGET_LINK_LIBRARIES( MorphologicalErosionOnLabelImageFilter ${ITK_LIBRARIES} )

ADD_EXECUTABLE( MorphologicalDilationOnLabelImageFilter MorphologicalDilationOnLabelImageFilter.cxx )
TARGET_LINK_LIBRARIES( MorphologicalDilationOnLabelImageFilter ${ITK_LIBRARIES} )
