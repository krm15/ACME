/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 870 $  // Revision of last commit
  Date: $Date: 2009-11-23 13:44:26 -0500 (Mon, 23 Nov 2009) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#ifndef __itkCellSegmentationStatistics_txx
#define __itkCellSegmentationStatistics_txx

#include "itkCellSegmentationStatistics.h"

namespace itk
{
//  Software Guide : BeginCodeSnippet
template < class TFeatureImage, class TInputImage, class TSegmentImage >
CellSegmentationStatistics< TFeatureImage, TInputImage, TSegmentImage >
::CellSegmentationStatistics()
{
  m_ComputeSpatialMap = false;
  m_Shape = ShapeFilterType::New();
  m_RawImg = 0;

  m_LargestCellRadius = 5.0;
}


template < class TFeatureImage,class TInputImage, class TSegmentImage >
void
CellSegmentationStatistics< TFeatureImage, TInputImage, TSegmentImage >::
GenerateData()
{
  ImageSizeType inputSize = m_Segmentation->GetLargestPossibleRegion().GetSize();
  ImageSpacingType spacing = m_Segmentation->GetSpacing();
  ImagePointType origin = m_Segmentation->GetOrigin();

  // convert the image into a collection of objects
  typename ShapeConverterType::Pointer shapeConverter =
    ShapeConverterType::New();
  shapeConverter->SetInput ( m_Segmentation );
  shapeConverter->SetBackgroundValue ( 0 );
  shapeConverter->Update();

  typename ShapeLabelMapType::Pointer shapeLabelMap =
    shapeConverter->GetOutput();

  m_NumOfLabels = shapeLabelMap->GetNumberOfLabelObjects();
  //std::cout << m_NumOfLabels << std::endl;

  LabelGeometryFilterPointer geomFilter = LabelGeometryFilterType::New();
  geomFilter->SetInput( m_Segmentation );
  geomFilter->CalculateOrientedBoundingBoxOn();
  geomFilter->Update();

  m_LabelLookup.resize( m_NumOfLabels );
  m_Centroid.resize ( m_NumOfLabels );
  m_Size.resize ( m_NumOfLabels );
  m_PhysicalSize.resize ( m_NumOfLabels );
  m_BinaryMoments.resize ( m_NumOfLabels );
  m_EllipsoidDiameter.resize ( m_NumOfLabels );
  m_BinaryPrincipalAxis.resize ( m_NumOfLabels );
  m_Region.resize ( m_NumOfLabels );
  m_EquivalentSphericalRadius.resize ( m_NumOfLabels );
  m_Elongation.resize ( m_NumOfLabels );

  m_CenterOfGravity.resize ( m_NumOfLabels );
  m_Mean.resize ( m_NumOfLabels );
  m_Median.resize ( m_NumOfLabels );
  m_Moments.resize ( m_NumOfLabels );
  m_SpatialMap.resize( m_NumOfLabels );

  m_OrientedBoundingBoxSize.resize( m_NumOfLabels );

  typename StatConverterType::Pointer statConverter;
  typename StatLabelMapType::Pointer statLabelMap ;
  if ( m_RawImg )
    {
    statConverter = StatConverterType::New();
    statConverter->SetInput ( m_Segmentation );
    statConverter->SetFeatureImage ( m_RawImg );
    statConverter->SetBackgroundValue ( 0 );
    try
      {
      statConverter->Update();
      }
    catch ( itk::ExceptionObject e )
      {
      std::cerr << "Error: " << e << std::endl;
      }
    statLabelMap = statConverter->GetOutput();
    }

  unsigned int label = 0;
  for( unsigned int i = 0; i < m_NumOfLabels; i++ )
    {
    m_LabelLookup[label] = shapeLabelMap->GetLabels()[i];

    vnlVectorType zero ( ImageDimension, 0 );
    const ShapeLabelObjectType *shapeObject
      = shapeLabelMap->GetLabelObject ( m_LabelLookup[label] );

    m_Size[label] = shapeObject->Size();
    m_PhysicalSize[label] = shapeObject->GetPhysicalSize();
    m_Centroid[label] = shapeObject->GetCentroid().GetVnlVector();
    m_BinaryMoments[label] = shapeObject->GetPrincipalMoments().GetVnlVector();
    m_EllipsoidDiameter[label] = shapeObject->GetEquivalentEllipsoidDiameter().GetVnlVector();
    m_BinaryPrincipalAxis[label] = shapeObject->GetPrincipalAxes().GetVnlMatrix();
    m_Region[label] = shapeObject->GetBoundingBox();
    m_EquivalentSphericalRadius[label] = shapeObject->GetEquivalentSphericalRadius();
    m_Elongation[label] = shapeObject->GetElongation();

    m_OrientedBoundingBoxSize[label] = geomFilter->GetOrientedBoundingBoxSize(label);

    if ( m_RawImg )
      {
      const StatLabelObjectType *statObject
        = statLabelMap->GetLabelObject ( m_LabelLookup[label] );
      m_CenterOfGravity[label] = statObject->GetCenterOfGravity().GetVnlVector();
      m_Mean[label] = statObject->GetMean();
      m_Median[label] = statObject->GetMedian();
      m_Moments[label] = statObject->GetPrincipalMoments().GetVnlVector();
      }

    // Compute major and minor axis here

    label++;
    }

  if ( m_ComputeSpatialMap )
  {
  vnlVectorType d;

    for( unsigned int i = 0; i < m_NumOfLabels; i++ )
    {
      for( unsigned int j = i; j < m_NumOfLabels; j++ )
      {
        d = ( m_Centroid[i] - m_Centroid[j] );
        if ( d.magnitude() < 2*m_LargestCellRadius )
        {
          m_SpatialMap[i].push_back( j );
          if ( i != j)
          {
            m_SpatialMap[j].push_back( i );
          }
        }
      }
    }
  }

}

template < class TFeatureImage,class TInputImage, class TSegmentImage >
void
CellSegmentationStatistics< TFeatureImage, TInputImage, TSegmentImage >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name: " << GetNameOfClass( ) << std::endl;
}

} /* end namespace itk */

#endif
