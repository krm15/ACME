/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 597 $  // Revision of last commit
  Date: $Date: 2009-08-25 16:43:34 -0400 (Tue, 25 Aug 2009) $  // Date of last commit
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
#ifndef __itkCellSegmentationStatistics_h
#define __itkCellSegmentationStatistics_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkCastImageFilter.h"
#include "itkImageRegion.h"
#include "itkRegion.h"
#include "itkIndex.h"
#include "itkSize.h"
#include <map>

#include "itkLabelObject.h"
#include "itkShapeLabelObject.h"
#include "itkStatisticsLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkShapeRelabelImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelGeometryImageFilter.h"

#include "itkNumericTraits.h"
#include "itkVector.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <list>

namespace itk
{
template < class TFeatureImage, class TInputImage, class TSegmentImage >
class ITK_EXPORT CellSegmentationStatistics : public Object
{
  public:
    typedef CellSegmentationStatistics           Self;
    typedef Object Superclass;
    typedef SmartPointer<Self>                               Pointer;
    typedef SmartPointer<const Self>                         ConstPointer;

    itkStaticConstMacro ( ImageDimension, unsigned int,
                          TSegmentImage::ImageDimension );

    /** Method for creation through object factory */
    itkNewMacro ( Self );

    /** Run-time type information */
    itkTypeMacro ( CellSegmentationStatistics, Object );

    /** Display */
    void PrintSelf ( std::ostream& os, Indent indent ) const;

    typedef TFeatureImage                           FeatureImageType;
    typedef typename FeatureImageType::Pointer      FeatureImagePointer;
    typedef typename FeatureImageType::ConstPointer FeatureImageConstPointer;

    typedef TInputImage                             ImageType;
    typedef typename ImageType::Pointer             ImagePointer;
    typedef typename ImageType::ConstPointer        ImageConstPointer;
    typedef typename ImageType::PixelType           ImagePixelType;
    typedef typename ImageType::RegionType          ImageRegionType;
    typedef typename ImageType::SizeType            ImageSizeType;
    typedef typename ImageSizeType::SizeValueType   ImageSizeValueType;
    typedef typename ImageType::SpacingType         ImageSpacingType;
    typedef typename ImageType::IndexType           ImageIndexType;
    typedef typename ImageType::PointType           ImagePointType;

    typedef TSegmentImage                           SegmentImageType;
    typedef typename SegmentImageType::Pointer      SegmentImagePointer;
    typedef typename SegmentImageType::ConstPointer SegmentImageConstPointer;
    typedef typename SegmentImageType::PixelType    SegmentImagePixelType;
    typedef typename SegmentImageType::RegionType   SegmentImageRegionType;
    typedef typename SegmentImageType::SizeType     SegmentImageSizeType;
    typedef typename SegmentImageSizeType::SizeValueType SegmentImageSizeValueType;
    typedef typename SegmentImageType::SpacingType  SegmentImageSpacingType;
    typedef typename SegmentImageType::IndexType    SegmentImageIndexType;
    typedef typename SegmentImageIndexType::IndexValueType
                                                    SegmentImageIndexValueType;
    typedef typename SegmentImageType::PointType    SegmentImagePointType;

    typedef unsigned int                                  LabelType;
    typedef ShapeLabelObject< LabelType, ImageDimension > ShapeLabelObjectType;
    typedef StatisticsLabelObject< LabelType, ImageDimension >
                                                          StatLabelObjectType;
    typedef LabelMap< ShapeLabelObjectType >              ShapeLabelMapType;
    typedef LabelMap< StatLabelObjectType >               StatLabelMapType;
    typedef LabelImageToShapeLabelMapFilter< SegmentImageType, ShapeLabelMapType >
                                                          ShapeConverterType;
    typedef LabelImageToStatisticsLabelMapFilter< SegmentImageType,
      FeatureImageType, StatLabelMapType >                StatConverterType;

    typedef ShapeLabelMapFilter< ShapeLabelMapType > ShapeFilterType;
    typedef ShapeRelabelImageFilter< SegmentImageType >     RelabelType;
    typedef LabelMapToLabelImageFilter< ShapeLabelMapType, SegmentImageType >
      MapToImageType;

    typedef LabelGeometryImageFilter< SegmentImageType > LabelGeometryFilterType;
    typedef typename LabelGeometryFilterType::Pointer LabelGeometryFilterPointer;

    typedef typename LabelGeometryFilterType::LabelPointType LabelPointType;

    typedef Vector< double, ImageDimension > VectorType;
    typedef Matrix< double, ImageDimension, ImageDimension > MatrixType;

    typedef vnl_matrix< double > vnlMatrixType;
    typedef vnl_vector< double > vnlVectorType;
    typedef typename std::list< SegmentImagePixelType > NeighborListType;

    itkGetConstMacro ( LargestCellRadius, double );
    itkSetMacro ( LargestCellRadius, double );
    itkGetConstMacro ( NumOfLabels, unsigned int );
    itkSetMacro ( NumOfLabels, unsigned int );

    itkSetMacro(ComputeSpatialMap, bool);
    itkGetConstReferenceMacro(ComputeSpatialMap, bool);
    itkBooleanMacro(ComputeSpatialMap);


    void SetRawImage ( FeatureImagePointer raw )
    {
      m_RawImg = raw;
      this->Modified();
    }

    void SetInput ( SegmentImagePointer seg )
    {
      m_Segmentation = seg;
      this->Modified();
    }

    std::vector< unsigned int > m_LabelLookup; // Size vectors
    std::vector< vnlVectorType > m_CenterOfGravity; // Centroid vectors
    std::vector< vnlVectorType > m_Centroid; // Centroid vectors
    std::vector< double > m_Mean; // Intensity means
    std::vector< double > m_Median; // Intensity median
    std::vector< unsigned int > m_Size; // Size vectors
    std::vector< double > m_PhysicalSize; // Physical size vectors
    std::vector< double > m_EquivalentSphericalRadius; // Physical size vectors
    std::vector< double > m_Elongation; // Elongation
    std::vector< vnlVectorType > m_BinaryMoments; // Principal moments
    std::vector< vnlVectorType > m_Moments; // Principal moments
    std::vector< vnlVectorType > m_EllipsoidDiameter; // Principal moments
    std::vector< vnlMatrixType > m_BinaryPrincipalAxis; // Principal axis
    std::vector< ImageRegionType > m_Region; // Region
    std::vector< NeighborListType > m_SpatialMap; // Spatial map based on nearest neighbor search
    std::vector< LabelPointType > m_OrientedBoundingBoxSize; // Oriented bounding box dimensions

    void GenerateData();


  protected:
    CellSegmentationStatistics();
    ~CellSegmentationStatistics() {}

    typename ShapeFilterType::Pointer m_Shape;
    typename RelabelType::Pointer     m_Relabel;

    FeatureImagePointer m_RawImg;
    SegmentImagePointer m_Segmentation;

    unsigned int m_NumOfLabels;
    double m_LargestCellRadius;
    bool m_ComputeSpatialMap;

  private:
    CellSegmentationStatistics ( Self& );   // intentionally not implemented
    void operator= ( const Self& );       // intentionally not implemented
  };

} /* namespace itk */

#include "itkCellSegmentationStatistics.txx"
#endif
