/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 663 $  // Revision of last commit
  Date: $Date: 2009-09-15 18:47:59 -0400 (Tue, 15 Sep 2009) $  // Date of last commit
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

#ifndef __itkLabelMapOverlapMeasures_h_
#define __itkLabelMapOverlapMeasures_h_

#include "itkCastImageFilter.h"
#include "itkLabelObject.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "itkImageRegionConstIterator.h"
#include "itkVector.h"
#include <map>
#include <list>

namespace itk {

template < class TInputImage >
class ITK_EXPORT LabelMapOverlapMeasures : public ImageToImageFilter<
  TInputImage, TInputImage >
{
public:

  typedef LabelMapOverlapMeasures Self;
  typedef ImageToImageFilter< TInputImage,TInputImage > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( LabelMapOverlapMeasures, ImageToImageFilter );

  /** Display */
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** Input image typedefs */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PixelType      InputPixelType;
  typedef typename InputImageType::RegionType     InputRegionType;
  typedef typename InputImageType::SizeType       InputSizeType;
  typedef typename InputSizeType::SizeValueType   InputSizeValueType;
  typedef typename InputImageType::SpacingType    InputSpacingType;
  typedef typename InputImageType::IndexType      InputIndexType;
  typedef typename InputIndexType::IndexValueType InputIndexValueType;
  typedef typename InputImageType::PointType      InputPointType;

  typedef BinaryThresholdImageFilter< InputImageType, InputImageType > BinaryThresholdFilterType;
  typedef typename BinaryThresholdFilterType::Pointer BinaryThresholdFilterPointer;
  typedef HausdorffDistanceImageFilter<InputImageType, InputImageType> HausdorffDistanceFilterType;
  typedef typename HausdorffDistanceFilterType::Pointer HausdorffDistanceFilterPointer;

  typedef ShapeLabelObject< InputPixelType, ImageDimension > ShapeLabelObjectType;
  typedef LabelMap< ShapeLabelObjectType > ShapeLabelMapType;
  typedef typename ShapeLabelMapType::Pointer ShapeLabelMapPointer;
  typedef LabelImageToShapeLabelMapFilter< InputImageType, ShapeLabelMapType >
    ShapeConverterType;
  typedef typename ShapeConverterType::Pointer ShapeConverterPointer;

  typedef ImageRegionConstIterator< InputImageType > ImageConstIteratorType;

  typedef vnl_matrix< double > DoubleMatrixType;
  typedef vnl_matrix< int > IntMatrixType;
  typedef vnl_vector< unsigned int > VectorType;
  typedef vnl_vector< double > DoubleVectorType;

  typedef itk::Vector< unsigned int, 2 > PairVectorType;
  typedef itk::Vector< double, ImageDimension > MeasurementVectorType;

  typedef std::map< double, PairVectorType, std::less< double > > DiceMapType; // std::greater< double >

  itkGetMacro( NumberOfMatches, unsigned int );
  itkSetMacro( NumberOfMatches, unsigned int );
  itkGetMacro( NumberOfSegLabels, unsigned int );
  itkSetMacro( NumberOfSegLabels, unsigned int );
  itkGetMacro( NumberOfGroundTruthLabels, unsigned int );
  itkSetMacro( NumberOfGroundTruthLabels, unsigned int );
  itkGetMacro( NumberOfOverSegmentations, unsigned int );
  itkSetMacro( NumberOfOverSegmentations, unsigned int );
  itkGetMacro( NumberOfUnderSegmentations, unsigned int );
  itkSetMacro( NumberOfUnderSegmentations, unsigned int );
  itkGetMacro( AverageValueOfDiceMetric, double );
  itkSetMacro( AverageValueOfDiceMetric, double );
  itkGetMacro( AverageValueOfL2Metric, double );
  itkSetMacro( AverageValueOfL2Metric, double );

  // Set these thresholds as inputs
  itkGetMacro( ConfusionThreshold, double );
  itkSetMacro( ConfusionThreshold, double );
  itkGetMacro( ThresholdForUndersegmentation, double );
  itkSetMacro( ThresholdForUndersegmentation, double );

  std::map< unsigned int, unsigned int > m_GroundTruthLabels;
  std::map< unsigned int, unsigned int > m_SegLabels;
  IntMatrixType m_OverlapArea;
  VectorType m_GroundTruthArea, m_SegmentedArea;

protected:

  LabelMapOverlapMeasures();
  ~LabelMapOverlapMeasures(){}
  void GenerateData();

  void IdentifyNumberOfCells();
  double EvaluateHaussdorf( unsigned int iR, unsigned int iC );
  void CreateConfusionMatrix();
  void AssignMatchingLabels();
  void IdentifyMismatch();

  unsigned int m_NumberOfGroundTruthLabels;
  unsigned int m_NumberOfSegLabels;
  unsigned int m_NumberOfMatches;
  unsigned int m_NumberOfOverSegmentations;
  unsigned int m_NumberOfUnderSegmentations;
  double m_AverageValueOfDiceMetric;
  double m_AverageValueOfL2Metric;
  double m_ConfusionThreshold;
  double m_ThresholdForUndersegmentation;

  DiceMapType m_Dice; // Dice

private:

  LabelMapOverlapMeasures(Self&);   // intentionally not implemented
  void operator=(const Self&);   // intentionally not implemented
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelMapOverlapMeasures.txx"
#endif

#endif
