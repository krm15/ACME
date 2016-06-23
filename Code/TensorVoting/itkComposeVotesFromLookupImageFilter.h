/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 585 $  // Revision of last commit
  Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
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
#ifndef __itkComposeVotesFromLookupImageFilter_h
#define __itkComposeVotesFromLookupImageFilter_h

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkInPlaceImageFilter.h"
#include "itkNumericTraits.h"
#include <vector>
#include "itkAddNonScalarImageFilter.txx"

#include <itkAzimuthElevationToCartesianTransform.h>
#include "itkGenerateRotationMatrixHelper.h"

namespace itk
{

/** \class ComposeVotesFromLookupImageFilter
 * This class composes votes given a lookup image of voxels with same orientation.
 * It is templated over the input image type and output image type.
 *
 * \ingroup Operators
 */
template <class TInputImage >
class ITK_EXPORT ComposeVotesFromLookupImageFilter :
public InPlaceImageFilter< TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef ComposeVotesFromLookupImageFilter Self;
  typedef InPlaceImageFilter< TInputImage > Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  // This is purposely not provided since this is an abstract class.
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComposeVotesFromLookupImageFilter, InPlaceImageFilter );

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PixelType      PixelType;
  typedef typename InputImageType::IndexType      IndexType;
  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename InputImageType::PointType      PointType;
  typedef typename InputImageType::SpacingType    SpacingType;

  typedef ImageRegionConstIterator< InputImageType > ConstInputIteratorType;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension > MatrixType;
  typedef std::list< IndexType > IdType;
  typedef Image< MatrixType, ImageDimension >     OutputImageType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef ImageRegionIterator< OutputImageType >  OutputIteratorType;

  typedef GenerateRotationMatrixHelper< OutputImageType > RotationMatrixHelperType;
  typedef TensorWithOrientation< OutputImageType > OrientedTensorGeneratorType;
  typedef typename OrientedTensorGeneratorType::Pointer OrientedTensorGeneratorPointer;
  typedef Image< double, ImageDimension > DoubleImageType;
  typedef typename DoubleImageType::Pointer DoubleImagePointer;
  
  typedef AddNonScalarImageFilter< OutputImageType > AddNonScalarFilterType;
  typedef typename AddNonScalarFilterType::Pointer AddNonScalarFilterPointer;

  typedef AzimuthElevationToCartesianTransform< double, ImageDimension > CoordinateTransformType;
  typedef typename CoordinateTransformType::Pointer CoordinateTransformPointer;
 
    
  void SetVotingField( OutputImagePointer votingField )
  {
    m_VotingField = votingField;
  }

  void SetSaliencyImage( DoubleImagePointer saliencyImage )
  {
    m_SaliencyImage = saliencyImage;
  }

  void SetOutputImage( OutputImagePointer output )
  {
    m_Output = output;
  }



protected:
  ComposeVotesFromLookupImageFilter();
  virtual ~ComposeVotesFromLookupImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void EnlargeOutputRequestedRegion(DataObject *output);
  void GenerateInputRequestedRegion();
  void BeforeThreadedGenerateData();
  void AfterThreadedGenerateData();
  void ThreadedGenerateData(const RegionType& windowRegion, ThreadIdType threadId);

  void OverlapRegion( OutputImagePointer A, OutputImagePointer B, RegionType& rA, RegionType& rB );
  void ComputeVote( double saliency, unsigned int threadId, OutputImagePointer orientedVotingField );

  RotationMatrixHelperType rMatrixHelper;

  OutputImagePointer m_VotingField;
  DoubleImagePointer m_SaliencyImage;
  OutputImagePointer m_Output;
  std::vector<OutputImagePointer> m_ThreadImage;
  unsigned int m_ValidThreads;

private:
  ComposeVotesFromLookupImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComposeVotesFromLookupImageFilter.txx"
#endif

#endif /* __itkComposeVotesFromLookupImageFilter_h */
