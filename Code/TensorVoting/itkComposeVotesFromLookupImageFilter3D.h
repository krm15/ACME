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
#ifndef __itkComposeVotesFromLookupImageFilter3D_h
#define __itkComposeVotesFromLookupImageFilter3D_h

#include "itkComposeVotesFromLookupImageFilter.h"
#include <itkAzimuthElevationToCartesianTransform.h>
#include "itkGenerateRotationMatrixHelper.h"

namespace itk
{

/** \class ComposeVotesFromLookupImageFilter3D
 * This class composes votes given a lookup image of voxels with same orientation.
 * It is templated over the input image type and output image type.
 *
 * \ingroup Operators
 */
template <class TInputImage >
class ITK_EXPORT ComposeVotesFromLookupImageFilter3D :
public ComposeVotesFromLookupImageFilter< TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef ComposeVotesFromLookupImageFilter3D Self;
  typedef ComposeVotesFromLookupImageFilter< TInputImage > Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  // This is purposely not provided since this is an abstract class.
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComposeVotesFromLookupImageFilter3D,
                ComposeVotesFromLookupImageFilter );

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
  typedef typename SizeType::SizeValueType        SizeValueType;
  typedef typename InputImageType::PointType      PointType;
  typedef typename InputImageType::SpacingType    SpacingType;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension > MatrixType;
  typedef Image< MatrixType, ImageDimension >     OutputImageType;
  typedef AzimuthElevationToCartesianTransform< double, ImageDimension > CoordinateTransformType;
  typedef typename CoordinateTransformType::Pointer CoordinateTransformPointer;
  typedef GenerateRotationMatrixHelper< OutputImageType > RotationMatrixHelperType;

protected:
  ComposeVotesFromLookupImageFilter3D();
  virtual ~ComposeVotesFromLookupImageFilter3D() {}
  void ComputeRotationMatrix( MatrixType& R,  PointType& pt_sph );

  RotationMatrixHelperType rMatrixHelper;
  CoordinateTransformPointer m_Transform;

private:
  ComposeVotesFromLookupImageFilter3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComposeVotesFromLookupImageFilter3D.txx"
#endif

#endif /* __itkComposeVotesFromLookupImageFilter3D_h */
