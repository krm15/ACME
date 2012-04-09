/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStickFieldGenerator2D.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:32 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorWithOrientation_h
#define __itkTensorWithOrientation_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkImageToImageFilter.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

// New class written for matrix interpolation
#include "itkMatrixLinearInterpolateImageFunction.h"

namespace itk
{

/** \class StickFieldGeneratorWithOrientation2D
 * This calculator computes the stick field in 2D for a given angle.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template < class TInputImage, class TOutputImage = TInputImage >
class ITK_EXPORT TensorWithOrientation :
	public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef TensorWithOrientation Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage >  Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorWithOrientation, ImageToImageFilter );

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

	typedef TInputImage                        ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::PixelType      PixelType;
  typedef typename ImageType::IndexType      IndexType;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename ImageType::PointType      PointType;
  typedef typename ImageType::SpacingType    SpacingType;

  /** Type definition for the output image. */
  typedef TOutputImage                       OutputImageType;
  typedef typename ImageType::Pointer        OutputImagePointer;

  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
	typedef Vector< double, ImageDimension > VectorType;
  typedef ImageRegionConstIterator< ImageType > InputIteratorType;
  typedef ImageRegionIterator< OutputImageType > OutputIteratorType;
  typedef MatrixLinearInterpolateImageFunction< ImageType, double > 
    InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  typedef vnl_vector< double > vnlVectorType;
  typedef vnl_matrix< double > vnlMatrixType;

  void SetRotationMatrix( MatrixType m )
  {
    m_RotationMatrix = m;
    m_RotationMatrixDefined = true;
  }
  itkGetMacro( RotationMatrix, MatrixType );
  itkSetMacro( Center, PointType );
  itkGetMacro( Center, PointType );
  itkSetMacro( OutputSpacing, SpacingType );
  itkGetMacro( OutputSpacing, SpacingType );
  itkSetMacro( OutputRegion, RegionType );
  itkGetMacro( OutputRegion, RegionType );

protected:
  TensorWithOrientation();
  virtual ~TensorWithOrientation() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
	void GenerateData();

  SpacingType m_OutputSpacing;
  RegionType m_OutputRegion;
  bool m_RotationMatrixDefined;
  MatrixType m_RotationMatrix;
  PointType m_Center;
  InterpolatorPointer interpolator;

private:
  TensorWithOrientation(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorWithOrientation.txx"
#endif

#endif /* __itkTensorWithOrientation_h */
