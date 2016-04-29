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
#ifndef __itkGPUTensorWithOrientation_h
#define __itkGPUTensorWithOrientation_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkImageToImageFilter.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "itkTensorWithOrientation.h"
#include "itkGPUImageToImageFilter.h"
#include "itkObjectFactoryBase.h"
#include "itkOpenCLUtil.h"

#include "itkGPUDataManager.h"
#include "itkGPUKernelManager.h"

// New class written for matrix interpolation
#include "itkMatrixLinearInterpolateImageFunction.h"

namespace itk
{

/** \class StickFieldGeneratorWithOrientation2D
 * This calculator reorients the voting field given a rotation matrix.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */

template< typename TInputImage, typename TOutputImage = TInputImage >
class GPUTensorWithOrientation :
  public GPUImageToImageFilter< TInputImage, TOutputImage, TensorWithOrientation< TInputImage, TOutputImage > >
{
public:
  /** Standard class typedefs. */
  typedef GPUTensorWithOrientation Self;
  typedef TensorWithOrientation< TInputImage, TOutputImage >  CPUSuperclass;
  typedef GPUImageToImageFilter< TInputImage, TOutputImage, CPUSuperclass >  GPUSuperclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information support. */
  itkTypeMacro(GPUTensorWithOrientation, GPUImageToImageFilter );

  itkStaticConstMacro ( ImageDimension, unsigned int, TInputImage::ImageDimension );

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   ImageConstPointer;
  typedef typename InputImageType::PixelType      PixelType;
  typedef typename InputImageType::IndexType      IndexType;
  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType        SizeValueType;
  typedef typename InputImageType::PointType      PointType;
  typedef typename InputImageType::SpacingType    SpacingType;

  /** Type definition for the output image. */
  typedef TOutputImage                       OutputImageType;
  typedef typename OutputImageType::Pointer  OutputImagePointer;

  /** GPU type definition */
  typedef typename GPUTraits< TInputImage >::Type  GPUInputImage;
  typedef typename GPUInputImage::Pointer               GPUInputImagePointer;
  typedef typename GPUTraits< TOutputImage >::Type GPUOutputImage;
  typedef typename GPUOutputImage::Pointer              GPUOutputImagePointer;
  typedef GPUImageDataManager<GPUInputImage>            GPUInputManagerType;
  typedef GPUImageDataManager<GPUOutputImage>           GPUOutputManagerType;

  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef Vector< double, ImageDimension > VectorType;
  typedef ImageRegionConstIterator< InputImageType > InputIteratorType;
  typedef ImageRegionIterator< OutputImageType > OutputIteratorType;
  typedef MatrixLinearInterpolateImageFunction< InputImageType, double > InterpolatorType;
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
  GPUTensorWithOrientation();

  virtual ~GPUTensorWithOrientation() {}

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  /** Standard GPU pipeline method. */
  void GPUGenerateData() ITK_OVERRIDE;

  SpacingType m_OutputSpacing;
  RegionType m_OutputRegion;
  bool m_RotationMatrixDefined;
  MatrixType m_RotationMatrix;
  PointType m_Center;
  InterpolatorPointer interpolator;

  typename GPUKernelManager::Pointer m_GPUKernelManager;

  /** GPU kernel handle for GPUInterpolate() */
  int m_InterpolateGPUKernelHandle;

private:
  GPUTensorWithOrientation(const Self&) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGPUTensorWithOrientation.txx"
#endif

#endif /* __itkGPUTensorWithOrientation_h */
