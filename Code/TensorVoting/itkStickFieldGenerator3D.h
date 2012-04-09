/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStickFieldGenerator3D.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:32 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkStickFieldGenerator3D_h
#define __itkStickFieldGenerator3D_h

#include "itkVoteFieldBase.h"
#include "itkStickFieldGenerator2D.h"
#include "itkMatrixLinearInterpolateImageFunction.h"

namespace itk
{

/** \class StickFieldGenerator3D
 * This calculator computes the stick field in 3D.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template <class TOutputImage>
class ITK_EXPORT StickFieldGenerator3D : public 
VoteFieldBase<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef StickFieldGenerator3D         Self;
  typedef VoteFieldBase<TOutputImage>   Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(StickFieldGenerator3D, VoteFieldBase);

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef TOutputImage                       ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::PixelType      PixelType;
  typedef typename ImageType::IndexType      IndexType;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename ImageType::PointType      PointType;
  typedef typename ImageType::SpacingType    SpacingType;

  typedef Matrix< double, ImageDimension-1, ImageDimension-1 > Matrix2DType;
  typedef Image< Matrix2DType, ImageDimension-1 > Image2DType;
  typedef typename Image2DType::Pointer Image2DPointer;
  typedef typename Image2DType::IndexType Index2DType;
  typedef typename Image2DType::PointType Point2DType;
  typedef typename Image2DType::RegionType Region2DType;
  typedef typename Image2DType::SpacingType Spacing2DType;

  typedef StickFieldGenerator2D< Image2DType > Stick2DGeneratorType;
  typedef typename Stick2DGeneratorType::Pointer Stick2DGeneratorPointer;
  typedef ImageRegionIteratorWithIndex< Image2DType > Iterator2DType;

  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef MatrixLinearInterpolateImageFunction< Image2DType, double > InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  /** Compute the stick field. */
  void ComputeStickField(void);

protected:
  StickFieldGenerator3D();
  virtual ~StickFieldGenerator3D() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  MatrixType RotationMatrixAboutZ( double theta );
  double FindTheta( double x, double y );

  InterpolatorPointer interpolator;

private:
  StickFieldGenerator3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStickFieldGenerator3D.txx"
#endif

#endif /* __itkStickFieldGenerator3D_h */
