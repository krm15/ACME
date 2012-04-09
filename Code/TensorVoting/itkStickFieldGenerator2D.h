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
#ifndef __itkStickFieldGenerator2D_h
#define __itkStickFieldGenerator2D_h

#include "itkVoteFieldBase.h"

namespace itk
{

/** \class StickFieldGenerator2D
 * This calculator computes the stick field in 2D.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template <class TOutputImage>
class ITK_EXPORT StickFieldGenerator2D : public 
VoteFieldBase<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef StickFieldGenerator2D         Self;
  typedef VoteFieldBase<TOutputImage>   Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(StickFieldGenerator2D, VoteFieldBase);

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

  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;

  /** Compute the stick field. */
  void ComputeStickField(void);

protected:
  StickFieldGenerator2D();
  virtual ~StickFieldGenerator2D() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  StickFieldGenerator2D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStickFieldGenerator2D.txx"
#endif

#endif /* __itkStickFieldGenerator2D_h */
