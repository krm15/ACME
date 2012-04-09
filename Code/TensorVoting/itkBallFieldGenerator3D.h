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
#ifndef __itkBallFieldGenerator3D_h
#define __itkBallFieldGenerator3D_h

#include "itkStickFieldGenerator3D.h"
#include "itkRotateFieldGenerator.h"
#include "itkGenerateRotationMatrixHelper.h"

namespace itk
{
/** \class BallFieldGenerator3D
 * This calculator computes the ball field in 3D.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template <class TOutputImage>
class ITK_EXPORT BallFieldGenerator3D : public 
RotateFieldGenerator<TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef BallFieldGenerator3D          Self;
  typedef RotateFieldGenerator<TOutputImage>  Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BallFieldGenerator2D, RotateFieldGenerator);

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef typename Superclass::ImageType  ImageType;
  typedef typename Superclass::SpacingType SpacingType;
  typedef typename Superclass::PointType  PointType;
  typedef typename Superclass::RegionType  RegionType;
  typedef typename Superclass::VectorType VectorType;
  typedef typename Superclass::MatrixType MatrixType;

  typedef StickFieldGenerator2D<ImageType> StickGeneratorType;
  typedef typename StickGeneratorType::Pointer StickGeneratorPointer;

  typedef GenerateRotationMatrixHelper< ImageType >  
    RotationMatrixHelperType;

  /** Compute the stick field. */
  void ComputeBallField(void);

protected:
  BallFieldGenerator3D();
  virtual ~BallFieldGenerator3D() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  void Initialize();

  RotationMatrixHelperType rMatrixHelper;

private:
  BallFieldGenerator3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBallFieldGenerator3D.txx"
#endif

#endif /* __itkBallFieldGenerator3D_h */
