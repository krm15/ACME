/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPlateFieldGenerator3D.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:32 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPlateFieldGenerator3D_h
#define __itkPlateFieldGenerator3D_h

#include "itkStickFieldGenerator3D.h"
#include "itkRotateFieldGenerator.h"
#include "itkGenerateRotationMatrixHelper.h"

namespace itk
{
/** \class PlateFieldGenerator3D
 * This calculator computes the ball field in 3D.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template < class TOutputImage >
class ITK_EXPORT PlateFieldGenerator3D : public 
RotateFieldGenerator<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef PlateFieldGenerator3D         Self;
  typedef RotateFieldGenerator<TOutputImage>  Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PlateFieldGenerator3D, RotateFieldGenerator);

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef typename Superclass::ImageType  ImageType;
  typedef typename Superclass::ImagePointer ImagePointer;
  typedef typename Superclass::SpacingType SpacingType;
  typedef typename Superclass::PointType  PointType;
  typedef typename Superclass::RegionType  RegionType;
  typedef typename Superclass::VectorType VectorType;
  typedef typename Superclass::MatrixType MatrixType;
  typedef typename Superclass::OrientedStickGeneratorType 
    OrientedTensorGeneratorType;
  typedef typename Superclass::OrientedStickGeneratorPointer
    OrientedTensorGeneratorPointer;

  typedef StickFieldGenerator3D<ImageType> StickGeneratorType;
  typedef typename StickGeneratorType::Pointer StickGeneratorPointer;

  typedef GenerateRotationMatrixHelper< ImageType >  
    RotationMatrixHelperType;

  /** Compute the plate field. */
  void ComputePlateField();

protected:
  PlateFieldGenerator3D();
  virtual ~PlateFieldGenerator3D() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void Initialize();

  RotationMatrixHelperType rMatrixHelper;

private:
  PlateFieldGenerator3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPlateFieldGenerator3D.txx"
#endif

#endif /* __itkPlateFieldGenerator3D_h */
