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
#ifndef __itkTensorVoting2D_h
#define __itkTensorVoting2D_h

#include "itkTensorVoting.h"
#include "itkBallFieldGenerator2D.h"
#include "itkStickFieldGenerator2D.h"

namespace itk
{

/** \class TensorVoting2D
 * This calculator performs tensor voting in 2D.
 * It is templated over the input image type.
 *
 * \ingroup Operators
 */
template <class TInputImage >
class ITK_EXPORT TensorVoting2D : public TensorVoting< TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef TensorVoting2D              Self;
  typedef TensorVoting< TInputImage > Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TensorVoting2D, TensorVoting );

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef typename Superclass::ImageType     ImageType;
  typedef typename Superclass::ImageConstPointer 
    ImageConstPointer;
  typedef typename Superclass::SpacingType SpacingType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::IndexType      IndexType;
  typedef typename Superclass::RegionType     RegionType;
  typedef typename Superclass::SizeType       SizeType;
  typedef typename Superclass::SizeValueType   SizeValueType;
  typedef typename Superclass::PointType      PointType;

  typedef typename Superclass::VectorType VectorType;
  typedef typename Superclass::MatrixType MatrixType;

  typedef BallFieldGenerator2D< ImageType > BallGeneratorType;
  typedef typename BallGeneratorType::Pointer BallGeneratorPointer;
	typedef StickFieldGenerator2D< ImageType > StickGeneratorType;
	typedef typename StickGeneratorType::Pointer StickGeneratorPointer;

protected:
  TensorVoting2D();
  virtual ~TensorVoting2D() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void InitializeVotingFields(void);
  virtual void ComputeThetaAndRotationAxis( VectorType& u, MatrixType& R );
  virtual void ComputeTensorVoting( IndexType& index, PixelType& p, 
    int threadId );

private:
  TensorVoting2D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorVoting2D.txx"
#endif

#endif /* __itkTensorVoting2D_h */
