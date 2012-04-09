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
#ifndef __itkRotateFieldGenerator_h
#define __itkRotateFieldGenerator_h

#include "itkTensorWithOrientation.h"
#include "itkVoteFieldBase.h"


namespace itk
{

/** \class RotateFieldGenerator
 * This calculator computes the ball field in 2D.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template <class TOutputImage >
class ITK_EXPORT RotateFieldGenerator : public 
VoteFieldBase<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef RotateFieldGenerator          Self;
  typedef VoteFieldBase<TOutputImage>   Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RotateFieldGenerator, VoteFieldBase);

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef TOutputImage                       ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SpacingType    SpacingType;

  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;

	typedef TensorWithOrientation< ImageType > OrientedStickGeneratorType;
  typedef typename OrientedStickGeneratorType::Pointer 
    OrientedStickGeneratorPointer;

  void SetInput( ImagePointer v )
  {
    m_VotingField = v;
  }

protected:
  RotateFieldGenerator();
  virtual ~RotateFieldGenerator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

//   virtual void Initialize(){}
  virtual void Update(ImagePointer image, MatrixType R, 
    SpacingType space, RegionType region, double weight );

  ImagePointer m_VotingField;

private:
  RotateFieldGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRotateFieldGenerator.txx"
#endif

#endif /* __itkRotateFieldGenerator_h */
