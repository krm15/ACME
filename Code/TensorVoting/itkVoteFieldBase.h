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
#ifndef __itkVoteFieldBase_h
#define __itkVoteFieldBase_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkTensorWithOrientation.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class VoteFieldBase
 * This calculator computes the ball field in 2D.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template <class TOutputImage >
class ITK_EXPORT VoteFieldBase : public Object
{
public:
  /** Standard class typedefs. */
  typedef VoteFieldBase                 Self;
  typedef Object                        Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VoteFieldBase, Object);

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

  typedef vnl_vector< double > vnlVectorType;
  typedef vnl_matrix< double > vnlMatrixType;
  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;

  itkSetMacro( Sigma, double );
  itkGetMacro( Sigma, double );
  itkSetMacro( Spacing, SpacingType );
  itkGetMacro( Spacing, SpacingType );

  ImagePointer GetOutputImage()
  {
    return m_Output;
  }

protected:
  VoteFieldBase();
  virtual ~VoteFieldBase() {};
  void AllocateOutput();
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void Initialize(){}
  virtual void Postprocess( ImagePointer image, double& count);

  SpacingType m_Spacing;
  PointType m_Origin;
  double m_Sigma;
  ImagePointer m_Output;

private:
  VoteFieldBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVoteFieldBase.txx"
#endif

#endif /* __itkVoteFieldBase_h */
