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

#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkImageToImageFilter.h"
#include "itkTensorWithOrientation.h"
#include "itkNumericTraits.h"
#include <vector>
#include "itkTensorToSaliencyImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBallFieldGenerator2D.h"
#include "itkStickFieldGenerator2D.h"
#include "itkComposeVotesFromLookupImageFilter.h"
#include "itkThreadSafeMersenneTwisterRandomVariateGenerator.h"

// #include "itkImageFileWriter.h"

namespace itk
{

/** \class TensorVoting2D
 * This calculator performs tensor voting in 3D.
 * It is templated over the input image type.
 *
 * \ingroup Operators
 */
template <class TInputImage >
class ITK_EXPORT TensorVoting2D :
       public ImageToImageFilter< TInputImage, TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef TensorVoting2D              Self;
  typedef ImageToImageFilter<TInputImage, TInputImage > Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

/** Method for creation through the object factory. */
itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TensorVoting2D, ImageToImageFilter );

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Type definition for the output image. */
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

  typedef StickFieldGenerator2D< InputImageType >   StickGeneratorType;
  typedef typename StickGeneratorType::Pointer      StickGeneratorPointer;
  typedef BallFieldGenerator2D< InputImageType >    BallGeneratorType;
  typedef typename BallGeneratorType::Pointer       BallGeneratorPointer;

  typedef GenerateRotationMatrixHelper< InputImageType > RotationMatrixHelperType;

  /** Type definition for the output image. */
  typedef Image< bool, ImageDimension >      TokenImageType;
  typedef typename TokenImageType::Pointer   TokenImagePointer;
  typedef ImageRegionIteratorWithIndex< TokenImageType > TokenIteratorType;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionConstIteratorWithIndex< InputImageType > ConstIteratorType;
  typedef ImageRegionIteratorWithIndex< InputImageType > IteratorType;

  typedef std::list< IndexType > IdType;
  typedef Image<IdType, ImageDimension> InternalImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;
  typedef ImageRegionIteratorWithIndex< InternalImageType > InternalIteratorType;
  
  typedef Image< double, ImageDimension > DoubleImageType;
  typedef typename DoubleImageType::Pointer DoubleImagePointer;
  typedef Image< VectorType, ImageDimension > VectorImageType;
  typedef TensorToSaliencyImageFilter< InputImageType, VectorImageType > SaliencyFilterType;
  typedef VectorIndexSelectionCastImageFilter< VectorImageType, DoubleImageType > IndexFilterType;
  typedef ImageRegionIteratorWithIndex< DoubleImageType > DoubleIteratorType;
  
  typedef ComposeVotesFromLookupImageFilter< InternalImageType >
    ComposeVotesFilterType;
  typedef typename ComposeVotesFilterType::Pointer ComposeVotesFilterPointer;

//   typedef Image< unsigned short, ImageDimension > OutputImageType;
//   typedef typename OutputImageType::Pointer OutputImagePointer;
//   typedef ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;
//   typedef ImageFileWriter< OutputImageType > WriterType;
//   typedef typename WriterType::Pointer WriterPointer;

  typedef Statistics::ThreadSafeMersenneTwisterRandomVariateGenerator
    RandomGeneratorType;
  typedef typename RandomGeneratorType::Pointer RandomGeneratorPointer;

  itkSetMacro( Sigma, double );
  itkGetMacro( Sigma, double );
  itkSetMacro( UseSparseVoting, bool );
  itkGetMacro( UseSparseVoting, bool );

  void SetTokenImage( TokenImagePointer token )
  {
    m_TokenImage = token;
  }

  void SetStickSaliencyImage( DoubleImagePointer saliency )
  {
    m_StickSaliencyImage = saliency;
  }

  void SetBallSaliencyImage( DoubleImagePointer saliency )
  {
    m_BallSaliencyImage = saliency;
  }

  void SetEigenMatrixImage( InputImagePointer eigen )
  {
    m_EigenMatrixImage = eigen;
  }

  TokenImagePointer GetTokenImage()
  {
    return m_TokenImage;
  }


protected:
  TensorVoting2D();
  virtual ~TensorVoting2D() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void ComputeLookup();
  void OverlapRegion( InputImagePointer A, InputImagePointer B, RegionType& rA, RegionType& rB );

  void InitializeVotingFields();
  void GenerateData();
  void InitializeLookupImages();
  double ComputeTheta( VectorType& u );

  double m_Sigma;
  bool m_UseSparseVoting;
  RegionType m_Region;

  TokenImagePointer m_TokenImage;
  InternalImagePointer m_LookupStick;
  InternalImagePointer m_LookupBall;

  DoubleImagePointer m_StickSaliencyImage;
  DoubleImagePointer m_BallSaliencyImage;

  InputImagePointer m_EigenMatrixImage;
  InputImagePointer m_OrientedVotingField;
  InputImagePointer m_Output;
  std::vector< InputImagePointer > m_VotingField;

private:
  TensorVoting2D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorVoting2D.txx"
#endif

#endif /* __itkTensorVoting2D_h */
