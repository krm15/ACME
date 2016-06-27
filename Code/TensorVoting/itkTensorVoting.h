/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 585 $  // Revision of last commit
  Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#ifndef __itkTensorVoting_h
#define __itkTensorVoting_h

#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkImageToImageFilter.h"
#include "itkTensorWithOrientation.h"
#include "itkNumericTraits.h"
#include <vector>
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkThreadSafeMersenneTwisterRandomVariateGenerator.h"
#include "itkSymmetricEigenAnalysis.h"

// #include "itkImageFileWriter.h"

namespace itk
{

/** \class TensorVoting
 * This abstract class performs tensor voting in nD.
 * It is templated over the input image type.
 *
 * \ingroup Operators
 */
template <class TInputImage >
class ITK_EXPORT TensorVoting :
public ImageToImageFilter< TInputImage, TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef TensorVoting              Self;
  typedef ImageToImageFilter<TInputImage, TInputImage >
                                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  // This is purposely not provided since this is an abstract class.
  //   itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TensorVoting, ImageToImageFilter );

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

  /** Type definition for the output image. */
  typedef Image< bool, ImageDimension >      TokenImageType;
  typedef typename TokenImageType::Pointer   TokenImagePointer;
  typedef ImageRegionIteratorWithIndex< TokenImageType > TokenIteratorType;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionConstIteratorWithIndex< InputImageType > ConstIteratorType;
  typedef ImageRegionIteratorWithIndex< InputImageType > IteratorType;

  typedef Image< VectorType, ImageDimension > VectorImageType;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  typedef ImageRegionIterator< VectorImageType > VectorIteratorType;

  typedef std::list< IndexType > IdType;
  typedef Vector< IdType, ImageDimension > VectorInternalType;
  typedef Image<VectorInternalType, ImageDimension> InternalImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;
  typedef ImageRegionIteratorWithIndex< InternalImageType > InternalIteratorType;

  typedef Image< double, ImageDimension > DoubleImageType;
  typedef typename DoubleImageType::Pointer DoubleImagePointer;
  typedef VectorIndexSelectionCastImageFilter< VectorImageType, DoubleImageType > IndexFilterType;
  typedef ImageRegionIteratorWithIndex< DoubleImageType > DoubleIteratorType;

  typedef SymmetricEigenAnalysis< MatrixType, VectorType, MatrixType >
    EigenCalculatorType;

  itkSetMacro( Sigma, double );
  itkGetMacro( Sigma, double );
  itkSetMacro( UseSparseVoting, bool );
  itkGetMacro( UseSparseVoting, bool );
  itkSetMacro( LowMemoryFilter, bool );
  itkGetMacro( LowMemoryFilter, bool );

  itkSetObjectMacro( TokenImage, TokenImageType );
  itkSetObjectMacro( SaliencyImage, VectorImageType );
  itkSetObjectMacro( EigenMatrixImage, InputImageType );

protected:
  TensorVoting();
  virtual ~TensorVoting() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();

  void OverlapRegion( InputImagePointer A, InputImagePointer B, RegionType& rA, RegionType& rB );
  double ComputeTheta( VectorType& u );
  void InitializeLookupImages();

  virtual void InitializeVotingFields() = 0;
  virtual void ComputeLookup() = 0;
  virtual void IntegrateVotes() = 0;

  double m_Sigma;
  bool m_UseSparseVoting;
  bool m_LowMemoryFilter;
  RegionType m_Region;

  TokenImagePointer m_TokenImage;
  InternalImagePointer m_Lookup;
  VectorImagePointer m_SaliencyImage;
  InputImagePointer m_EigenMatrixImage;
  InputImagePointer m_Output;
  std::vector< InputImagePointer > m_VotingField;
  EigenCalculatorType m_EigenCalculator;

private:
  TensorVoting(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorVoting.txx"
#endif

#endif /* __itkTensorVoting_h */
