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

#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkImageToImageFilter.h"
#include "itkTensorWithOrientation.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkNumericTraits.h"
#include <vector>


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
  typedef TInputImage                        ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::PixelType      PixelType;
  typedef typename ImageType::IndexType      IndexType;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename ImageType::PointType      PointType;
  typedef typename ImageType::SpacingType    SpacingType;

  typedef TInputImage                        InputImageType;
  typedef typename ImageType::Pointer        InputImagePointer;
  typedef typename ImageType::RegionType     OutputImageRegionType;

  /** Type definition for the output image. */
  typedef Image< bool, ImageDimension >      TokenImageType;
  typedef typename TokenImageType::Pointer   TokenImagePointer;
  typedef ImageRegionIteratorWithIndex< TokenImageType > 
    TokenIteratorType;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionConstIteratorWithIndex< ImageType > 
    ConstIteratorType;
  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;

	typedef TensorWithOrientation< ImageType > OrientedTensorGeneratorType;
  typedef typename OrientedTensorGeneratorType::Pointer 
    OrientedTensorGeneratorPointer;
  typedef SymmetricEigenAnalysis< MatrixType, VectorType, MatrixType > 
    EigenCalculatorType;


  itkSetMacro( Sigma, double );
  itkGetMacro( Sigma, double );
  itkSetMacro( UseSparseVoting, bool );
  itkGetMacro( UseSparseVoting, bool );

  void SetTokenImage( TokenImagePointer token )
  {
    m_TokenImage = token;
  }

  TokenImagePointer GetTokenImage()
  {
    return m_TokenImage;
  }

protected:
  TensorVoting();
  virtual ~TensorVoting() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ComputeOrientedField( double saliency, PointType& iCenter, 
    MatrixType& R, unsigned int i, int threadId );
  void ComputeVote( ImagePointer field, double saliency, int threadId );
  void OverlapRegion( ImagePointer A, ImagePointer B, 
    RegionType& rA, RegionType& rB );

  void EnlargeOutputRequestedRegion(DataObject *output);
  void GenerateInputRequestedRegion();
  void BeforeThreadedGenerateData();
  void AfterThreadedGenerateData();
  void PadImage(const OutputImageRegionType& windowRegion, int threadId);
  void ThreadedGenerateData(const OutputImageRegionType& windowRegion, 
    ThreadIdType threadId);

  virtual void InitializeVotingFields(void);
  virtual void ComputeThetaAndRotationAxis( VectorType& u, 
    MatrixType& R ) = 0;
  virtual void ComputeTensorVoting( IndexType& index, 
    PixelType& p, int threadId ) = 0;

  double m_Sigma;
  bool m_UseSparseVoting;
  RegionType m_Region;

  TokenImagePointer m_TokenImage;
  ImagePointer m_Output;
  EigenCalculatorType m_EigenCalculator;
  std::vector< ImagePointer > m_VotingField;
  std::vector<ImagePointer> m_ThreadImage;

private:
  TensorVoting(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorVoting.txx"
#endif

#endif /* __itkTensorVoting_h */
