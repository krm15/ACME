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
#ifndef __itkTensorToSaliencyImageFilter_h
#define __itkTensorToSaliencyImageFilter_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkImageToImageFilter.h"
#include "itkSymmetricEigenAnalysis.h"

namespace itk
{

/** \class TensorToSaliencyImageFilter
 * This calculator computes the saliency image given a tensor image
 *
 * \ingroup Operators
 */
template <class TInputImage, class TOutputImage >
class ITK_EXPORT TensorToSaliencyImageFilter :
public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef TensorToSaliencyImageFilter          Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TensorToSaliencyImageFilter, ImageToImageFilter );

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef TInputImage                        ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::PointType      PointType;
  typedef typename ImageType::SpacingType    SpacingType;

  typedef TOutputImage                       OutputImageType;
  typedef typename OutputImageType::Pointer  OutputImagePointer;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension > MatrixType;
  typedef Image< MatrixType, ImageDimension > MatrixImageType;
  typedef typename MatrixImageType::Pointer MatrixImagePointer;
  typedef ImageRegionConstIterator< ImageType > ConstIteratorType;
  typedef ImageRegionIterator< OutputImageType > OutputIteratorType;

  typedef ImageRegionIterator< MatrixImageType > MatrixIteratorType;
	
  typedef SymmetricEigenAnalysis< MatrixType, VectorType, MatrixType > EigenCalculatorType;

  MatrixImagePointer GetEigenMatrix()
    {
    return m_EigenMatrix;
    }
	
  itkGetConstMacro ( ComputeEigenMatrix, bool );
  itkSetMacro ( ComputeEigenMatrix, bool );
	
protected:
  TensorToSaliencyImageFilter();
  virtual ~TensorToSaliencyImageFilter() {}

  void EnlargeOutputRequestedRegion(DataObject *output);
  void GenerateInputRequestedRegion();
  void BeforeThreadedGenerateData();
  void AfterThreadedGenerateData(){}
  void ThreadedGenerateData(const RegionType& windowRegion,
    ThreadIdType threadId);

  void PrintSelf(std::ostream& os, Indent indent) const;

  EigenCalculatorType m_EigenCalculator;
	
  bool m_ComputeEigenMatrix;
  MatrixImagePointer m_EigenMatrix;

private:
  TensorToSaliencyImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorToSaliencyImageFilter.txx"
#endif

#endif /* __itkTensorToSaliencyImageFilter_h */
