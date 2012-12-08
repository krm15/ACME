/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMembraneScaleSelectionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/06/12 22:59:15 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMembraneScaleSelectionImageFilter_h
#define __itkMembraneScaleSelectionImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
/** \class MembraneScaleSelectionImageFilter
 * \brief A filter to enhance 3D vascular structures
 *
 * The plate measure is based on the analysis of the the Hessian
 * eigen system. The filter takes an
 * image of a Hessian pixels ( SymmetricSecondRankTensor pixels ) and
 * produces an enhanced image. The Hessian input image can be produced using
 * itkHessianSmoothedRecursiveGaussianImageFilter.
 *
 *
 * \par References
 *  Mosaliganti, K, F. Janoos, Gelas, A, Machiraju, R and Megason, S (2009). Membrane Enhancing
 *  Diffusion: A Scale Space Representation of Membrane Structures. ISBI 2010
 *
 * \sa MultiScaleMembraneScaleSelectionImageFilter
 * \sa Hessian3DToMembranenessMeasureImageFilter
 * \sa HessianSmoothedRecursiveGaussianImageFilter
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 *
 * \ingroup IntensityImageFilters TensorObjects
 *
 */

template < typename  TImage, typename TOutputImage >
class ITK_EXPORT MembraneScaleSelectionImageFilter : public
ImageToImageFilter< TImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MembraneScaleSelectionImageFilter Self;

  typedef ImageToImageFilter< TImage, TOutputImage > Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef TImage  ImageType;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::PixelType   PixelType;
  typedef typename ImageType::RegionType  RegionType;

  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;

  /** Image dimension = 3. */
  itkStaticConstMacro( ImageDimension, unsigned int, ImageType::ImageDimension);

  itkStaticConstMacro(PixelDimension, unsigned int,
    PixelType::Dimension);

  typedef FixedArray< double, PixelDimension > EigenValueArrayType;
  typedef SymmetricSecondRankTensor< double, ImageDimension > HessianMatrixType;

  typedef Matrix< double, ImageDimension, ImageDimension > EigenMatrixType;
  typedef Image< EigenMatrixType, ImageDimension>  EigenMatrixImageType;
  typedef typename EigenMatrixImageType::Pointer EigenMatrixPointer;

  typedef SymmetricEigenAnalysis< HessianMatrixType, EigenValueArrayType, EigenMatrixType > EigenCalculatorType;
  typedef ImageRegionConstIteratorWithIndex< ImageType > ConstIteratorType;
  typedef ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(DoubleConvertibleToOutputCheck,
                  (Concept::Convertible<double, PixelType>));
  /** End concept checking */
#endif

  /** Set/Get macros for Alpha */
  itkSetMacro(SigmaMin, double);
  itkGetMacro(SigmaMin, double);

  /** Set/Get macros for Beta */
  itkSetMacro(SigmaMax, double);
  itkGetMacro(SigmaMax, double);

protected:
  MembraneScaleSelectionImageFilter();
  ~MembraneScaleSelectionImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ComputeAbsoluteMagnitudeOrdering( EigenValueArrayType& eigenValue, EigenMatrixType& eigenMatrix );
  double GetStructureMeasure( double& Lambda1, double& Lambda2, double& Lambda3 );

  virtual void BeforeThreadedGenerateData();
  virtual void AfterThreadedGenerateData();
  virtual void ThreadedGenerateData(const RegionType& windowRegion,
                                    ThreadIdType threadId);
  void GenerateInputRequestedRegion();
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));


private:
  MembraneScaleSelectionImageFilter(const Self&);
  void operator=(const Self&); //purposely not implemented

  double m_SigmaMin;
  double m_SigmaMax;
  OutputImagePointer  m_SigmaImage;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMembraneScaleSelectionImageFilter.txx"
#endif

#endif
