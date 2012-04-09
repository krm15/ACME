/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStructMeasureImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/06/12 22:59:15 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkStructMeasureImageFilter_h
#define __itkStructMeasureImageFilter_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
/** \class StructMeasureImageFilter
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
 * \sa MultiScaleStructMeasureImageFilter
 * \sa Hessian3DToMembranenessMeasureImageFilter
 * \sa HessianSmoothedRecursiveGaussianImageFilter
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 *
 * \ingroup IntensityImageFilters TensorObjects
 *
 */

template < typename  TPixel >
class ITK_EXPORT StructMeasureImageFilter : public
ImageToImageFilter< Image< SymmetricSecondRankTensor< double, 3 >, 3 >,
                                                  Image< TPixel, 3 > >
{
public:
  /** Standard class typedefs. */
  typedef StructMeasureImageFilter Self;

  typedef ImageToImageFilter<
          Image< SymmetricSecondRankTensor< double, 3 >, 3 >,
          Image< TPixel, 3 > >                 Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename InputImageType::RegionType  InputRegionType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::OutputImageRegionType 
    OutputImageRegionType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef TPixel                               OutputPixelType;


  /** Image dimension = 3. */
  itkStaticConstMacro( ImageDimension,
    unsigned int, ::itk::GetImageDimension<InputImageType>::ImageDimension);

  itkStaticConstMacro(InputPixelDimension, unsigned int,
    InputPixelType::Dimension);

  typedef  FixedArray< double, InputPixelDimension > EigenValueArrayType;
  typedef SymmetricSecondRankTensor< double, ImageDimension > HessianMatrixType;

  typedef Matrix< double, ImageDimension, ImageDimension > EigenMatrixType;
  typedef Image< EigenMatrixType, ImageDimension>  EigenMatrixImageType;
  typedef typename EigenMatrixImageType::Pointer EigenMatrixPointer;

  typedef SymmetricEigenAnalysis< HessianMatrixType, EigenValueArrayType, EigenMatrixType > EigenCalculatorType;
  typedef ImageRegionIteratorWithIndex< EigenMatrixImageType > EigenIteratorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get macros for Alpha */
  itkSetMacro(Alpha, double);
  itkGetMacro(Alpha, double);

  /** Set/Get macros for Beta */
  itkSetMacro(Beta, double);
  itkGetMacro(Beta, double);

  /** Set/Get macros for Gamma */
  itkSetMacro(Gamma, double);
  itkGetMacro(Gamma, double);

  /** Set/Get macros for C */
  itkSetMacro(C, double);
  itkGetMacro(C, double);

  /** Set/Get macros for Gamma */
  itkSetMacro(ObjectType, unsigned char);
  itkGetMacro(ObjectType, unsigned char);

  /** Macro to scale the vesselness measure with the
      largest eigenvalue or not */
  itkSetMacro( ScaleMembranenessMeasure, bool );
  itkGetMacro( ScaleMembranenessMeasure, bool );
	itkSetMacro( StructureIsPositive, bool );
  itkGetMacro( StructureIsPositive, bool );
  itkSetMacro( ComputeEigenVectors, bool );
  itkGetMacro( ComputeEigenVectors, bool );
  itkBooleanMacro( ScaleMembranenessMeasure );
  itkBooleanMacro( StructureIsPositive );
  itkBooleanMacro( ComputeEigenVectors );

  EigenMatrixPointer GetEigenMatrix()
    {
    return m_EigenMatrixImage;
    }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(DoubleConvertibleToOutputCheck,
                  (Concept::Convertible<double, OutputPixelType>));
  /** End concept checking */
#endif

protected:
  StructMeasureImageFilter();
  ~StructMeasureImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ComputeLambdaOrdering( double& Lambda1, double& Lambda2,
    double& Lambda3, EigenValueArrayType& eigenValue );
  double GetStructureMeasure( double& Lambda1, double& Lambda2, double& Lambda3 );

  virtual void BeforeThreadedGenerateData();
  virtual void AfterThreadedGenerateData();
  virtual void ThreadedGenerateData(const OutputImageRegionType& windowRegion,
                                    ThreadIdType threadId);

  /** HoughTransformRadialVotingImageFilter needs the entire input. Therefore
   * it must provide an implementation GenerateInputRequestedRegion().
   * \sa ProcessObject::GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** HoughTransformRadialVotingImageFilter's produces all the output.
   * Therefore, it must provide an implementation of
   * EnlargeOutputRequestedRegion.
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));


private:
  StructMeasureImageFilter(const Self&);
  void operator=(const Self&); //purposely not implemented

  EigenMatrixPointer  m_EigenMatrixImage;

  double m_Alpha;
  double m_Beta;
  double m_Gamma;
  double m_C;
  bool   m_ScaleMembranenessMeasure;
	bool   m_StructureIsPositive;
  bool   m_ComputeEigenVectors;
  unsigned char m_ObjectType;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStructMeasureImageFilter.txx"
#endif

#endif
