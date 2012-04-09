/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiscaleStructMeasureImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/04/01 23:13:46 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiscaleStructMeasureImageFilter_h
#define __itkMultiscaleStructMeasureImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkStructMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
/**\class MultiscaleStructMeasureImageFilter
 * \brief A filter to enhance 3D vascular structures using Hessian
 *         eigensystem in a multiscale framework
 *
 * The plate measure is based on the analysis of the the Hessian
 * eigen system. The filter takes an
 * image of any pixel type and generates a Hessian image pixels at different
 * scale levels. The plate measure is computed from the Hessian image
 * at each scale level and the best response is selected.  The plate
 * measure is computed using StructMeasureImageFilter.
 *
 * Minimum and maximum sigma value can be set using SetMinSigma and SetMaxSigma
 * methods respectively. The number of scale levels is set using
 * SetNumberOfSigmaSteps method. Exponentially distributed scale levels are
 * computed within the bound set by the minimum and maximum sigma values
 *
 *
 * \par References
 *  Mosaliganti, K, F. Janoos, Gelas, A, Machiraju, R and Megason, S (2009). Membrane Enhancing
 *  Diffusion: A Scale Space Representation of Membrane Structures. ISBI 2010
 *
 * \sa MultiscaleStructMeasureImageFilter
 * \sa Hessian3DToMembranenessMeasureImageFilter
 * \sa HessianSmoothedRecursiveGaussianImageFilter
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 *
 * \ingroup IntensityImageFilters TensorObjects
 *
 */
template <class TInputImage,
          class TOutputImage = TInputImage >
class ITK_EXPORT MultiscaleStructMeasureImageFilter
: public
ImageToImageFilter< TInputImage,TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiscaleStructMeasureImageFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>              Superclass;

  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename OutputImageType::RegionType OutputRegionType;

  /** Image dimension = 3. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                   ::itk::GetImageDimension<InputImageType>::ImageDimension);

  /** Update image buffer that holds the best vesselness response */
  typedef HessianRecursiveGaussianImageFilter< InputImageType > HessianFilterType;
  typedef typename HessianFilterType::Pointer HessianFilterPointer;

  typedef StructMeasureImageFilter< double > StructFilterType;
  typedef typename StructFilterType::Pointer StructFilterPointer;
  typedef typename StructFilterType::OutputImageType StructOutputImageType;

  // Define image of matrix pixel type
  typedef Matrix< double, ImageDimension, ImageDimension> EigenMatrixType;
  typedef Image< EigenMatrixType, ImageDimension>  EigenMatrixImageType;
  typedef typename EigenMatrixImageType::Pointer EigenMatrixPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get macros for Alpha */
  itkSetMacro(SigmaMin, double);
  itkGetMacro(SigmaMin, double);

  /** Set/Get macros for Beta */
  itkSetMacro(SigmaMax, double);
  itkGetMacro(SigmaMax, double);

  /** Set/Get macros for Number of Scales */
  itkSetMacro(NumberOfSigmaSteps, int);
  itkGetMacro(NumberOfSigmaSteps, int);

  /** Set/Get macros for object type */
  itkSetMacro(ObjectType, unsigned char);
  itkGetMacro(ObjectType, unsigned char);

  StructFilterPointer GetStructFilter()
    {
    return m_StructFilter;
    }

  EigenMatrixPointer GetEigenMatrix()
    {
    return m_EigenMatrixImage;
    }

protected:
  MultiscaleStructMeasureImageFilter();
  ~MultiscaleStructMeasureImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Generate Data */
  void GenerateData( void );

private:
  void UpdateMaximumResponse();

  double ComputeSigmaValue( int scaleLevel );

  void   AllocateUpdateBuffer();

  //purposely not implemented
  MultiscaleStructMeasureImageFilter(const Self&);
  void operator=(const Self&); //purposely not implemented

  double                             m_SigmaMin;
  double                             m_SigmaMax;
  int                                m_NumberOfSigmaSteps;
  unsigned char                      m_ObjectType;

  StructFilterPointer                m_StructFilter;
  HessianFilterPointer               m_HessianFilter;
  EigenMatrixPointer                 m_EigenMatrixImage;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiscaleStructMeasureImageFilter.txx"
#endif

#endif
