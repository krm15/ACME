/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleLoGImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/04/01 23:13:46 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiScaleLoGImageFilter_h
#define __itkMultiScaleLoGImageFilter_h


#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
/**\class MultiScaleLoGImageFilter
 * \brief A filter to enhance 3D blob structures using LoG filter
 */
template <class TInputImage, class TOutputImage = TInputImage >
class ITK_EXPORT MultiScaleLoGImageFilter
: public ImageToImageFilter< TInputImage,TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiScaleLoGImageFilter                       Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>   Superclass;

  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  typedef TInputImage                                    InputImageType;
  typedef typename InputImageType::Pointer               InputImagePointer;
  typedef typename InputImageType::ConstPointer          InputImageConstPointer;
  typedef typename InputImageType::PixelType             InputPixelType;

  typedef TOutputImage                                   OutputImageType;
  typedef typename OutputImageType::Pointer              OutputImagePointer;
  typedef typename OutputImageType::PixelType            OutputPixelType;

  /** Image dimension = 3. */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  /** Update image buffer that holds the best laplacian response */
  typedef LaplacianRecursiveGaussianImageFilter< InputImageType, OutputImageType >
                                                        LaplacianFilterType;
  typedef typename LaplacianFilterType::Pointer         LaplacianFilterPointer;

  typedef ImageRegionIterator<OutputImageType>         IteratorType;

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


protected:
  MultiScaleLoGImageFilter();
  ~MultiScaleLoGImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData( void );

private:
  void UpdateMaximumResponse();
  double ComputeSigmaValue( int scaleLevel );

  //purposely not implemented
  MultiScaleLoGImageFilter(const Self&);
  void operator=(const Self&); //purposely not implemented

  double                             m_SigmaMin;
  double                             m_SigmaMax;
  int                                m_NumberOfSigmaSteps;

  LaplacianFilterPointer             m_LaplacianFilter;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleLoGImageFilter.txx"
#endif

#endif
