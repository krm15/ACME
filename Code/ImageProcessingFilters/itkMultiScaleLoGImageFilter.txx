/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleLoGImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.13 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiScaleLoGImageFilter_txx
#define __itkMultiScaleLoGImageFilter_txx

#include "itkMultiScaleLoGImageFilter.h"
#include "vnl/vnl_math.h"

#define EPSILON  1e-03

namespace itk
{

/**
 * Constructor
 */
template <typename TInputImage, typename TOutputImage >
MultiScaleLoGImageFilter
<TInputImage,TOutputImage>
::MultiScaleLoGImageFilter()
{
  m_SigmaMin = 0.2;
  m_SigmaMax = 2.0;
  m_NumberOfSigmaSteps = 4;

  m_LaplacianFilter = LaplacianFilterType::New();
}

template <typename TInputImage, typename TOutputImage >
void
MultiScaleLoGImageFilter
<TInputImage,TOutputImage>
::GenerateData()
{
  // Allocate the output
  this->GetOutput()->SetBufferedRegion(
                 this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();
  this->GetOutput()->FillBuffer( -1000 );

  InputImageConstPointer input = this->GetInput();

  m_LaplacianFilter->SetInput( input );
  m_LaplacianFilter->SetNormalizeAcrossScale( false );

  double sigma = m_SigmaMin;
  int scaleLevel = 1;
  while ( sigma <= m_SigmaMax )
    {
    std::cout << "Computing LoG for sigma = " << sigma << std::endl;
    m_LaplacianFilter->SetSigma( sigma );
    m_LaplacianFilter->Update();
    this->UpdateMaximumResponse();
    sigma  = this->ComputeSigmaValue( scaleLevel );
    scaleLevel++;
    }
}

template <typename TInputImage, typename TOutputImage >
void
MultiScaleLoGImageFilter
<TInputImage,TOutputImage>
::UpdateMaximumResponse()
{

  IteratorType oit(this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion());
  oit.GoToBegin();

  IteratorType it(m_LaplacianFilter->GetOutput(),
    m_LaplacianFilter->GetOutput()->GetLargestPossibleRegion());
  it.GoToBegin();

  while( !oit.IsAtEnd() )
    {
    if( oit.Value() < it.Value() )
      {
      oit.Value() = it.Value();
      }
    ++oit;
    ++it;
    }
}

template <typename TInputImage, typename TOutputImage >
double
MultiScaleLoGImageFilter
<TInputImage,TOutputImage>
::ComputeSigmaValue( int ScaleLevel )
{
  double stepSize =
     ( vcl_log( m_SigmaMax )  - vcl_log( m_SigmaMin) ) / m_NumberOfSigmaSteps;

  if( stepSize <= 1e-10 )
    {
    stepSize = 1e-10;
    }

  return ( vcl_exp( vcl_log (m_SigmaMin) + stepSize * ScaleLevel) );
}

template <typename TInputImage, typename TOutputImage >
void
MultiScaleLoGImageFilter
<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "SigmaMin:  " << m_SigmaMin << std::endl;
  os << indent << "SigmaMax:  " << m_SigmaMax  << std::endl;
}


} // end namespace itk

#endif
