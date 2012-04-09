/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiscaleStructMeasureImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.13 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiscaleStructMeasureImageFilter_txx
#define __itkMultiscaleStructMeasureImageFilter_txx

#include "itkMultiscaleStructMeasureImageFilter.h"
#include "vnl/vnl_math.h"

#define EPSILON  1e-03

namespace itk
{

/**
 * Constructor
 */
template <typename TInputImage, typename TOutputImage >
MultiscaleStructMeasureImageFilter
< TInputImage,TOutputImage >
::MultiscaleStructMeasureImageFilter()
{
  m_SigmaMin = 0.2;
  m_SigmaMax = 2.0;
  m_NumberOfSigmaSteps = 10;
  m_ObjectType = 0;

  m_HessianFilter    = HessianFilterType::New();
  m_StructFilter      = StructFilterType::New();
  m_EigenMatrixImage = EigenMatrixImageType::New();

  //Turn off membraneness measure scaling
  m_StructFilter->SetScaleMembranenessMeasure( false );
}

template <typename TInputImage, typename TOutputImage >
void
MultiscaleStructMeasureImageFilter
<TInputImage,TOutputImage>
::AllocateUpdateBuffer()
{
  EigenMatrixType p;
  p.Fill( 0 );

  m_EigenMatrixImage->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
  m_EigenMatrixImage->CopyInformation( this->GetInput() );
  m_EigenMatrixImage->Allocate();
  m_EigenMatrixImage->FillBuffer( p );
}

template <typename TInputImage, typename TOutputImage >
void
MultiscaleStructMeasureImageFilter
<TInputImage,TOutputImage>
::GenerateData()
{
  // Allocate the output
  this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();

  // Allocate the buffer
  AllocateUpdateBuffer();

  this->m_HessianFilter->SetInput( this->GetInput() );
  this->m_HessianFilter->SetNormalizeAcrossScale( true );

  double sigma = m_SigmaMin;
  int scaleLevel = 1;
  while ( sigma <= m_SigmaMax )
    {
    std::cout << "Computing structure with sigma= " << sigma << std::endl;
    m_HessianFilter->SetSigma( sigma );
    m_StructFilter->SetInput ( m_HessianFilter->GetOutput() );
    m_StructFilter->SetObjectType( m_ObjectType );
    m_StructFilter->Update();

    this->UpdateMaximumResponse();
    sigma  = this->ComputeSigmaValue( scaleLevel );
    scaleLevel++;
    }
}

template <typename TInputImage, typename TOutputImage >
void
MultiscaleStructMeasureImageFilter
<TInputImage,TOutputImage>
::UpdateMaximumResponse()
{
  OutputRegionType region = this->GetOutput()->GetLargestPossibleRegion();
  ImageRegionIterator< OutputImageType > oit( this->GetOutput(), region );

  ImageRegionIterator<StructOutputImageType>
    it(m_StructFilter->GetOutput(), region );

  ImageRegionIterator<EigenMatrixImageType>
    mit( m_EigenMatrixImage, region );

  ImageRegionIterator<EigenMatrixImageType>
    hit( m_StructFilter->GetEigenMatrix(), region );

  oit.GoToBegin();
  it.GoToBegin();
  mit.GoToBegin();
  hit.GoToBegin();

  while( !oit.IsAtEnd() )
    {
    if( oit.Value() < static_cast< OutputPixelType >( it.Value() ) )
      {
      oit.Value() = static_cast< OutputPixelType >( it.Value() );
      mit.Set( hit.Get() );
      }
    ++oit;
    ++it;
    ++hit;
    ++mit;
    }
}

template <typename TInputImage, typename TOutputImage >
double
MultiscaleStructMeasureImageFilter
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
MultiscaleStructMeasureImageFilter
<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "SigmaMin:  " << m_SigmaMin << std::endl;
  os << indent << "SigmaMax:  " << m_SigmaMax  << std::endl;
}


} // end namespace itk

#endif
