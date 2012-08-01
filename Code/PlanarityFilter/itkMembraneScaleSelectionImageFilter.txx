/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMembraneScaleSelectionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 20:59:44 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMembraneScaleSelectionImageFilter_txx
#define __itkMembraneScaleSelectionImageFilter_txx

#include "itkMembraneScaleSelectionImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

#define EPSILON  1e-03

namespace itk
{
template < typename  TImage, typename TOutputImage >
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::MembraneScaleSelectionImageFilter()
{
  m_SigmaMin = 0.4;
  m_SigmaMax = 1.4;
}


template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >::
ComputeAbsoluteMagnitudeOrdering( EigenValueArrayType& eigenValue, EigenMatrixType& eigenMatrix )
{
  double Lambda;

  for ( unsigned int i=1; i <=2; i++ )
  {
    for ( unsigned int j = 0; j <= 1; j++ )
    {
      if ( vnl_math_abs( eigenValue[j] ) > vnl_math_abs( eigenValue[j+1] ) )
      {
        Lambda = eigenValue[j];
        eigenValue[j] = eigenValue[j+1];
        eigenValue[j+1] = Lambda;

        for( unsigned int k = 0; k < 3; k++ )
        {
          Lambda = eigenMatrix[j][k];
          eigenMatrix[j][k] = eigenMatrix[j+1][k];
          eigenMatrix[j+1][k] = Lambda;
        }
      }
    }
  }
}


template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    ImagePointer image = const_cast< ImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::BeforeThreadedGenerateData()
{}

template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::AfterThreadedGenerateData()
{}


template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::ThreadedGenerateData(const RegionType& windowRegion, ThreadIdType threadId)
{
  ImageConstPointer input = this->GetInput();
  OutputImagePointer output = this->GetOutput();

  EigenCalculatorType m_EigenCalculator;
  m_EigenCalculator.SetDimension( ImageDimension );

  ConstIteratorType iIt( input, windowRegion );
  OutputIteratorType oIt( output, windowRegion );

  iIt.GoToBegin();
  oIt.GoToBegin();

  // walk the region of eigen values and get the vesselness measure
  EigenValueArrayType eigenValue;
  EigenMatrixType     eigenMatrix;
  double sigma;

  while ( !iIt.IsAtEnd() )
    {
    // Get the eigen value in ascending order
    // Preserves sign of eigenValue but eigenVector ordering is not useful
    m_EigenCalculator.SetOrderEigenMagnitudes(0);
    m_EigenCalculator.ComputeEigenValuesAndVectors( iIt.Get(), eigenValue, eigenMatrix );

    // Recreate magnitude ordering with sign
    ComputeAbsoluteMagnitudeOrdering( eigenValue, eigenMatrix );

    sigma = m_SigmaMin + (m_SigmaMax - m_SigmaMin) * vcl_abs( eigenMatrix[2][2] );

    oIt.Set( sigma );
    ++iIt;
    ++oIt;
    }
}


template < typename  TImage, typename TOutputImage >
void
MembraneScaleSelectionImageFilter< TImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
