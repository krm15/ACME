/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 720 $  // Revision of last commit
  Date: $Date: 2009-09-30 20:23:54 -0400 (Wed, 30 Sep 2009) $  // Date of last commit
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

#ifndef __itkGradientWeightedDistanceImageFilter_txx
#define __itkGradientWeightedDistanceImageFilter_txx

#include "itkGradientWeightedDistanceImageFilter.h"

namespace itk
{
//  Software Guide : BeginCodeSnippet
template< class TFeatureImage,class TInputImage,class TSegmentImage >
GradientWeightedDistanceImageFilter< TFeatureImage, TInputImage, TSegmentImage >
::GradientWeightedDistanceImageFilter()
{
  m_UseLevelSet = true;
  m_LargestCellRadius = 5.0;
  m_NucleiSigma = 0.4;
  m_Alpha = 1;
  m_Beta = 0.0;

  m_ForegroundMap = NULL;
  m_DistanceMap = NULL;
  m_Gradient = NULL;

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TInputImage::New() );
}

template< class TFeatureImage,class TInputImage,class TSegmentImage >
void
GradientWeightedDistanceImageFilter< TFeatureImage, TInputImage, TSegmentImage >
::DistanceMap()
{
  MaurerPointer m_Maurer = MaurerType::New();
  m_Maurer->SetInput ( m_ForegroundMap );
  m_Maurer->SetSquaredDistance ( 0 );
  m_Maurer->SetUseImageSpacing ( 1 );
  m_Maurer->SetInsideIsPositive ( 0 );
  m_Maurer->Update();

  ImagePointer image;
  if ( m_UseLevelSet )
    {
    // Nuclei only
    AbsFilterPointer  m_absFilter = AbsFilterType::New();
    m_absFilter->SetInput ( m_Maurer->GetOutput() );
    m_absFilter->Update();
    image = m_absFilter->GetOutput();
    image->DisconnectPipeline();
    }
  else
    {
      // Membranes only
      image = m_Maurer->GetOutput();
      image->DisconnectPipeline();
    }

  IteratorType It( image, image->GetLargestPossibleRegion() );
  It.GoToBegin();
  while( !It.IsAtEnd() )
    {
    // Set pixels on or inside the membrane to 0
    if ( ( !m_UseLevelSet ) && ( It.Get() < 0 ) )
      {
      It.Set( 0 );
      }

    if ( vcl_abs( It.Get() ) > m_LargestCellRadius )
      {
      It.Set( m_LargestCellRadius );
      }
    ++It;
    }

  m_DistanceMap = image;
}

template< class TFeatureImage,class TInputImage,class TSegmentImage >
void
GradientWeightedDistanceImageFilter< TFeatureImage, TInputImage, TSegmentImage >::
GenerateData()
{
  DistanceMap();

  ImagePointer outputImg;

  if ( m_Beta > 0 )
  {
    MultiScaleLoGFilterPointer log = MultiScaleLoGFilterType::New();
    log->SetInput( this->GetInput () );
    log->SetSigmaMin( 0.4  );
    log->SetSigmaMax( 2.0 );
    log->SetNumberOfSigmaSteps( 5 );
    log->Update();

    RescaleFilterPointer rescale = RescaleFilterType::New();
    rescale->SetInput( log->GetOutput() );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 1 );
    rescale->Update();

    outputImg = rescale->GetOutput();
    outputImg->DisconnectPipeline();
  }
  else
  {
    outputImg = ImageType::New();
    outputImg->CopyInformation( this->GetInput() );
    outputImg->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
    outputImg->Allocate();
    outputImg->FillBuffer( 1.0 );
  }

  double gmax, gmin;

  // Compute the gradient magnitude
  if ( m_Alpha > 0 )
  {
    GradientFilterPointer m_gradientMagnitude = GradientFilterType::New();
    m_gradientMagnitude->SetInput ( this->GetInput () );
    m_gradientMagnitude->SetSigma ( m_NucleiSigma );
    m_gradientMagnitude->Update();
    m_Gradient = m_gradientMagnitude->GetOutput();
    m_Gradient->DisconnectPipeline();

    MinMaxCalculatorPointer minMax3 = MinMaxCalculatorType::New();
    minMax3->SetImage( m_Gradient );
    minMax3->ComputeMaximum();
    gmax = static_cast< double >( minMax3->GetMaximum() );

    MinMaxCalculatorPointer minMax4 = MinMaxCalculatorType::New();
    minMax4->SetImage( m_Gradient );
    minMax4->ComputeMinimum();
    gmin = static_cast< double >( minMax4->GetMinimum() );
  }
  else
  {
    m_Gradient = outputImg;
    gmax = 1;
    gmin = 0;
  }

  IteratorType It ( outputImg,
                    outputImg->GetLargestPossibleRegion() );
  IteratorType dIt ( m_DistanceMap,
                    m_DistanceMap->GetLargestPossibleRegion() );
  IteratorType gIt ( m_Gradient,
                     m_Gradient->GetLargestPossibleRegion() );

  double p, q;
  gIt.GoToBegin();
  dIt.GoToBegin();
  It.GoToBegin();
  while ( !It.IsAtEnd() )
    {
    p = vcl_exp( m_Alpha * (gmax - static_cast<double>(gIt.Get()) )/ ( gmax - gmin ) );
    q = vcl_exp( - m_Beta * static_cast<double>( It.Get() ) );

    It.Set ( dIt.Get() * p * q );

    ++gIt;
    ++dIt;
    ++It;
    }

    this->GraftOutput ( outputImg );
}

template< class TFeatureImage,class TInputImage,class TSegmentImage >
void
GradientWeightedDistanceImageFilter< TFeatureImage, TInputImage, TSegmentImage >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name:   " << this->GetNameOfClass( ) << std::endl;
  os << indent << "Use LevelSet: " << m_UseLevelSet << std::endl;
}

} /* end namespace itk */

#endif
