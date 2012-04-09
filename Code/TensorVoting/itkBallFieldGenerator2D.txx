/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStickFieldGenerator2D.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 23:24:23 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBallFieldGenerator2D_txx
#define __itkBallFieldGenerator2D_txx

#include "itkBallFieldGenerator2D.h"

namespace itk
{
/**
 * Constructor
 */
template<class TOutputImage>
BallFieldGenerator2D<TOutputImage>
::BallFieldGenerator2D()
{}

template<class TOutputImage>
void
BallFieldGenerator2D<TOutputImage>
::Initialize(void)
{
  this->AllocateOutput();

  // Compute stick at a highest resolution and use that
  // to interpolate the ball field
  double sp = NumericTraits<double>::max();
  SpacingType spacing;
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    if (sp > static_cast< double >( this->m_Spacing[i] ) )
    {
       sp = this->m_Spacing[i];
    }
  }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    spacing[i] = sp;

  if ( !this->m_VotingField )
  {
    StickGeneratorPointer stick = StickGeneratorType::New();
    stick->SetSigma ( this->m_Sigma );
    stick->SetSpacing( spacing );
    stick->ComputeStickField(); 
    this->m_VotingField = stick->GetOutputImage();
    this->m_VotingField->DisconnectPipeline();
  }
}

template<class TOutputImage>
void
BallFieldGenerator2D<TOutputImage>
::ComputeBallField(void)
{
  this->Initialize();

  RegionType region = this->m_Output->GetLargestPossibleRegion();

  MatrixType R;
	double m_Theta;
  double count = 0, weight = 1;
	// Add the contributions of each oriented voting field
	for( unsigned int i = 0; i < 360; i++, count += 1 )
	{
		m_Theta = i*vnl_math::pi/180;

    R[0][0] = vcl_cos( m_Theta );
    R[0][1] = -vcl_sin( m_Theta );
    R[1][0] = vcl_sin( m_Theta );
    R[1][1] = vcl_cos( m_Theta );

    this->Update( this->m_Output, R, this->m_Spacing, region, weight );
  }

  this->Postprocess( this->m_Output, count );
}


template<class TOutputImage>
void
BallFieldGenerator2D<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
