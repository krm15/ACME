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
#ifndef __itkBallFieldGenerator3D_txx
#define __itkBallFieldGenerator3D_txx

#include "itkBallFieldGenerator3D.h"

namespace itk
{
/**
 * Constructor
 */
template<class TOutputImage>
BallFieldGenerator3D<TOutputImage>
::BallFieldGenerator3D()
{}

template<class TOutputImage>
void
BallFieldGenerator3D<TOutputImage>
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
BallFieldGenerator3D<TOutputImage>
::ComputeBallField(void)
{
  this->Initialize();

  RegionType region = this->m_Output->GetLargestPossibleRegion();

  PointType p;
  MatrixType R;
  VectorType u;
	double m_Theta, m_Phi, weight, strength = 0;
	// Add the contributions of each oriented voting field
	double k = 5, l = 5;
  for( double i = 0; i <= 180; i += k )
  {
  	for( double j = 0; j < 360; j += l )
    {
      m_Theta = i * vnl_math::pi/180 - vnl_math::pi_over_2;
      m_Phi   = j * vnl_math::pi/180;
      weight = vcl_cos( m_Theta );
      strength += weight;

      u[0] = vcl_cos( m_Theta ) * vcl_cos( m_Phi );
      u[1] = vcl_cos( m_Theta ) * vcl_sin( m_Phi );
      u[2] = vcl_sin( m_Theta );

      rMatrixHelper.ComputeRotationMatrix( u, R );

      this->Update( this->m_Output, R, this->m_Spacing, region, weight );
	  }
  }
  this->Postprocess( this->m_Output, strength );
}


template<class TOutputImage>
void
BallFieldGenerator3D<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
