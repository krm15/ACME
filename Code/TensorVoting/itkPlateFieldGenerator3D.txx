/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStickFieldGenerator3D.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 23:24:23 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPlateFieldGenerator3D_txx
#define __itkPlateFieldGenerator3D_txx

#include "itkPlateFieldGenerator3D.h"

namespace itk
{

/**
 * Constructor
 */
template<class TOutputImage>
PlateFieldGenerator3D<TOutputImage>
::PlateFieldGenerator3D()
{}

template<class TOutputImage>
void
PlateFieldGenerator3D<TOutputImage>
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
PlateFieldGenerator3D<TOutputImage>
::ComputePlateField(void)
{
  this->Initialize();
  RegionType region = this->m_VotingField->GetLargestPossibleRegion();
  SpacingType sp = this->m_VotingField->GetSpacing();

  ImagePointer temp = ImageType::New();
  temp->SetRegions( region );
  temp->SetOrigin( this->m_Origin );
  temp->SetSpacing( sp );
  temp->Allocate();

  MatrixType R;
	double m_Theta;
	unsigned int k = 5;
  double count = 0, weight = 1;
	// Add the contributions of each oriented voting field
	for( unsigned int i = 0; i < 360/k; i++, count += 1 )
	{
		m_Theta = i*k*vnl_math::pi/180;
    rMatrixHelper.ComputeRotationMatrixWithAxis( 1, 0, 0, m_Theta, R);
		this->Update( temp, R, sp, region, weight );
	}
  this->Postprocess( temp, count );

  // Rotate along y by 90 to orient the normal along z
  // The stick rotated in x-y plane
	m_Theta = -vnl_math::pi_over_2;
  rMatrixHelper.ComputeRotationMatrixWithAxis( 0, 1, 0, m_Theta, R);

  OrientedTensorGeneratorPointer orientedTensor = 
    OrientedTensorGeneratorType::New();
	orientedTensor->SetInput( temp );
	orientedTensor->SetRotationMatrix( R );
  orientedTensor->SetOutputSpacing( this->m_Spacing );
  orientedTensor->SetOutputRegion(   
    this->m_Output->GetLargestPossibleRegion() );
	orientedTensor->Update();

	this->m_Output = orientedTensor->GetOutput();
}


template<class TOutputImage>
void
PlateFieldGenerator3D<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
