/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 585 $  // Revision of last commit
  Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
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

#ifndef __itkTensorVoting2D_txx
#define __itkTensorVoting2D_txx

#include "itkTensorVoting2D.h"

namespace itk
{

/**
 * Constructor
 */
template< class TInputImage >
TensorVoting2D<TInputImage >
::TensorVoting2D()
{}

template< class TInputImage >
void
TensorVoting2D< TInputImage >
::InitializeVotingFields(void)
{
  Superclass::InitializeVotingFields();

  SpacingType m_Spacing = this->GetInput()->GetSpacing();
  double sp = NumericTraits<double>::max();
  SpacingType spacing;
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    if (sp > static_cast< double >( m_Spacing[i] ) )
    {
      sp = m_Spacing[i];
    }
  }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    spacing[i] = sp;

  // Compute the stick and ball voting fields
  StickGeneratorPointer stick = StickGeneratorType::New();
  stick->SetSigma ( this->m_Sigma );
  stick->SetSpacing( spacing );
  stick->ComputeStickField();
  this->m_VotingField.push_back( stick->GetOutputImage() );

  BallGeneratorPointer ball = BallGeneratorType::New();
  ball->SetSigma ( this->m_Sigma );
  ball->SetSpacing( spacing );
  ball->ComputeBallField();
  this->m_VotingField.push_back( ball->GetOutputImage() );

  double radius = vcl_floor( vcl_sqrt( -vcl_log(0.01) 
    * (this->m_Sigma) * (this->m_Sigma) ) );
  SizeType size;
  SizeValueType rad;
  IndexType index;
  for( unsigned int i = 0; i< ImageDimension; i++ )
    {
    rad = static_cast< SizeValueType >( radius/m_Spacing[i] );
    size[i] = 2*rad + 1;
    index[i] = 0;
    }

  this->m_Region.SetSize( size );
  this->m_Region.SetIndex( index );
}


template< class TInputImage >
void
TensorVoting2D< TInputImage >
::ComputeThetaAndRotationAxis( VectorType& u, MatrixType& R )
{
  double theta;
  if ( u[0] == 0 )
  {
    if ( u[1] > 0 )
      theta = vnl_math::pi_over_2;
    else
      theta = -vnl_math::pi_over_2;
  }
  else
  {
    theta = vcl_atan(  vcl_abs( u[1] / u[0] ) );
  }

  // Determine the quadrant of theta
  if ( u[0] > 0 )
  {
    if ( u[1] < 0 )
    {
      theta = - theta;
    }
  }
  else
  {
    if ( u[1] < 0 )
    {
      theta += vnl_math::pi;
    }
    else
    {
      theta = vnl_math::pi - theta;
    }
  }

  theta = vnl_math::pi_over_2 - theta;

  R[0][0] = vcl_cos( theta );
  R[0][1] = -vcl_sin( theta );
  R[1][0] = vcl_sin( theta );
  R[1][1] = vcl_cos( theta );
}


template< class TInputImage >
void
TensorVoting2D< TInputImage >
::ComputeTensorVoting( IndexType& index, PixelType& p, int threadId )
{
  VectorType eigenVector, u, v;
  MatrixType eigenMatrix, R;
  PixelType pixel;
  PointType pt;

  this->m_Output->TransformIndexToPhysicalPoint( index, pt );

  // Compute eigen values and eigen matrix
  this->m_EigenCalculator.SetOrderEigenValues( 1 );
  this->m_EigenCalculator.ComputeEigenValuesAndVectors( p, eigenVector, eigenMatrix );

  if ( vcl_abs( eigenVector[1] - eigenVector[0] ) > 0.001 )
  {
    // Compute orientation (theta in 2D)
    u = eigenMatrix[ImageDimension-1];
    ComputeThetaAndRotationAxis( u, R );
    this->ComputeOrientedField( eigenVector[1] - eigenVector[0], pt, R, 0, threadId );
  }

  if ( vcl_abs( eigenVector[0] ) > 0.001 )
  {
    R.SetIdentity();
    this->ComputeOrientedField( eigenVector[0], pt, R, 1, threadId );
  }
}


template< class TInputImage >
void
TensorVoting2D< TInputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
