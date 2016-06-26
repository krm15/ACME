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

#ifndef __itkTensorVoting_txx
#define __itkTensorVoting_txx

#include "itkTensorVoting.h"

namespace itk
{

/**
 * Constructor
 */
template< class TInputImage >
TensorVoting<TInputImage >
::TensorVoting()
{
  m_Sigma               = 5.0;
  m_UseSparseVoting     = 0;
  m_TokenImage          = ITK_NULLPTR;
  m_Lookup              = ITK_NULLPTR;
  m_SaliencyImage       = ITK_NULLPTR;
  m_EigenMatrixImage    = ITK_NULLPTR;
  m_Output              = ITK_NULLPTR;

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TInputImage::New() );
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::InitializeLookupImages()
{
  IdType emptyList;
  emptyList.clear();
  VectorInternalType vec;
  vec.Fill( emptyList );

  SpacingType nSpacing;
  PointType nOrigin;
  SizeType nSize;
  IndexType nIndex;

  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    nSpacing[i] = 1.0;
    nOrigin[i] = -90.0;
    nSize[i] = 180;
    nIndex[i] = 0;
  }

  nOrigin[ ImageDimension-1 ] = 1.0;
  nSize[ ImageDimension-1 ] = 1;

  RegionType nRegion;
  nRegion.SetSize( nSize );
  nRegion.SetIndex( nIndex );

  m_Lookup = InternalImageType::New();
  m_Lookup->SetRegions( nRegion );
  m_Lookup->SetOrigin( nOrigin );
  m_Lookup->SetSpacing( nSpacing );
  m_Lookup->Allocate();
  m_Lookup->FillBuffer( vec );
}


template< class TInputImage >
double
TensorVoting< TInputImage >
::ComputeTheta( VectorType& u )
{
  double theta;
  if ( u[0] == 0 )
  {
    theta = vnl_math::pi_over_2;
  }
  else
  {
    theta = vcl_atan(  ( u[1] /vcl_abs( u[0] ) ) );
  }

  // Determine the quadrant of theta
  if ( u[0] < 0 )
  {
    theta = vnl_math::pi - theta;
  }

  theta = vnl_math::pi_over_2 - theta;

  return theta;
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::GenerateData()
{
  InitializeLookupImages();

  ImageConstPointer input = this->GetInput();
  RegionType region = input->GetLargestPossibleRegion();
  
  // A zero tensor
  MatrixType ZeroTensor;
  ZeroTensor.Fill( 0 );

  // Allocate the output image of voted tensors
  m_Output = InputImageType::New();
  m_Output->SetRegions( region );
  m_Output->CopyInformation( input );
  m_Output->Allocate();
  m_Output->FillBuffer( ZeroTensor );
  
  // Initialize the voting fields
  this->InitializeVotingFields();
  std::cout << "Voting fields initialized" << std::endl;

  // Token image of type bool specifying sparse or dense tokens
  if ( !m_TokenImage )
  {    
    m_TokenImage = TokenImageType::New();
    m_TokenImage->SetRegions( region );
    m_TokenImage->SetOrigin( input->GetOrigin() );
    m_TokenImage->SetSpacing( input->GetSpacing() );
    m_TokenImage->Allocate();
    if ( m_UseSparseVoting )
      m_TokenImage->FillBuffer( false );
    else
      m_TokenImage->FillBuffer( true );
  }

  // Fill orientation lookup image with lists of similarly oriented voxels
  ComputeLookup();
  std::cout << "Computing lookup finished" << std::endl;

  // Fill orientation lookup image with lists of similarly oriented tokens
  // Change ComposeVoteFilter to take a m_VotingField as input and m_Lookup image
  ComposeVotesFilterPointer compose = ComposeVotesFilterType::New();
  compose->SetInput( m_Lookup );
  compose->SetVotingField( m_VotingField );
  compose->SetSaliencyImage( m_SaliencyImage );
  compose->SetOutputImage( m_Output );
  compose->SetNumberOfThreads( this->GetNumberOfThreads() );
  compose->Update();
  std::cout << "Integration completed" << std::endl;

  // Graft the output image
  this->GraftOutput( m_Output );
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
