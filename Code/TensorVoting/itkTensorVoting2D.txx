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
{
  m_Sigma = 5.0;
  m_UseSparseVoting = 0;
  m_OrientedVotingField = ITK_NULLPTR;
  m_StickSaliencyImage = ITK_NULLPTR;
  m_BallSaliencyImage = ITK_NULLPTR;
  m_EigenMatrixImage = ITK_NULLPTR;

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TInputImage::New() );
}


template< class TInputImage >
void
TensorVoting2D< TInputImage >
::InitializeLookupImages()
{
  IdType emptyList;
  emptyList.clear();

  SpacingType nSpacing;
  nSpacing[0] = nSpacing[1] = 1.0;

  PointType nOrigin;
  nOrigin[0] = -90.0;
  nOrigin[2] = 1.0;

  SizeType nSize;
  nSize[0] = 180;
  nSize[1] = 1;

  IndexType nIndex;
  nIndex[0] = 0;
  nIndex[1] = 0;

  RegionType nRegion;
  nRegion.SetSize( nSize );
  nRegion.SetIndex( nIndex );

  m_LookupStick = InternalImageType::New();
  m_LookupStick->SetRegions( nRegion );
  m_LookupStick->SetOrigin( nOrigin );
  m_LookupStick->SetSpacing( nSpacing );
  m_LookupStick->Allocate();
  m_LookupStick->FillBuffer( emptyList );

  m_LookupBall = InternalImageType::New();
  m_LookupBall->SetRegions( nRegion );
  m_LookupBall->SetOrigin( nOrigin );
  m_LookupBall->SetSpacing( nSpacing );
  m_LookupBall->Allocate();
  m_LookupBall->FillBuffer( emptyList );
}

template< class TInputImage >
void
TensorVoting2D< TInputImage >
::InitializeVotingFields()
{
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
  {
    spacing[i] = sp;
  }

  // Compute the stick and ball voting fields
  StickGeneratorPointer stick = StickGeneratorType::New();
  stick->SetSigma ( this->m_Sigma );
  stick->SetSpacing( spacing );
  stick->ComputeStickField();
  this->m_VotingField.push_back( stick->GetOutputImage() );

  BallGeneratorPointer ball = BallGeneratorType::New();
  ball->SetSigma ( this->m_Sigma );
  ball->SetInput( stick->GetOutputImage() );
  ball->SetSpacing( spacing );
  ball->ComputeBallField();
  this->m_VotingField.push_back( ball->GetOutputImage() );

  double radius = vcl_floor( vcl_sqrt( -vcl_log(0.01) * (m_Sigma) * (m_Sigma) ) );

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
::ComputeLookup()
{
  RegionType region = this->GetInput()->GetLargestPossibleRegion();

  if  ( ( !m_EigenMatrixImage   ) ||
        ( !m_StickSaliencyImage ) ||
        ( !m_BallSaliencyImage  ) )
  {
    typename SaliencyFilterType::Pointer saliencyFilter = SaliencyFilterType::New();
    saliencyFilter->SetInput( this->GetInput() );
    saliencyFilter->SetComputeEigenMatrix( 1 );
    saliencyFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    saliencyFilter->Update();
    typename VectorImageType::Pointer saliencyImage = saliencyFilter->GetOutput();
    this->m_EigenMatrixImage = saliencyFilter->GetEigenMatrix();

    typename IndexFilterType::Pointer componentExtractor1 = IndexFilterType::New();
    componentExtractor1->SetInput( saliencyImage );
    componentExtractor1->SetIndex( 1 );
    componentExtractor1->Update();
    this->m_StickSaliencyImage = componentExtractor1->GetOutput();

    typename IndexFilterType::Pointer componentExtractor2 = IndexFilterType::New();
    componentExtractor2->SetInput( saliencyImage );
    componentExtractor2->SetIndex( 0 );
    componentExtractor2->Update();
    this->m_BallSaliencyImage = componentExtractor2->GetOutput();
  }
    
  IndexType index, index2;
  PixelType p;
  VectorType u, v;
  double stickSaliency, ballSaliency;
  MatrixType eigenMatrix, R;
  double theta;
  PointType pt;
  bool token;

  SizeType nSize = m_LookupBall->GetLargestPossibleRegion().GetSize();

  DoubleIteratorType sIt( m_StickSaliencyImage, region );
  DoubleIteratorType bIt( m_BallSaliencyImage, region );
  sIt.GoToBegin();
  bIt.GoToBegin();

  // Compute an image of lists with similar pixel types
  ConstIteratorType It( this->GetInput(), region );
  IteratorType eIt( m_EigenMatrixImage, region );
  TokenIteratorType tIt( m_TokenImage, region );
  It.GoToBegin();
  eIt.GoToBegin();
  tIt.GoToBegin();

  while( !It.IsAtEnd() )
  {
    index = It.GetIndex();
    p = It.Get();
    eigenMatrix = eIt.Get();
    token = tIt.Get();

    stickSaliency = sIt.Get();
    ballSaliency  = bIt.Get();

    if ( ( stickSaliency > 0.001 ) && ( token ) )
    {
      // Compute orientation (theta in 2D)
      // Eigen values were sorted based on magnitude.
      // So last vector is the normal to the surface
      u = eigenMatrix[ImageDimension-1]; // e2 is the normal
      if ( u[1] < 0 )
      {
        u = -u;
      }

      pt[0] = ComputeTheta( u );
      pt[1] = 1.0;

      m_LookupStick->TransformPhysicalPointToIndex( pt, index2 );
      IdType *ll = &( m_LookupStick->GetPixel( index2 ) );
      ll->push_back( index );
    }

    RandomGeneratorPointer rand = RandomGeneratorType::New();
    rand->Initialize();


    if ( ( ballSaliency > 0.001 ) && ( token ) )
    {
      for(unsigned int i = 0; i < ImageDimension; i++)
      {
        index2[i] = rand->GetIntegerVariate( nSize[i]-1 );
      }
      IdType *ll = &( m_LookupBall->GetPixel( index2 ) );
      ll->push_back( index );
    }

    ++sIt;
    ++bIt;

    ++It;
    ++eIt;
    ++tIt;
  }
}


template< class TInputImage >
double
TensorVoting2D< TInputImage >
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
TensorVoting2D< TInputImage >
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
  ComposeVotesFilterPointer composerStick = ComposeVotesFilterType::New();
  composerStick->SetInput( m_LookupStick );
  composerStick->SetVotingField( m_VotingField[0] );
  composerStick->SetSaliencyImage( m_StickSaliencyImage );
  composerStick->SetOutputImage( m_Output );
  composerStick->SetNumberOfThreads( this->GetNumberOfThreads() );
  composerStick->Update();
  std::cout << "Stick voting completed" << std::endl;

  // Add Ball voting field
  ComposeVotesFilterPointer composerBall = ComposeVotesFilterType::New();
  composerBall->SetInput( m_LookupBall );
  composerBall->SetVotingField( m_VotingField[1] );
  composerBall->SetVotingFieldTypeAsBall( true );
  composerBall->SetSaliencyImage( m_BallSaliencyImage );
  composerBall->SetOutputImage( m_Output );
  composerBall->SetNumberOfThreads( this->GetNumberOfThreads() );
  composerBall->Update();
  std::cout << "Ball voting completed" << std::endl;
  std::cout << "Integration completed" << std::endl;

  // Graft the output image
  this->GraftOutput( m_Output );
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
