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

  BallGeneratorPointer ball = BallGeneratorType::New();
  ball->SetSigma ( this->m_Sigma );
  ball->SetInput( stick->GetOutputImage() );
  ball->SetSpacing( spacing );
  ball->ComputeBallField();

  this->m_VotingField.push_back( ball->GetOutputImage() );
  this->m_VotingField.push_back( stick->GetOutputImage() );

  double radius = vcl_floor( vcl_sqrt( -vcl_log(0.01) * (this->m_Sigma) * (this->m_Sigma) ) );

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
::IntegrateVotes()
{
  // Fill orientation lookup image with lists of similarly oriented tokens
  // Change ComposeVoteFilter to take a m_VotingField as input and m_Lookup image
  ComposeVotesFilterPointer compose = ComposeVotesFilterType::New();
  compose->SetInput( this->m_Lookup );
  compose->SetVotingField( this->m_VotingField );
  compose->SetSaliencyImage( this->m_SaliencyImage );
  compose->SetOutputImage( this->m_Output );
  compose->SetNumberOfThreads( this->GetNumberOfThreads() );
  compose->Update();
}



template< class TInputImage >
void
TensorVoting2D< TInputImage >
::ComputeLookup()
{
  RegionType region = this->GetInput()->GetLargestPossibleRegion();

  IndexType index, index2;
  PixelType p;
  VectorType u, saliency;
  MatrixType eigenMatrix, R;
  PointType pt_sph;
  bool token;

  if ( this->m_LowMemoryFilter )
  {
    saliency.Fill( 0.0 );

    // Allocate saliency image
    this->m_SaliencyImage = VectorImageType::New();
    this->m_SaliencyImage->SetRegions( region );
    this->m_SaliencyImage->CopyInformation( this->GetInput() );
    this->m_SaliencyImage->Allocate();
    this->m_SaliencyImage->FillBuffer( saliency );
  }
  else
  {
    if  ( ( !this->m_EigenMatrixImage ) || ( !this->m_SaliencyImage ) )
    {
      typename SaliencyFilterType::Pointer saliencyFilter = SaliencyFilterType::New();
      saliencyFilter->SetInput( this->GetInput() );
      saliencyFilter->SetComputeEigenMatrix( 1 );
      saliencyFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
      saliencyFilter->Update();
      this->m_SaliencyImage = saliencyFilter->GetOutput();
      this->m_SaliencyImage->DisconnectPipeline();
      this->m_EigenMatrixImage = saliencyFilter->GetEigenMatrix();
      this->m_EigenMatrixImage->DisconnectPipeline();
    }
  }

  SizeType nSize = this->m_Lookup->GetLargestPossibleRegion().GetSize();


  // Compute an image of lists with similar pixel types
  ConstIteratorType   It( this->GetInput(), region );
  TokenIteratorType  tIt( this->m_TokenImage, region );
  VectorIteratorType sIt( this->m_SaliencyImage, region );
  It.GoToBegin();
  tIt.GoToBegin();
  sIt.GoToBegin();

  RandomGeneratorPointer rand = RandomGeneratorType::New();
  rand->Initialize();

  while( !It.IsAtEnd() )
  {
    index = It.GetIndex();
    p = It.Get();
    token = tIt.Get();

    if ( !this->m_LowMemoryFilter )
    {
      eigenMatrix = this->m_EigenMatrixImage->GetPixel( index );
      saliency = sIt.Get();
    }
    else
    {
      // Compute saliency and eigenMatrix
      this->m_EigenCalculator.ComputeEigenValuesAndVectors( p,
        saliency, eigenMatrix );

      for( unsigned int j = ImageDimension-1; j > 0; j-- )
      {
        saliency[j] = vcl_abs( saliency[j] - saliency[j-1] );
      }
      sIt.Set( saliency );
    }

    if ( ( saliency[0] > 0.001 ) && ( token ) )
    {
      for(unsigned int i = 0; i < ImageDimension; i++)
      {
        index2[i] = rand->GetIntegerVariate( nSize[i]-1 );
      }
      IdType *ll = &( this->m_Lookup->GetPixel( index2 )[0] );
      ll->push_back( index );
    }

    for (unsigned int j = 1; j < ImageDimension; j++ )
    {
      if ( ( saliency[j] > 0.001 ) && ( token ) )
      {
        u = eigenMatrix[1];
        if ( u[0] < 0 )
        {
          u = -u;
        }

        pt_sph[0] = Superclass::ComputeTheta( u );
        pt_sph[1] = 1.0;

        this->m_Lookup->TransformPhysicalPointToIndex( pt_sph, index2 );
        IdType *ll = &( this->m_Lookup->GetPixel( index2 )[j] );
        ll->push_back( index );
      }
    }

    ++sIt;
    ++It;
    ++tIt;
  }
}

} // end namespace itk

#endif
