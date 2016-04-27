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
  m_Sigma = 5.0;
  m_UseSparseVoting = 0;
  m_EigenCalculator.SetDimension( ImageDimension );

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TInputImage::New() );
}

template< class TInputImage >
void
TensorVoting< TInputImage >
::InitializeVotingFields(void)
{
  // A zero tensor
  MatrixType ZeroTensor;
  ZeroTensor.Fill( 0 );

  // Allocate the output image
  this->m_Output = this->GetOutput();
  this->m_Output->FillBuffer( ZeroTensor );
}

template< class TInputImage >
void
TensorVoting< TInputImage >
::ComputeOrientedField( double saliency, PointType& iCenter,
  MatrixType& R, unsigned int i, int threadId )
{
    if (i == 2)
        std::cout << "Yahoo" << std::endl;

  // Generate local field
//  OrientedTensorGeneratorPointer orientedTensor =
//    OrientedTensorGeneratorType::New();
//  orientedTensor->SetInput( m_VotingField[i] );
//  orientedTensor->SetCenter( iCenter );
//  orientedTensor->SetRotationMatrix( R );
//  orientedTensor->SetOutputSpacing( m_Output->GetSpacing() );
//  orientedTensor->SetOutputRegion( m_Region );
//  orientedTensor->Update();
//
//  ComputeVote( orientedTensor->GetOutput(), saliency, threadId );
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::ComputeVote( ImagePointer field, double saliency, int threadId )
{
  // Vote on the output
  RegionType region1, region2;
  OverlapRegion( field, m_ThreadImage[threadId], region2, region1 );

  IteratorType oIt( m_ThreadImage[threadId], region1 );
  IteratorType vIt( field, region2 );

  oIt.GoToBegin();
  vIt.GoToBegin();
  while( !vIt.IsAtEnd() )
  {
    if (saliency > 0)
       oIt.Set( oIt.Get() + vIt.Get()*saliency );
    else
       vIt.Set( oIt.Get() + vIt.Get() );

    ++vIt;
    ++oIt;
  }
}


// Compute the region R in A that overlaps with B
template< class TInputImage >
void
TensorVoting< TInputImage >
::OverlapRegion( ImagePointer A, ImagePointer B,
  RegionType& rA, RegionType& rB )
{
  SizeType sizeA, sizeB, s;
  sizeA = A->GetLargestPossibleRegion().GetSize();
  sizeB = B->GetLargestPossibleRegion().GetSize();

  IndexType sIndexA, sIndexB;
  IndexType tIndexA, tIndexB;

  A->TransformPhysicalPointToIndex( B->GetOrigin(), tIndexA );
  B->TransformPhysicalPointToIndex( A->GetOrigin(), tIndexB );

  VectorType p = A->GetOrigin() - B->GetOrigin();

  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    if ( p[i] > 0.0 )
    {
      sIndexA[i] = 0;
      sIndexB[i] = tIndexB[i];
      s[i] = sizeA[i];
      if ( s[i] > static_cast< SizeValueType >( sizeB[i] - sIndexB[i] - 1 ) )
      {
        s[i] = sizeB[i] - sIndexB[i];
      }
    }
    else
    {
      sIndexB[i] = 0;
      sIndexA[i] = tIndexA[i];
      s[i] = sizeB[i];
      if ( s[i] > static_cast< SizeValueType >(
        sizeA[i] - sIndexA[i] - 1 ) )
      {
        s[i] = sizeA[i] - sIndexA[i];
      }
    }
  }

  rA.SetIndex( sIndexA );
  rA.SetSize( s );
  rB.SetIndex( sIndexB );
  rB.SetSize( s );
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    InputImagePointer image =
      const_cast< InputImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template< class TInputImage >
void
TensorVoting< TInputImage >
::BeforeThreadedGenerateData()
{
  ImageConstPointer input = this->GetInput();
  RegionType region = input->GetLargestPossibleRegion();

  this->InitializeVotingFields();

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

  m_ThreadImage.resize( this->GetNumberOfThreads() );
  for( unsigned int i = 0; i < this->GetNumberOfThreads(); i++ )
  {
    m_ThreadImage[i] = 0;
  }
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::AfterThreadedGenerateData()
{
  for( unsigned int i = 0; i < this->GetNumberOfThreads(); i++ )
  {
    if ( m_ThreadImage[i] )
    {
      ComputeVote( this->m_Output, 0, i );
    }
  }
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::PadImage(const OutputImageRegionType& windowRegion,
  int threadId)
{
  IndexType start = windowRegion.GetIndex();
  SizeType size = windowRegion.GetSize();
  SpacingType spacing = this->GetInput()->GetSpacing();

  double radius = vcl_floor( vcl_sqrt( -vcl_log(0.01) *
    m_Sigma*m_Sigma ) );

  SizeValueType rad;
  RegionType region;
  IndexType index2;
  for(unsigned int i = 0; i< ImageDimension; i++)
  {
    rad = static_cast< SizeValueType >( radius/spacing[i] );
    index2[i] = 0;
    start[i] = start[i] - rad;
    size[i] = size[i] + 2 * rad;
  }
  region.SetIndex( index2 );
  region.SetSize( size );

  PointType origin;
  this->GetInput()->TransformIndexToPhysicalPoint( start, origin );

//   std::cout << region << ' ' << origin << std::endl;

  PixelType zeroTensor;
  zeroTensor.Fill( 0.0 );

  ImagePointer image = ImageType::New();
  image->SetRegions( region );
  image->SetOrigin( origin );
  image->SetSpacing( this->GetInput()->GetSpacing() );
  image->Allocate();
  image->FillBuffer( zeroTensor );

  m_ThreadImage[threadId] = image ;

}


template< class TInputImage >
void
TensorVoting< TInputImage >
::ThreadedGenerateData(const OutputImageRegionType& windowRegion,
  ThreadIdType threadId)
{
  PadImage(windowRegion, threadId);
  ImageConstPointer input = this->GetInput();

  // Iterate through the input image
  IndexType index;
  PixelType p;

  TokenIteratorType tIt( m_TokenImage, windowRegion );
  ConstIteratorType It( input, windowRegion );
  tIt.GoToBegin();
  It.GoToBegin();
  while( !It.IsAtEnd() )
  {
    if ( tIt.Get() )
    {
      index = It.GetIndex();
      p = It.Get();
      this->ComputeTensorVoting( index, p, static_cast<int>( threadId ) );
    }
    ++It;
    ++tIt;
  }
}


template< class TInputImage >
void
TensorVoting< TInputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
}


} // end namespace itk

#endif
