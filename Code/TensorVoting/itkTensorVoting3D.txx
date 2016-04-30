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

#ifndef __itkTensorVoting3D_txx
#define __itkTensorVoting3D_txx

#include "itkTensorVoting3D.h"

namespace itk
{

/**
 * Constructor
 */
template< class TInputImage >
TensorVoting3D<TInputImage >
::TensorVoting3D()
{
  m_Sigma = 5.0;
  m_UseSparseVoting = 0;
  m_OrientedVotingField = ITK_NULLPTR;

  IdType emptyList;
  emptyList.clear();

  SpacingType nSpacing;
  nSpacing[0] = nSpacing[1] = nSpacing[2] = 0.01;

  PointType nOrigin;
  nOrigin[0] = nOrigin[1] = -1.0;
  nOrigin[2] = 0;

  SizeType nSize;
  nSize[0] = nSize[1] = 201;
  nSize[2] = 101;

  IndexType nIndex;
  nIndex[0] = nIndex[1] = nIndex[2] = 0;

  RegionType nRegion;
  nRegion.SetSize( nSize );
  nRegion.SetIndex( nIndex );

  m_Lookup = InternalImageType::New();
  m_Lookup->SetRegions( nRegion );
  m_Lookup->SetOrigin( nOrigin );
  m_Lookup->SetSpacing( nSpacing );
  m_Lookup->Allocate();
  m_Lookup->FillBuffer( emptyList );

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TInputImage::New() );
}

template< class TInputImage >
void
TensorVoting3D< TInputImage >
::InitializeVotingFields(void)
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
    spacing[i] = sp;

  // Compute the stick and ball voting fields
  StickGeneratorPointer stick = StickGeneratorType::New();
  stick->SetSigma ( this->m_Sigma );
  stick->SetSpacing( spacing );
  stick->ComputeStickField();
  this->m_VotingField.push_back( stick->GetOutputImage() );

  PlateGeneratorPointer plate = PlateGeneratorType::New();
  plate->SetSigma ( this->m_Sigma );
  plate->SetInput( stick->GetOutputImage() );
  plate->SetSpacing( spacing );
  plate->ComputePlateField();
  this->m_VotingField.push_back( plate->GetOutputImage() );

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

// Compute the region R in A that overlaps with B
template< class TInputImage >
void
TensorVoting3D< TInputImage >
::OverlapRegion( InputImagePointer A, InputImagePointer B,
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
TensorVoting3D< TInputImage >
::ComputeVote( double saliency )
{
  // Vote on the output
  RegionType region1, region2;
  OverlapRegion( m_OrientedVotingField, m_Output, region1, region2 );

  IteratorType vIt( m_OrientedVotingField, region1 );
  IteratorType oIt( m_Output, region2 );

  oIt.GoToBegin();
  vIt.GoToBegin();
  while( !vIt.IsAtEnd() )
  {
    oIt.Set( oIt.Get() + vIt.Get()*saliency );
    ++vIt;
    ++oIt;
  }
}


template< class TInputImage >
void
TensorVoting3D< TInputImage >
::ComputeOrientedField( MatrixType& R, unsigned int i )
{
  // Generate local field
  OrientedTensorGeneratorPointer orientedTensor = OrientedTensorGeneratorType::New();
  orientedTensor->SetInput( m_VotingField[i] );
  orientedTensor->SetRotationMatrix( R );
  orientedTensor->SetOutputSpacing( m_Output->GetSpacing() );
  orientedTensor->SetOutputRegion( m_Region );
  orientedTensor->Update();
  this->m_OrientedVotingField = orientedTensor->GetOutput();
  this->m_OrientedVotingField->DisconnectPipeline();
}


template< class TInputImage >
void
TensorVoting3D< TInputImage >
::ComputeLookup()
{
  IndexType index, index2;
  PixelType p;
  VectorType u, v;
  double stickSaliency;
  MatrixType eigenMatrix, R;
  PointType pt;
  IdType ll;

  RegionType region = this->GetInput()->GetLargestPossibleRegion();

  typename SaliencyFilterType::Pointer saliencyFilter = SaliencyFilterType::New();
  saliencyFilter->SetInput( this->GetInput() );
  saliencyFilter->SetComputeEigenMatrix( 1 );
  saliencyFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  saliencyFilter->Update();
  typename VectorImageType::Pointer saliencyImage = saliencyFilter->GetOutput();
  InputImagePointer eigenMatrixImage = saliencyFilter->GetEigenMatrix();

  typename IndexFilterType::Pointer componentExtractor = IndexFilterType::New();
  componentExtractor->SetInput( saliencyImage );
  componentExtractor->SetIndex( 2 );
  componentExtractor->Update();
  this->m_StickSaliencyImage = componentExtractor->GetOutput();

  // Compute an image of lists with similar pixel types
  ConstIteratorType It( this->GetInput(), region );
  DoubleIteratorType sIt( m_StickSaliencyImage, region );
  IteratorType eIt( eigenMatrixImage, region );

  It.GoToBegin();
  sIt.GoToBegin();
  eIt.GoToBegin();
  while( !It.IsAtEnd() )
    {
    index = It.GetIndex();
    p = It.Get();
    stickSaliency = sIt.Get();
    eigenMatrix = eIt.Get();

    if ( stickSaliency > 0.001 )
      {
      // Compute orientation (theta in 2D)
      u = eigenMatrix[ImageDimension-1];
      if ( u[2] < 0 )
        u = -u;

      for(unsigned int i = 0; i < ImageDimension; i++)
          pt[i] = u[i];

      m_Lookup->TransformPhysicalPointToIndex( pt, index2 );
      ll = m_Lookup->GetPixel( index2 );
      ll.push_back( index );
      m_Lookup->SetPixel( index2, ll );
      }

    ++It;
    ++sIt;
    ++eIt;
    }
}

template< class TInputImage >
void
TensorVoting3D< TInputImage >
::ComposeOutput()
{
  IndexType index, index2;
  VectorType u;
  double stickSaliency;
  MatrixType R;
  PointType pt, origin, origin2;
  IdType ll;

  // Iterate through the list
  InternalIteratorType iIt( m_Lookup, m_Lookup->GetLargestPossibleRegion() );
  iIt.GoToBegin();
  while( !iIt.IsAtEnd() )
    {
    ll = iIt.Get();
    if ( ll.size() > 0 )
      {
      index2 = iIt.GetIndex();
      m_Lookup->TransformIndexToPhysicalPoint( index2, pt );
      for(unsigned int i = 0; i < ImageDimension; i++)
        u[i] = pt[i];
      rMatrixHelper.ComputeRotationMatrix( u, R );

      // Compute the rotated voting field
      this->ComputeOrientedField( R, 0 );
      origin = m_OrientedVotingField->GetOrigin();

      // Iterate the list,
      for (typename IdType::const_iterator iter = ll.begin(); iter != ll.end(); ++iter)
        {
        // Update the origin of the voting field
        index = *iter;
        m_Output->TransformIndexToPhysicalPoint( index, pt );

        for( unsigned int i = 0; i < ImageDimension; i++ )
          origin2[i] = origin[i] + pt[i];

        m_OrientedVotingField->SetOrigin( origin2 );
        stickSaliency = m_StickSaliencyImage->GetPixel( index );
        // Add to the output with salience
        ComputeVote( stickSaliency );
        }
      }
    ++iIt;
    }
}


template< class TInputImage >
void
TensorVoting3D< TInputImage >
::GenerateData()
{
  ImageConstPointer input = this->GetInput();
  RegionType region = input->GetLargestPossibleRegion();

  // A zero tensor
  MatrixType ZeroTensor;
  ZeroTensor.Fill( 0 );

  // Allocate the output image
  m_Output = InputImageType::New();
  m_Output->SetRegions( region );
  m_Output->CopyInformation( input );
  m_Output->Allocate();
  m_Output->FillBuffer( ZeroTensor );

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

  itk::TimeProbe cputimer1, cputimer2;
  cputimer1.Start();

  // Fill Lookup with lists of all similar voxels
  ComputeLookup();

  cputimer1.Stop();
  std::cout << "ComputeLookup took " << cputimer1.GetMean() << " seconds" << std::endl;

  cputimer2.Start();

  // Compose the voting fields
  ComposeOutput();

  cputimer2.Stop();
  std::cout << "ComposeOutput took " << cputimer2.GetMean() << " seconds" << std::endl;

  // Copy the padded output image to output
  this->GraftOutput( m_Output );

  std::cout << "Grafting complete..." << std::endl;
}

template< class TInputImage >
void
TensorVoting3D< TInputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
