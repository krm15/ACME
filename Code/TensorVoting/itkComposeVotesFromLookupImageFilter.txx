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

#ifndef __itkComposeVotesFromLookupImageFilter_txx
#define __itkComposeVotesFromLookupImageFilter_txx

#include "itkComposeVotesFromLookupImageFilter.h"

namespace itk
{

/**
 * Constructor
 */
template< class TInputImage >
ComposeVotesFromLookupImageFilter<TInputImage >
::ComposeVotesFromLookupImageFilter()
{
  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );
  this->Superclass::SetNthOutput ( 0, TInputImage::New() );
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
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
ComposeVotesFromLookupImageFilter< TInputImage >
::BeforeThreadedGenerateData()
{
  m_ThreadImage.resize( this->GetNumberOfThreads() );
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::AfterThreadedGenerateData()
{  
  for( unsigned int i = 0; i < this->GetNumberOfThreads(); i++ )
    {
    if ( m_ThreadImage[i] )
      {
      OutputIteratorType tIt( m_ThreadImage[i], m_ThreadImage[i]->GetLargestPossibleRegion() );
      OutputIteratorType oIt( m_Output, m_Output->GetLargestPossibleRegion() );

      oIt.GoToBegin();
      tIt.GoToBegin();
      while( !tIt.IsAtEnd() )
        {
        oIt.Set( oIt.Get() + tIt.Get() );
        ++tIt;
        ++oIt;
        }
      }
    }
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::ComputeVote( double saliency, unsigned int threadId, OutputImagePointer orientedVotingField )
{
  // Vote on the output
  RegionType region1, region2;
  OverlapRegion( orientedVotingField, m_ThreadImage[threadId], region1, region2 );

  OutputIteratorType vIt( orientedVotingField, region1 );
  OutputIteratorType oIt( m_ThreadImage[threadId], region2 );

  oIt.GoToBegin();
  vIt.GoToBegin();
  while( !vIt.IsAtEnd() )
  {
    oIt.Set( oIt.Get() + vIt.Get()*saliency );
    ++vIt;
    ++oIt;
  }
}


// Compute the region R in A that overlaps with B
template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::OverlapRegion( OutputImagePointer A, OutputImagePointer B,
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
ComposeVotesFromLookupImageFilter< TInputImage >
::ThreadedGenerateData(const RegionType& windowRegion,
  ThreadIdType threadId)
{
  // A zero tensor
  MatrixType ZeroTensor;
  ZeroTensor.Fill( 0 );
  
  m_ThreadImage[threadId] = OutputImageType::New();
  m_ThreadImage[threadId]->SetRegions( m_Output->GetLargestPossibleRegion() );
  m_ThreadImage[threadId]->SetSpacing( m_Output->GetSpacing() );
  m_ThreadImage[threadId]->CopyInformation( m_Output );
  m_ThreadImage[threadId]->Allocate();
  m_ThreadImage[threadId]->FillBuffer( ZeroTensor );
  
  InputImageConstPointer input = this->GetInput();

  IndexType index, index2;
  VectorType u;
  double stickSaliency;
  MatrixType R;
  PointType pt, origin, origin2;
  IdType ll;

  // Iterate through the list
  ConstInputIteratorType iIt( input, windowRegion );
  iIt.GoToBegin();
  while( !iIt.IsAtEnd() )
    {
    ll = iIt.Get();
    if ( ll.size() > 0 )
      {
      index2 = iIt.GetIndex();
      input->TransformIndexToPhysicalPoint( index2, pt );
      for(unsigned int i = 0; i < ImageDimension; i++)
        {
        u[i] = pt[i];
        // cos(theta)cos(phi), cos(theta)sin(phi), sin(theta)
        }

      rMatrixHelper.ComputeRotationMatrix( u, R );

      // Compute the rotated voting field
      OutputImagePointer orientedVotingField;
//      ComputeOrientedField( R, 0, orientedVotingField );

      OrientedTensorGeneratorPointer orientedTensor = OrientedTensorGeneratorType::New();
      orientedTensor->SetInput( m_VotingField );
      orientedTensor->SetRotationMatrix( R );
      orientedTensor->SetOutputSpacing( m_Output->GetSpacing() );
      orientedTensor->SetOutputRegion( m_VotingField->GetLargestPossibleRegion() );
      orientedTensor->Update();
      orientedVotingField = orientedTensor->GetOutput();

      origin = orientedVotingField->GetOrigin();

      // Iterate the list,
      for (typename IdType::const_iterator iter = ll.begin(); iter != ll.end(); ++iter)
        {
        // Update the origin of the voting field
        index = *iter;
        m_Output->TransformIndexToPhysicalPoint( index, pt );
        for( unsigned int i = 0; i < ImageDimension; i++ )
          origin2[i] = origin[i] + pt[i];

        orientedVotingField->SetOrigin( origin2 );
        stickSaliency = m_StickSaliencyImage->GetPixel( index );
        // Add to the output with salience
        ComputeVote( stickSaliency, threadId, orientedVotingField );
        }
      }
    ++iIt;
    }    
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}


} // end namespace itk

#endif
