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
  this->Superclass::SetNthOutput ( 0, InputImageType::New() );
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion( output );
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
  
  const InputImageType *outputPtr = this->GetOutput();
  const ImageRegionSplitterBase * splitter = this->GetImageRegionSplitter();
  m_ValidThreads = splitter->GetNumberOfSplits( outputPtr->GetRequestedRegion(), this->GetNumberOfThreads() );
  m_ThreadImage.resize( m_ValidThreads );
}


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::AfterThreadedGenerateData()
{
  AddNonScalarFilterPointer addFilter = AddNonScalarFilterType::New();
  addFilter->SetInput( m_ThreadImage[0] );
  addFilter->SetImageVector( m_ThreadImage );
  addFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  addFilter->Update();

  OutputIteratorType tIt( addFilter->GetOutput(), m_OutputImage->GetLargestPossibleRegion() );
  OutputIteratorType oIt( m_OutputImage, m_OutputImage->GetLargestPossibleRegion() );
  oIt.GoToBegin();
  tIt.GoToBegin();
  while( !tIt.IsAtEnd() )
  {
    oIt.Set( oIt.Get() + tIt.Get() );
    ++tIt;
    ++oIt;
  }
  
  for( unsigned int i = 0; i < m_ValidThreads; i++ )
    {
    m_ThreadImage[i] = ITK_NULLPTR;
    }
  m_ThreadImage.clear();    
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
  m_ThreadImage[threadId]->SetRegions( m_OutputImage->GetLargestPossibleRegion() );
  m_ThreadImage[threadId]->CopyInformation( m_OutputImage );
  m_ThreadImage[threadId]->Allocate();
  m_ThreadImage[threadId]->FillBuffer( ZeroTensor );
  
  InputImageConstPointer input = this->GetInput();

  IntegrateBallVotingField( threadId );
  
  IndexType index, index2;
  double saliency;
  MatrixType R;
  PointType pt_sph, origin, origin2;
  IdType ll;
  OutputImagePointer orientedVotingField;
  double theta;

  ConstInputIteratorType iIt( input, input->GetLargestPossibleRegion() );
  iIt.GoToBegin();
  unsigned int counter = 0;
  while( !iIt.IsAtEnd() )
  {
    index2 = iIt.GetIndex();
    if ( counter%m_ValidThreads == static_cast<unsigned int>(threadId) )
    {
      for( unsigned int j = 1; j < ImageDimension; j++ )
      {
        ll = iIt.Get()[j];
        if ( ll.size() > 0 )
        {
          input->TransformIndexToPhysicalPoint( index2, pt_sph );
          ComputeRotationMatrix( R, pt_sph );

          // Compute the rotated voting field
          OrientedTensorGeneratorPointer orientedTensor = OrientedTensorGeneratorType::New();
          orientedTensor->SetInput( m_VotingField[j] );
          orientedTensor->SetRotationMatrix( R );
          orientedTensor->SetOutputSpacing( m_OutputImage->GetSpacing() );
          orientedTensor->SetOutputRegion( m_VotingField[j]->GetLargestPossibleRegion() );
          orientedTensor->Update();
          orientedVotingField = orientedTensor->GetOutput();
          orientedVotingField->DisconnectPipeline();

          origin = orientedVotingField->GetOrigin();

          // Iterate the list,
          for (typename IdType::const_iterator iter = ll.begin();
               iter != ll.end(); ++iter)
          {
            // Update the origin of the voting field
            index = *iter;
            m_OutputImage->TransformIndexToPhysicalPoint( index, pt_sph );
            for( unsigned int i = 0; i < ImageDimension; i++ )
            {
              origin2[i] = origin[i] + pt_sph[i];
            }
            orientedVotingField->SetOrigin( origin2 );
            saliency = m_SaliencyImage->GetPixel( index )[j];

            // Add to the output with saliency
            ComputeVote( saliency, threadId, orientedVotingField );
          }
        }
      }
    }
    ++iIt;
    counter++;
  }
} 


template< class TInputImage >
void
ComposeVotesFromLookupImageFilter< TInputImage >
::IntegrateBallVotingField( ThreadIdType threadId )
{
  IndexType index, index2;
  VectorType u;
  double saliency;

  PointType pt_cart, pt_sph, origin, origin2;
  IdType ll;

  MatrixType R;
  R.SetIdentity();

  InputImageConstPointer input = this->GetInput();

  // Compute the rotated voting field
  OrientedTensorGeneratorPointer orientedTensor = OrientedTensorGeneratorType::New();
  orientedTensor->SetInput( m_VotingField[0] );
  orientedTensor->SetRotationMatrix( R );
  orientedTensor->SetOutputSpacing( m_OutputImage->GetSpacing() );
  orientedTensor->SetOutputRegion( m_VotingField[0]->GetLargestPossibleRegion() );
  orientedTensor->Update();
  OutputImagePointer orientedVotingField = orientedTensor->GetOutput();
  orientedVotingField->DisconnectPipeline();
  origin = orientedVotingField->GetOrigin();

  ConstInputIteratorType iIt( input, input->GetLargestPossibleRegion() );
  iIt.GoToBegin();
  unsigned int counter = 0;

  while( !iIt.IsAtEnd() )
  {
    index2 = iIt.GetIndex();
    ll = iIt.Get()[0];
    if ( ( counter%m_ValidThreads == static_cast<unsigned int>(threadId) ) &&
         ( ll.size() > 0 ) )
    {
      // Iterate the list,
      for (typename IdType::const_iterator iter = ll.begin(); iter != ll.end(); ++iter)
      {
        // Update the origin of the voting field
        index = *iter;
        m_OutputImage->TransformIndexToPhysicalPoint( index, pt_cart );
        for( unsigned int i = 0; i < ImageDimension; i++ )
        {
          origin2[i] = origin[i] + pt_cart[i];
        }
        orientedVotingField->SetOrigin( origin2 );
        saliency = m_SaliencyImage->GetPixel( index )[0];

        // Add to the output with saliency
        ComputeVote( saliency, threadId, orientedVotingField );
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


