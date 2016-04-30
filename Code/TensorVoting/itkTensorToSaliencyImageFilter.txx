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
#ifndef __itkTensorToSaliencyImageFilter_txx
#define __itkTensorToSaliencyImageFilter_txx

#include "itkTensorToSaliencyImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{

/**
 * Constructor
 */
template< class TInputImage, class TOutputImage >
TensorToSaliencyImageFilter< TInputImage, TOutputImage >
::TensorToSaliencyImageFilter()
{
  m_EigenCalculator.SetDimension( ImageDimension );

  m_EigenMatrix = 0;
  m_ComputeEigenMatrix = false;

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TOutputImage::New() );
}


template< class TInputImage, class TOutputImage >
void
TensorToSaliencyImageFilter< TInputImage, TOutputImage >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template< class TInputImage, class TOutputImage >
void
TensorToSaliencyImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    ImagePointer image = const_cast< ImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}


template< class TInputImage, class TOutputImage >
void
TensorToSaliencyImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  ImageConstPointer input = this->GetInput();
  RegionType region       = input->GetLargestPossibleRegion();
  SpacingType spacing     = input->GetSpacing();
  PointType origin        = input->GetOrigin();

  // A zero vector and tensor
  VectorType ZeroVector;
  ZeroVector.Fill( 0.0 );

  MatrixType ZeroMatrix;
  ZeroMatrix.Fill( 0.0 );

  // Allocate the output image
  OutputImagePointer m_Output = this->GetOutput();
  m_Output->SetRegions( region );
  m_Output->SetOrigin( origin );
  m_Output->SetSpacing( spacing );
  m_Output->Allocate();
  m_Output->FillBuffer( ZeroVector );

  if ( m_ComputeEigenMatrix )
    {
    m_EigenMatrix = MatrixImageType::New();
    m_EigenMatrix->SetRegions( region );
    m_EigenMatrix->SetOrigin( origin );
    m_EigenMatrix->SetSpacing( spacing );
    m_EigenMatrix->Allocate();
    m_EigenMatrix->FillBuffer( ZeroMatrix );\
    }
}


template< class TInputImage, class TOutputImage >
void
TensorToSaliencyImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const RegionType& windowRegion, ThreadIdType threadId)
{
  ImageConstPointer input = this->GetInput();
  OutputImagePointer m_Output = this->GetOutput();

  // Iterate through the input image
  MatrixType eigenMatrix;
  VectorType eigenVector;
  ConstIteratorType It( input, windowRegion );
  OutputIteratorType oIt( m_Output, windowRegion );

  if ( m_ComputeEigenMatrix )
  {
    MatrixIteratorType mIt( m_EigenMatrix, windowRegion );
    mIt.GoToBegin();
    It.GoToBegin();
    oIt.GoToBegin();
    while( !It.IsAtEnd() )
    {
      // Compute eigen values and eigen matrix
      m_EigenCalculator.SetOrderEigenValues( 1 );
      m_EigenCalculator.ComputeEigenValuesAndVectors( It.Get(), eigenVector, eigenMatrix );
      mIt.Set( eigenMatrix );

      for( unsigned int i = ImageDimension-1; i > 0; i-- )
        eigenVector[i] -= eigenVector[i-1];

      oIt.Set( eigenVector );
      ++It;
      ++oIt;
      ++mIt;
    }
  }
  else
  {
    It.GoToBegin();
    oIt.GoToBegin();
    while( !It.IsAtEnd() )
    {
      // Compute eigen values and eigen matrix
      m_EigenCalculator.SetOrderEigenValues( 1 );
      m_EigenCalculator.ComputeEigenValuesAndVectors( It.Get(), eigenVector, eigenMatrix );

      for( unsigned int i = ImageDimension-1; i > 0; i-- )
        eigenVector[i] -= eigenVector[i-1];

      oIt.Set( eigenVector );
      ++It;
      ++oIt;
    }
  }
}


template< class TInputImage, class TOutputImage >
void
TensorToSaliencyImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
