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
#ifndef __itkGPUTensorWithOrientation_txx
#define __itkGPUTensorWithOrientation_txx

#include "itkGPUTensorWithOrientation.h"
#include "itkNumericTraits.h"

namespace itk
{

/**
 * Constructor
 */
template<typename TInputImage, typename TOutputImage>
GPUTensorWithOrientation<TInputImage, TOutputImage>
::GPUTensorWithOrientation()
{
  m_OutputSpacing.Fill( 1.0 );
  m_Center.Fill( 0.0 );

  m_RotationMatrix.SetIdentity();
  m_RotationMatrixDefined = false;
  interpolator = InterpolatorType::New();

  m_GPUKernelManager = GPUKernelManager::New();
}


template<typename TInputImage, typename TOutputImage>
void
GPUTensorWithOrientation<TInputImage, TOutputImage>
::GPUGenerateData(void)
{
  GPUInputImagePointer  inPtr =  dynamic_cast< GPUInputImage * >( this->ProcessObject::GetInput(0) );
  GPUOutputImagePointer otPtr =  dynamic_cast< GPUOutputImage * >( this->ProcessObject::GetOutput(0) );

  ImageConstPointer votingField = this->GetInput();
  PointType origin = votingField->GetOrigin();

  this->GetOutput()->SetOrigin( origin );
  this->GetOutput()->SetSpacing( m_OutputSpacing );
  this->GetOutput()->SetDirection( votingField->GetDirection() );
  this->GetOutput()->SetLargestPossibleRegion( m_OutputRegion );

  // A zero tensor
  MatrixType ZeroTensor;
  ZeroTensor.Fill( 0 );

  OutputImagePointer m_Output = OutputImageType::New();
  m_Output->SetRegions( m_OutputRegion );
  m_Output->SetOrigin( origin );
  m_Output->SetSpacing( m_OutputSpacing );
  m_Output->Allocate();
  m_Output->FillBuffer( ZeroTensor );

  if ( !m_RotationMatrixDefined )
    {
    std::cerr << "Rotation matrix not provided as input." << std::endl;
    }

  // Determine rotation matrix
  MatrixType iRot;
  iRot = m_RotationMatrix.GetInverse();

  // Set the interpolator input
  interpolator->SetInputImage( votingField );

  // Iterate across all the pixels in the voting field
  PixelType p, q;
  PointType p1, p3;
  vnlVectorType p2;
  IndexType index;
  OutputIteratorType It( m_Output, m_OutputRegion );
  It.GoToBegin();
  while( !It.IsAtEnd() )
    {
    index = It.GetIndex();
    m_Output->TransformIndexToPhysicalPoint( index, p1 );

    p2 = m_RotationMatrix.GetVnlMatrix() * p1.GetVnlVector();

    for( unsigned int i = 0; i < ImageDimension; i++ )
      p3[i] = p2[i];

    if ( interpolator->IsInsideBuffer( p3 ) )
      {
      p = interpolator->Evaluate( p3 );
      q = iRot.GetVnlMatrix() * p.GetVnlMatrix() * iRot.GetVnlMatrix().transpose();
      It.Set( q );
      }
    else
      {
      p.Fill( 0.0 );
      It.Set( p );
      }
		++It;
	}

  for( unsigned int i = 0; i < ImageDimension; i++ )
    origin[i] += m_Center[i];

  m_Output->SetOrigin( origin );

  this->GraftOutput( m_Output );
}


template<typename TInputImage, typename TOutputImage>
void
GPUTensorWithOrientation<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  GPUSuperclass::PrintSelf(os,indent);
  os << indent << "Rotation Matrix: " << this->m_RotationMatrix << std::endl;
}

} // end namespace itk

#endif
