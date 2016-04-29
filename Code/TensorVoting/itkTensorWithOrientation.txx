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
#ifndef __itkTensorWithOrientation_txx
#define __itkTensorWithOrientation_txx

#include "itkTensorWithOrientation.h"
#include "itkNumericTraits.h"

namespace itk
{

/**
 * Constructor
 */
template<class TInputImage, class TOutputImage>
TensorWithOrientation<TInputImage, TOutputImage>
::TensorWithOrientation()
{
  m_OutputSpacing.Fill( 1.0 );
  m_Center.Fill( 0.0 );

  m_RotationMatrix.SetIdentity();
  m_RotationMatrixDefined = false;
  interpolator = InterpolatorType::New();
}


template<class TInputImage, class TOutputImage>
void
TensorWithOrientation<TInputImage, TOutputImage>
::GenerateData(void)
{
  ImageConstPointer votingField = this->GetInput();
  PointType origin = votingField->GetOrigin();

  this->GetOutput()->SetOrigin( origin );
  this->GetOutput()->SetSpacing( m_OutputSpacing );
  this->GetOutput()->SetDirection( votingField->GetDirection() );
  this->GetOutput()->SetLargestPossibleRegion( m_OutputRegion );

  // A zero tensor
  MatrixType ZeroTensor;
  ZeroTensor.Fill( 0 );

  ImagePointer m_Output = ImageType::New();
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


template<class TInputImage, class TOutputImage>
void
TensorWithOrientation<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Rotation Matrix: " << this->m_RotationMatrix << std::endl;
}

} // end namespace itk

#endif
