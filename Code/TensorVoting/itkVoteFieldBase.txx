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
#ifndef __itkVoteFieldBase_txx
#define __itkVoteFieldBase_txx

#include "itkVoteFieldBase.h"

namespace itk
{
/**
 * Constructor
 */
template<class TOutputImage>
VoteFieldBase<TOutputImage>
::VoteFieldBase()
{
  m_Sigma = 5.0;
  m_Output = ImageType::New();

  m_Origin.Fill( 0.0 );
  m_Spacing.Fill( 1.0 );
}


template<class TOutputImage>
void
VoteFieldBase<TOutputImage>
::AllocateOutput(void)
{
  //Compute stick parameters
  double radius = vcl_floor( vcl_sqrt( -vcl_log(0.01) * 
    m_Sigma*m_Sigma ) );

  SizeValueType rad;
  SizeType size;
  IndexType index;
  for( unsigned int i = 0; i< ImageDimension; i++ )
    {
    rad = static_cast< SizeValueType >( radius/m_Spacing[i] );
    size[i] = 2*rad + 1;
    index[i] = 0;
    m_Origin[i] = -static_cast< double >( rad * m_Spacing[i] );
    }

  // Set the region
  RegionType region;
  region.SetSize( size );
  region.SetIndex( index );

  // A zero tensor
  MatrixType ZeroTensor;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      ZeroTensor[i][j] = 0;
      }

  // Initialize the output image
  m_Output->SetRegions( region );
  m_Output->SetOrigin( m_Origin );
  m_Output->SetSpacing( m_Spacing );
  m_Output->Allocate();
  m_Output->FillBuffer( ZeroTensor );
}


template<class TOutputImage>
void
VoteFieldBase<TOutputImage>
::Postprocess(ImagePointer image, double& count)
{
  // Normalization
  RegionType region = image->GetLargestPossibleRegion();
  IteratorType It1( image, region );
  It1.GoToBegin();
  while ( !It1.IsAtEnd() )
  {
    It1.Set( It1.Get()/count );
    ++It1;
  }
}


template<class TOutputImage>
void
VoteFieldBase<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
}

} // end namespace itk

#endif
