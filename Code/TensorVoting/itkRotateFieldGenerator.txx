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
#ifndef __itkRotateFieldGenerator_txx
#define __itkRotateFieldGenerator_txx

#include "itkRotateFieldGenerator.h"

namespace itk
{
/**
 * Constructor
 */
template<class TOutputImage>
RotateFieldGenerator<TOutputImage>
::RotateFieldGenerator()
{
  m_VotingField = 0;
}


template<class TOutputImage>
void
RotateFieldGenerator<TOutputImage>
::Update( ImagePointer image, MatrixType R, SpacingType space, 
  RegionType region, double weight )
{
  OrientedStickGeneratorPointer orientedStick = 
    OrientedStickGeneratorType::New();
  orientedStick->SetInput( this->m_VotingField );
  orientedStick->SetRotationMatrix( R );
  orientedStick->SetOutputSpacing( space );
  orientedStick->SetOutputRegion( region );
  orientedStick->Update();

  IteratorType It1( image, region );
  IteratorType It2( orientedStick->GetOutput(), region );
  It1.GoToBegin();
  It2.GoToBegin();
  while ( !It1.IsAtEnd() )
  {
    It1.Set( It1.Get() + It2.Get()*weight );
    ++It1;
    ++It2;
  }
}


template<class TOutputImage>
void
RotateFieldGenerator<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
