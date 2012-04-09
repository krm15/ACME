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
#ifndef __itkStickFieldGenerator2D_txx
#define __itkStickFieldGenerator2D_txx

#include "itkStickFieldGenerator2D.h"

namespace itk
{

/**
 * Constructor
 */
template<class TOutputImage>
StickFieldGenerator2D<TOutputImage>
::StickFieldGenerator2D()
{}


template<class TOutputImage>
void
StickFieldGenerator2D<TOutputImage>
::ComputeStickField(void)
{
  this->AllocateOutput();
  double sigma = this->m_Sigma;
  // Have to take care of c with m_Spacing
//   double c = (-16 * vcl_log(0.1)*(m_Sigma-1) * 
//     vnl_math::one_over_pi * vnl_math::one_over_pi)/ vcl_log(2);
  double c = 0.25;

  double theta, l, s, k, d;
  PointType pt;
  MatrixType p;
  IndexType index;

  // Iterate across all the pixels
  IteratorType It( this->m_Output, 
    this->m_Output->GetLargestPossibleRegion() );
  It.GoToBegin();
  while( !It.IsAtEnd() )
    {
    index = It.GetIndex();
    this->m_Output->TransformIndexToPhysicalPoint( index,pt );

    if (pt[0] != 0)
      {
      theta = vcl_atan( pt[1]/pt[0] );
      }
    else
      {
      if (pt[1] == 0)
        theta = 0;
      else
        theta = vnl_math::pi_over_2;
      }

    theta = vcl_abs( theta );

    p[0][0] = vcl_sin( 2*theta ) * vcl_sin( 2*theta );
    p[0][1] = vcl_sin( 2*theta ) * vcl_cos( 2*theta );
    p[1][0] = p[0][1];
    p[1][1] = vcl_cos( 2*theta ) * vcl_cos( 2*theta );

    if ( pt[0]*pt[1] >= 0 )
      {
      p[1][0] *= -1;
      p[0][1] *= -1;
      }

    l = vcl_sqrt( pt[0]*pt[0] + pt[1]*pt[1] );

    if (  l == 0 )
      {
      s = 0;
      k = 0;
      }
    else
      {
      if ( theta == 0 )
        {
        s = l;
        }
      else
        {
        s = theta*l/vcl_sin(theta);
        }
      k = 2*vcl_sin(theta)/l;
      }


    double x = ( static_cast<double>( l ) - 0.5*sigma ) / 2;
    double e = 1.0 / ( 1.0 + vcl_exp(- x ) );
    k *= e;

    d = vcl_exp(-(s*s + c*k*k) / (sigma*sigma) );

    if ( theta > 0.5*vnl_math::pi_over_2 )
      {
      d = 0;
      }


    // Eliminate out of radius artifacts caused by rotation
    double radius = vcl_floor( vcl_sqrt( -vcl_log(0.01) * 
      sigma * sigma ) );
    if ( l >= radius )
    {
      d = 0;
    }

    p[0][0] *= d;
    p[0][1] *= d;
    p[1][0] *= d;
    p[1][1] *= d;

    It.Set( p );
    ++It;
    }
}


template<class TOutputImage>
void
StickFieldGenerator2D<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
