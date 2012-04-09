/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStickFieldGenerator3D.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 23:24:23 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkStickFieldGenerator3D_txx
#define __itkStickFieldGenerator3D_txx

#include "itkStickFieldGenerator3D.h"

namespace itk
{

/**
 * Constructor
 */
template<class TOutputImage>
StickFieldGenerator3D<TOutputImage>
::StickFieldGenerator3D()
{
  interpolator = InterpolatorType::New();
}


template<class TOutputImage>
double
StickFieldGenerator3D<TOutputImage>
::FindTheta( double x, double y )
{
  double theta;

  if ( x == 0 )
  {
    if ( y == 0 )
    {
      theta = 0;
    }
    else
    {
      theta = vnl_math::pi_over_2;
    }
  }
  else
  {
    theta = vcl_atan(  vcl_abs( y / x ) );
  }

  // Determine the quadrant of theta
  if ( x > 0 )
  {
    if ( y < 0 )
    {
      theta = - theta;
    }
  }
  else
  {
    if ( y < 0 )
    {
      theta += vnl_math::pi;
    }
    else
    {
      theta = vnl_math::pi - theta;
    }
  }

  return theta;
}


template<class TOutputImage>
typename StickFieldGenerator3D<TOutputImage>::
MatrixType
StickFieldGenerator3D<TOutputImage>
::RotationMatrixAboutZ( double theta )
{
  MatrixType Rot;

  Rot[0][0] = vcl_cos(theta);
  Rot[0][1] = -vcl_sin(theta);
  Rot[0][2] = 0;
  Rot[1][0] = vcl_sin(theta);
  Rot[1][1] = vcl_cos(theta);
  Rot[1][2] = 0;
  Rot[2][0] = 0;
  Rot[2][1] = 0;
  Rot[2][2] = 1;

  return Rot;
}


template<class TOutputImage>
void
StickFieldGenerator3D<TOutputImage>
::ComputeStickField(void)
{
  this->AllocateOutput();

  // Compute stick at a highest resolution and use that
  // to interpolate the ball field
  double sp = NumericTraits<double>::max();
  Spacing2DType spacing;
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    if (sp > static_cast< double >( this->m_Spacing[i] ) )
    {
       sp = this->m_Spacing[i];
    }
  }

  for( unsigned int i = 0; i < ImageDimension-1; i++ )
    spacing[i] = sp;

  Stick2DGeneratorPointer stick = Stick2DGeneratorType::New();
  stick->SetSigma ( this->m_Sigma );
  stick->SetSpacing( spacing );
  stick->ComputeStickField();
  Image2DPointer stick2DField = stick->GetOutputImage();
  Region2DType region = stick2DField->GetLargestPossibleRegion();

  // Populate the 3D field fully
  IteratorType It( this->m_Output, 
    this->m_Output->GetLargestPossibleRegion() );

  // Set the interpolator input
  interpolator->SetInputImage( stick2DField );

  // temporary variables in 3D and 2D
  double theta;
  IndexType index;
  PointType p;
  MatrixType Rot, iRot, pixel;

  Point2DType p2D;
  Matrix2DType pixel2D;
  pixel2D.Fill( 0 );

  It.GoToBegin();
  while( !It.IsAtEnd() )
    {
    index = It.GetIndex();
    this->m_Output->TransformIndexToPhysicalPoint( index, p );

	  //Project onto xz plane (y=0)
	  theta = FindTheta( p[0], p[1] );

    // Rotation matrix that takes the output coordinates 
    // to 2D coordinates
    Rot = RotationMatrixAboutZ( -theta );

    // Find point in 2D plane
 		p2D[0] = Rot[0][0]*p[0]  + Rot[0][1]*p[1];
		p2D[1] = p[2];

    // interpolation here
    if ( interpolator->IsInsideBuffer( p2D ) )
      {
      pixel2D = interpolator->Evaluate( p2D );
      pixel.Fill( 0 );
      pixel[0][0] = pixel2D[0][0];
      pixel[0][2] = pixel2D[0][1];
      pixel[2][0] = pixel2D[1][0];
      pixel[2][2] = pixel2D[1][1];

      //Rotate the tensor by inverse(rot)
      iRot = Rot.GetInverse();

      pixel = iRot.GetVnlMatrix() * pixel.GetVnlMatrix() * 
        Rot.GetVnlMatrix();
      It.Set( pixel );
      }
    ++It;
    }
}


template<class TOutputImage>
void
StickFieldGenerator3D<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
