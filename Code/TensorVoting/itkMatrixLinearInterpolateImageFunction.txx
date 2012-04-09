/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMatrixLinearInterpolateImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2009-10-29 11:18:50 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMatrixLinearInterpolateImageFunction_txx
#define __itkMatrixLinearInterpolateImageFunction_txx

#include "itkMatrixLinearInterpolateImageFunction.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep>
const unsigned long
MatrixLinearInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep>
MatrixLinearInterpolateImageFunction< TInputImage, TCoordRep >
::MatrixLinearInterpolateImageFunction()
{}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void
MatrixLinearInterpolateImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename MatrixLinearInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
MatrixLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType& index) const
{
  unsigned int dim;  // index over dimension

  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  IndexType baseIndex;
  double distance[ImageDimension];

  for( dim = 0; dim < ImageDimension; dim++ )
    {
    baseIndex[dim] = static_cast< IndexValueType >( vcl_floor( index[dim] ) );
    distance[dim] = index[dim] - static_cast< double >( baseIndex[dim] );
    }

  /**
   * Interpolated value is the weighted sum of each of the surrounding
   * neighbors. The weight for each neighbor is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  OutputType output;
  output.Fill( 0.0 );

  double totalOverlap = 0.0;

  for( unsigned int counter = 0; counter < m_Neighbors; counter++ )
    {

    double overlap = 1.0;          // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour
    IndexType neighIndex;

    // get neighbor index and overlap fraction
    for( dim = 0; dim < ImageDimension; dim++ )
      {

      if ( upper & 1 )
        {
        neighIndex[dim] = baseIndex[dim] + 1;
// #ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
        // Take care of the case where the pixel is just
        // in the outer upper boundary of the image grid.
//         std::cout << "CenteredPixelCoordinates" << std::endl;
        if( neighIndex[dim] > this->m_EndIndex[dim] )
          {
          neighIndex[dim] = this->m_EndIndex[dim];
          }
// #endif
        overlap *= distance[dim];
        }
      else
        {
        neighIndex[dim] = baseIndex[dim];
// #ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
        // Take care of the case where the pixel is just
        // in the outer lower boundary of the image grid.
        if( neighIndex[dim] < this->m_StartIndex[dim] )
          {
          neighIndex[dim] = this->m_StartIndex[dim];
          }
// #endif
        overlap *= 1.0 - distance[dim];
        }

      upper >>= 1;

      }

    // get neighbor value only if overlap is not zero
    if( overlap )
      {
      const PixelType input = this->GetInputImage()->GetPixel( neighIndex );
      for(unsigned int k = 0; k < Dimension; k++ )
        {
        for(unsigned int l = 0; l < Dimension; l++ )
          {
          output[k][l] += overlap * static_cast<RealType>( input[k][l] );
          }
        }
      totalOverlap += overlap;
      }

    if( totalOverlap == 1.0 )
      {
      // finished
      break;
      }

    }

  return ( output );
}

} // end namespace itk

#endif
