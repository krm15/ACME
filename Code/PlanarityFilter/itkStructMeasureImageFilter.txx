/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStructMeasureImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 20:59:44 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkStructMeasureImageFilter_txx
#define __itkStructMeasureImageFilter_txx

#include "itkStructMeasureImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

#define EPSILON  1e-03

namespace itk
{
template < typename TPixel >
StructMeasureImageFilter< TPixel >
::StructMeasureImageFilter()
{
  m_Alpha = 0.5;
  m_Beta  = 0.5;
  m_Gamma = 8.0;
  m_C     = 10e-6;
  m_ObjectType = 0;

  m_ComputeEigenVectors = true;
  m_ScaleMembranenessMeasure  = false;//true
  m_StructureIsPositive  = true;
}


template < typename TPixel >
double
StructMeasureImageFilter< TPixel >::
GetStructureMeasure( double& Lambda1, double& Lambda2, double& Lambda3 )
{
  double A, B, S;
  double vesMeasure_1 = 0, vesMeasure_2 = 0;
  double vesMeasure_3 = 0, vesMeasure_4 = 0;
  double Lambda1Abs = vnl_math_abs( Lambda1 );
  double Lambda2Abs = vnl_math_abs( Lambda2 );
  double Lambda3Abs = vnl_math_abs( Lambda3 );

  double Lambda1Sqr = vnl_math_sqr( Lambda1 );
  double Lambda2Sqr = vnl_math_sqr( Lambda2 );
  double Lambda3Sqr = vnl_math_sqr( Lambda3 );

  double AlphaSqr = vnl_math_sqr( m_Alpha );
  double BetaSqr = vnl_math_sqr( m_Beta );
  double GammaSqr = vnl_math_sqr( m_Gamma );

  if ( m_ObjectType == 0 ) // plate
  {
    A  = Lambda2Abs / Lambda3Abs;//plate vs vessel ~ ball
    B  = vcl_sqrt ( vnl_math_abs( Lambda1 * Lambda2 )) / ( Lambda3Abs );//vessel vs plate vs ball
    S  = vcl_sqrt( Lambda1Sqr + Lambda2Sqr + Lambda3Sqr );

    vesMeasure_1  =
      ( vcl_exp(-0.5*(( vnl_math_sqr(A) ) / ( AlphaSqr ) )));

    vesMeasure_2  =
      vcl_exp ( -0.5 * ((vnl_math_sqr( B )) /  ( BetaSqr)));

    vesMeasure_3  =
      ( 1 - vcl_exp( -1.0 * (( vnl_math_sqr( S )) / ( 2.0 * ( GammaSqr)))));

    vesMeasure_4  =
      vcl_exp ( -1.0 * ( 2.0 * vnl_math_sqr( m_C )) / ( Lambda3Sqr) );
  }

  if ( m_ObjectType == 1 )  // blob
  {
    A  = vcl_sqrt( Lambda1Abs * Lambda2Abs ) / Lambda3Abs;
    B  = ( Lambda1Abs )/vcl_sqrt( Lambda2Abs * Lambda3Abs );
    S  = vcl_sqrt( Lambda1Sqr + Lambda2Sqr + Lambda3Sqr );

    vesMeasure_1  =
      vcl_exp( -0.5 * ( vnl_math_sqr( A ) ) / AlphaSqr );

    vesMeasure_2  =
      vcl_exp( -0.5 * ( vnl_math_sqr( B ) ) / BetaSqr );

    vesMeasure_3  =
      1 - vcl_exp( -0.5 * vnl_math_sqr( S ) / GammaSqr );

    vesMeasure_4  = 1.0;
  }

  if ( m_ObjectType == 2 )  // vesselness
  {
    A  = vcl_sqrt( Lambda2Abs ) / Lambda3Abs;
    B  = ( Lambda1Abs )/vcl_sqrt( Lambda2Abs * Lambda3Abs );
    S  = vcl_sqrt( Lambda1Sqr + Lambda2Sqr + Lambda3Sqr );

    vesMeasure_1  =
      1 - vcl_exp( -0.5 * ( vnl_math_sqr( A ) ) / AlphaSqr );

    vesMeasure_2  =
      vcl_exp( -0.5 * ( vnl_math_sqr( B ) ) / BetaSqr );

    vesMeasure_3  =
      1 - vcl_exp( -0.5 * vnl_math_sqr( S ) / GammaSqr );

    vesMeasure_4  = 1.0;
  }

  return vesMeasure_1 * vesMeasure_2 * vesMeasure_3;
}


template < typename TPixel >
void
StructMeasureImageFilter< TPixel >::
ComputeLambdaOrdering( double& Lambda1, double& Lambda2, double& Lambda3, EigenValueArrayType& eigenValue )
{
  double smallest, largest;

  // Find the smallest eigenvalue
  smallest = vnl_math_abs( eigenValue[0] );
  Lambda1 = eigenValue[0];

  for ( unsigned int i=1; i <=2; i++ )
    {
    if ( vnl_math_abs( eigenValue[i] ) < smallest )
      {
      Lambda1 = eigenValue[i];
      smallest = vnl_math_abs( eigenValue[i] );
      }
    }

  // Find the largest eigenvalue
  largest = vnl_math_abs( eigenValue[0] );
  Lambda3 = eigenValue[0];

  for ( unsigned int i = 1; i <= 2; i++ )
    {
    if (  vnl_math_abs( eigenValue[i] ) > largest )
      {
      Lambda3 = eigenValue[i];
      largest = vnl_math_abs( eigenValue[i] );
      }
    }

  //  find Lambda2 so that |Lambda1| < |Lambda2| < |Lambda3|
  Lambda2 = eigenValue[0];

  for ( unsigned int i=0; i <=2; i++ )
    {
    if ( ( eigenValue[i] != Lambda1 ) && ( eigenValue[i] != Lambda3 ) )
      {
      Lambda2 = eigenValue[i];
      break;
      }
    }
}


template < typename TPixel >
void
StructMeasureImageFilter< TPixel >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template < typename TPixel >
void
StructMeasureImageFilter< TPixel >
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

template < typename TPixel >
void
StructMeasureImageFilter< TPixel >
::BeforeThreadedGenerateData()
{
  InputRegionType region = this->GetInput()->GetLargestPossibleRegion();

  // Allocate Eigen matrix
  if ( m_ComputeEigenVectors )
  {
    m_EigenMatrixImage = EigenMatrixImageType::New();
    m_EigenMatrixImage->SetRegions( region );
    m_EigenMatrixImage->CopyInformation( this->GetInput() );
    m_EigenMatrixImage->Allocate();
  }

}

template < typename TPixel >
void
StructMeasureImageFilter< TPixel >
::AfterThreadedGenerateData()
{}


template < typename TPixel >
void
StructMeasureImageFilter< TPixel >
::ThreadedGenerateData(const OutputImageRegionType& windowRegion, ThreadIdType threadId)
{
  OutputImagePointer output = this->GetOutput();

  EigenCalculatorType m_EigenCalculator;
  m_EigenCalculator.SetDimension( ImageDimension );

  EigenIteratorType eit;
  if ( m_ComputeEigenVectors )
  {
    eit = EigenIteratorType( m_EigenMatrixImage, windowRegion );
    eit.GoToBegin();
  }

  ImageRegionConstIterator< InputImageType > it( this->GetInput(), windowRegion );
  ImageRegionIterator< OutputImageType > oit( output, windowRegion );

  // walk the region of eigen values and get the vesselness measure
  EigenValueArrayType eigenValue, h;
  EigenMatrixType     eigenMatrix, H;
  double vesselnessMeasure = 0;
  double Lambda1, Lambda2, Lambda3;

  oit.GoToBegin();
  it.GoToBegin();
  while ( !it.IsAtEnd() )
    {
    // Get the eigen value in ascending order
    // Preserves sign of eigenValue but ordering is not useful
    m_EigenCalculator.SetOrderEigenMagnitudes(0);
    m_EigenCalculator.ComputeEigenValuesAndVectors( it.Get(), eigenValue, H );
//    std::cout << eigenValue[0] << ' ' << eigenValue[1] << ' ' << eigenValue[2] << std::endl;

    // Get the eigen value in magnitude ascending order
    // Sign is lost but eigenVectors are useful
    m_EigenCalculator.SetOrderEigenMagnitudes(1);
    m_EigenCalculator.ComputeEigenValuesAndVectors( it.Get(), h, eigenMatrix );
//    std::cout << h[0] << ' ' << h[1] << ' ' << h[2] << std::endl;

    // Recreate magnitude ordering with sign
    ComputeLambdaOrdering( Lambda1, Lambda2, Lambda3, eigenValue );

//    std::cout << Lambda1 << ' ' << Lambda2 << ' ' << Lambda3 << std::endl << std::endl;

    if ( m_ComputeEigenVectors )
    {
      eit.Set( eigenMatrix );
      ++eit;
    }

    vesselnessMeasure = GetStructureMeasure( Lambda1, Lambda2, Lambda3 );

		if(  m_ScaleMembranenessMeasure )
		{
			oit.Set( static_cast< OutputPixelType >( vnl_math_abs( Lambda3 )*vesselnessMeasure ) );
		}
		else
		{
			oit.Set( static_cast< OutputPixelType >( vesselnessMeasure ) );
		}


		if ( m_StructureIsPositive )
		{
			if ( Lambda3 >= - EPSILON )
			{
				oit.Set( NumericTraits< OutputPixelType >::Zero );
			}
		}
		else
		{
			if ( Lambda3 <= EPSILON )
			{
				oit.Set( NumericTraits< OutputPixelType >::Zero );
			}
		}
//     std::cout << vesselnessMeasure << std::endl;
    ++it;
    ++oit;
    }
}


template < typename TPixel >
void
StructMeasureImageFilter< TPixel >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta:  " << m_Beta  << std::endl;
  os << indent << "Gamma: " << m_Gamma << std::endl;
  os << indent << "C: " << m_C << std::endl;
}

} // end namespace itk

#endif
