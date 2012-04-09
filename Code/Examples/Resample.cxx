/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 667 $  // Revision of last commit
  Date: $Date: 2009-09-16 13:12:21 -0400 (Wed, 16 Sep 2009) $  // Date of last commit
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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkImageRegionIterator.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " input output spacingFactorX spacingFactorY spacingFactorZ" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > InputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageFileWriter< InputImageType > WriterType;
  typedef itk::ResampleImageFilter< InputImageType, InputImageType > ResampleFilterType;
  typedef itk::LinearInterpolateImageFunction< InputImageType > InterpolatorType;
  typedef itk::AffineTransform< double, Dimension > TransformType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();

  float factor;

  InputImageType::PointType origin = reader->GetOutput()->GetOrigin();

  InputImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  InputImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();


  // No z sampling
  for( unsigned int i = 0; i < Dimension; i++ )
    {
    factor = atof(argv[3+i]);
    size[i] = size[i]/factor;
    spacing[i] = spacing[i]*factor;
    }
  std::cout << size << std::endl;

  // create the resample filter, transform and interpolator
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  InterpolatorType::Pointer interp = InterpolatorType::New();

  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform ( transform );
  resample->SetInterpolator ( interp );
  resample->SetInput ( reader->GetOutput() );
  resample->SetSize ( size );
  resample->SetOutputOrigin ( origin );
  resample->SetOutputSpacing ( spacing );
  resample->SetDefaultPixelValue ( 0 );
  resample->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput ( resample->GetOutput() );
  writer->SetFileName ( argv[2] );
  writer->SetUseCompression( true );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject e )
    {
    std::cerr << "Error: " << e << std::endl;
    }

  return EXIT_SUCCESS;
  }
