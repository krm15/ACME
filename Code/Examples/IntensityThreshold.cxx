/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1471 $  // Revision of last commit
  Date: $Date: 2010-05-18 12:53:51 -0400 (Tue, 18 May 2010) $  // Date of last commit
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

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " input output threshold" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;

  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< float, Dimension > InputImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::BinarySizeThresholdImageFilter< ImageType, ImageType > SizeThresholdFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType  > BinaryThresholdFilterType;

  ImageType::Pointer input;
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
    reader->Update();
    input = reader->GetOutput();
    input->DisconnectPipeline();
  }

  int t = atoi( argv[3] );

  // Binarize
  BinaryThresholdFilterType::Pointer thresh = BinaryThresholdFilterType::New();
  thresh->SetUpperThreshold( t );
  thresh->SetLowerThreshold( 0 );
  thresh->SetInsideValue( 0 );
  thresh->SetOutsideValue( 255 );
  thresh->SetInput( input );
  thresh->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( thresh->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

    return EXIT_SUCCESS;
}
