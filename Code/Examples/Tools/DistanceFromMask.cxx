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
#include "itkImageRegionIterator.h"
#include "itkDanielssonDistanceMapImageFilter.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputMask seeds outputFile" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned int, Dimension > ImageType;
  typedef itk::Image< float, Dimension > InputImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::DanielssonDistanceMapImageFilter< ImageType, InputImageType, ImageType > DistanceFilterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();
  ImageType::Pointer mask = reader->GetOutput();

  // Danielsson distance filter
  DistanceFilterType::Pointer dist = DistanceFilterType::New();
  dist->SetInput( mask );
  dist->UseImageSpacingOn();
  dist->UpdateLargestPossibleRegion();
  InputImageType::Pointer distanceMap = dist->GetOutput();

  // Open the seed file
  std::fstream infile;
  infile.open ( argv[2],std::ios::in );
  unsigned int numOfLines;
  infile >> numOfLines;

  std::fstream outFile;
  outFile.open ( argv[3],std::ios::out );
  outFile << numOfLines << std::endl;

  ImageType::PointType pt;
  ImageType::IndexType idx;
  unsigned int cellID;

  for ( unsigned int j = 0; j < numOfLines; j++ )
  {
    for ( unsigned int i = 0; i < Dimension; i++ )
    {
      infile >> pt[i];
    }
    infile >> cellID;

    mask->TransformPhysicalPointToIndex( pt, idx );
    outFile << cellID << ' ' << distanceMap->GetPixel( idx ) << std::endl;
  }
  infile.close();
  outFile.close();

  return EXIT_SUCCESS;
  }
