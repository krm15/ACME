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
#include "itkBinarydilateImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSizeThresholdImageFilter.h"

#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkRegionOfInterestImageFilter.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " input kernelSize output" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned int, Dimension >  ImageType;
  typedef ImageType::PixelType                   PixelType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  typedef itk::ImageRegionIterator< ImageType > IteratorType;

  typedef itk::BinaryBallStructuringElement < int, Dimension >  StructuringElementType;
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType > DilateFilterType;
  typedef itk::SizeThresholdImageFilter< ImageType, ImageType > SizeFilterType;
  typedef unsigned int                                  LabelType;
  typedef itk::ShapeLabelObject< LabelType, Dimension > ShapeLabelObjectType;
  typedef itk::LabelMap< ShapeLabelObjectType >              ShapeLabelMapType;
  typedef itk::LabelImageToShapeLabelMapFilter< ImageType, ShapeLabelMapType > ShapeConverterType;
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ROIFilterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();
  ImageType::Pointer input = reader->GetOutput();
  input->DisconnectPipeline();
  ImageType::RegionType inputRegion = input->GetLargestPossibleRegion();

  StructuringElementType ball;
  StructuringElementType::SizeType ballSize;

  // ballSize[0] = ballSize[1] = 5;
  // ballSize[2] = atoi( argv[2] );

  ballSize.Fill( atoi( argv[2] ) );

//  ballSize[0] *= 2;
//  ballSize[1] *= 2;
//  ballSize[2] /= 2;

  ball.SetRadius( ballSize );
  ball.CreateStructuringElement();

  ImageType::Pointer outputImg = ImageType::New();
  outputImg->CopyInformation( input );
  outputImg->SetRegions( input->GetLargestPossibleRegion() );
  outputImg->Allocate();
  outputImg->FillBuffer( 0 );

  // convert the image into a collection of objects
  ShapeConverterType::Pointer shapeConverter = ShapeConverterType::New();
  shapeConverter->SetInput ( input );
  shapeConverter->SetBackgroundValue ( 0 );
  shapeConverter->Update();
  ShapeLabelMapType::Pointer shapeLabelMap = shapeConverter->GetOutput();

  unsigned int m_NumberOfLabels = shapeLabelMap->GetNumberOfLabelObjects();

  ImageType::RegionType region;
  ImageType::IndexType  index;
  ImageType::SizeType   size;
  ROIFilterType::Pointer roi = ROIFilterType::New();
  roi->SetInput( input );

  unsigned int j;
  for( unsigned int i = 0; i < m_NumberOfLabels; i++ )
  {
    j = shapeLabelMap->GetLabels()[i];
    ShapeLabelObjectType *labelObject = shapeLabelMap->GetLabelObject ( j );

    region = labelObject->GetBoundingBox();
    region.PadByRadius( ballSize );
    region.Crop( inputRegion );

    roi->SetRegionOfInterest( region );
    roi->Update();

    DilateFilterType::Pointer dilate = DilateFilterType::New();
    dilate->SetInput( roi->GetOutput() );
    dilate->SetKernel( ball );
    dilate->SetForegroundValue( j );
    dilate->Update();

    IteratorType sIt( dilate->GetOutput(), dilate->GetOutput()->GetLargestPossibleRegion() );
    sIt.GoToBegin();
    IteratorType oIt( outputImg, region );
    oIt.GoToBegin();
    while( !sIt.IsAtEnd() )
    {
      if( sIt.Get() == j )
      {
        oIt.Set( j );
      }
      ++sIt;
      ++oIt;
    }
  }

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput ( outputImg );
  writer->SetFileName ( argv[3] );

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
