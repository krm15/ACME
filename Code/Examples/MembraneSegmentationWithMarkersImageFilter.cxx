/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 817 $  // Revision of last commit
  Date: $Date: 2009-11-07 17:21:12 -0500 (Sat, 07 Nov 2009) $  // Date of last commit
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
#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkGradientWeightedDistanceImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMorphologicalWatershedFromMarkersImageFilter2.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iMembraneImg iTensorSaliencyImg iSeedImage";
    std::cerr << " oSegmentImg threshold <fgImage> " << std::endl;
    return EXIT_FAILURE;
    }

  const int Dimension = 3;
  typedef itk::Image< unsigned char, Dimension >   FeatureImageType;
  typedef itk::Image< double, Dimension >          InputImageType;
  typedef itk::Image< int, Dimension >             SegmentImageType;
  typedef itk::ImageFileReader< FeatureImageType > FeatureReaderType;
  typedef itk::ImageFileReader< InputImageType >   InputReaderType;
  typedef itk::ImageFileReader< SegmentImageType >   SegmentReaderType;

  typedef itk::GradientWeightedDistanceImageFilter< FeatureImageType, InputImageType, SegmentImageType >
    DistanceFilterType;

  typedef itk::InvertIntensityImageFilter< InputImageType, InputImageType > RInvertType;

  typedef itk::MorphologicalWatershedFromMarkersImageFilter2< InputImageType, SegmentImageType > WatershedFilterType;

  typedef itk::ImageRegionIterator< SegmentImageType > SegmentIteratorType;

  typedef itk::RescaleIntensityImageFilter< InputImageType, FeatureImageType > RescaleFilterType;
  typedef itk::ImageRegionIterator< FeatureImageType > FeatureIteratorType;
  typedef itk::ImageRegionIterator< SegmentImageType > SegmentIteratorType;

  double m_Alpha = 1.0;
  double m_Beta = 0.0;
  double m_ThresholdMin = atof( argv[5] );

  // Read in the raw image
  FeatureReaderType::Pointer raw = FeatureReaderType::New();
  raw->SetFileName ( argv[1] );
  raw->Update();

  // Tensor saliency image
  FeatureReaderType::Pointer sreader = FeatureReaderType::New();
  sreader->SetFileName ( argv[2] );
  sreader->Update();

  SegmentImageType::Pointer m_FgImg = SegmentImageType::New();
  m_FgImg->CopyInformation( sreader->GetOutput() );
  m_FgImg->SetRegions( sreader->GetOutput()->GetLargestPossibleRegion() );
  m_FgImg->Allocate();
  m_FgImg->FillBuffer( 1 );

  if ( argc > 6 )
  {
    std::fstream inFile( argv[6],std::ios::in );
    if ( inFile.is_open() )
    {
      inFile.close();
      SegmentReaderType::Pointer reader = SegmentReaderType::New();
      reader->SetFileName ( argv[6] );
      reader->Update();
      m_FgImg = reader->GetOutput();
      m_FgImg->DisconnectPipeline();
    }
  }

  SegmentImageType::Pointer m_ForegroundImg = SegmentImageType::New();
  m_ForegroundImg->CopyInformation( sreader->GetOutput() );
  m_ForegroundImg->SetRegions( sreader->GetOutput()->GetLargestPossibleRegion() );
  m_ForegroundImg->Allocate();
  m_ForegroundImg->FillBuffer( 0 );

  SegmentIteratorType fIt ( m_ForegroundImg, m_ForegroundImg->GetLargestPossibleRegion() );
  SegmentIteratorType tIt ( m_FgImg, m_FgImg->GetLargestPossibleRegion() );
  FeatureIteratorType sIt ( sreader->GetOutput(), sreader->GetOutput()->GetLargestPossibleRegion() );

  fIt.GoToBegin();
  tIt.GoToBegin();
  sIt.GoToBegin();
  while ( !tIt.IsAtEnd() )
    {
    if ( tIt.Get() > 0 )
      {
      if ( sIt.Get() > m_ThresholdMin )
        {
        fIt.Set ( 1 );
        }
      else
        {
        fIt.Set ( 0 );
        }
      }
      else
      {
      fIt.Set( 1 );
      }
    ++tIt;
    ++sIt;
    ++fIt;
  }

  DistanceFilterType::Pointer distFilter = DistanceFilterType::New();
  distFilter->SetInput ( raw->GetOutput() );
  distFilter->SetUseLevelSet( false );
  distFilter->SetLargestCellRadius( 4.0 );
  distFilter->SetForeground ( m_ForegroundImg );
  distFilter->SetNucleiSigma ( 0.4 );
  distFilter->SetAlpha( m_Alpha );
  distFilter->SetBeta( m_Beta );
  distFilter->Update();
  std::cout << "Computed distance map" << std::endl;

  RInvertType::Pointer idistance = RInvertType::New();
  idistance->SetInput( distFilter->GetOutput() );
  idistance->SetMaximum( 10.0 );
  idistance->Update();
  std::cout << "Inverted distance map" << std::endl;

  SegmentImageType::Pointer m_MarkerImg;
  {
    SegmentReaderType::Pointer mReader = SegmentReaderType::New();
    mReader->SetFileName ( argv[3] );
    mReader->Update();
    m_MarkerImg = mReader->GetOutput();
    m_MarkerImg->DisconnectPipeline();
  }

  WatershedFilterType::Pointer wshed = WatershedFilterType::New();
  wshed->SetInput( idistance->GetOutput() );
  wshed->SetMarkWatershedLine( false );
  wshed->FullyConnectedOn();
  wshed->SetNumberOfThreads( 12 );
  wshed->SetMarkerImage( m_MarkerImg );
  wshed->SetForegroundImage( m_FgImg );
  wshed->Update();
  SegmentImageType::Pointer output = wshed->GetOutput();
  output->DisconnectPipeline();
  std::cout << "Computed watershed segmentation" << std::endl;

  typedef itk::ImageFileWriter< SegmentImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );//relabel->GetOutput()
  writer->SetFileName( argv[4] );
  writer->SetUseCompression( true );
  writer->Update();

  return EXIT_SUCCESS;
}
