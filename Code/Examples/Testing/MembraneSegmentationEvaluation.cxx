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
#include "itkBinaryThresholdImageFilter.h"
#include "itkGradientWeightedDistanceImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMorphologicalWatershedImageFilter2.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkShapeRelabelImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkTensorVoting3D.h"
#include "itkTensorToSaliencyImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkMultiscaleStructMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 9 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iMembraneImg iFgImg oSegmentImg <opt> alpha sigma omega threshold" << std::endl;
    return EXIT_FAILURE;
    }

  const int Dimension = 3;
  typedef itk::Image< unsigned char, Dimension >   FeatureImageType;
  typedef itk::Image< double, Dimension >          InputImageType;
  typedef itk::Image< int, Dimension >             SegmentImageType;
  typedef itk::ImageFileReader< FeatureImageType > FeatureReaderType;
  typedef itk::ImageFileReader< InputImageType >   InputReaderType;
  typedef itk::ImageFileReader< SegmentImageType >   SegmentReaderType;

  typedef itk::Matrix< double, Dimension, Dimension> MatrixPixelType;
  typedef itk::Image< MatrixPixelType, Dimension > TensorImageType;
  typedef itk::ImageFileReader< TensorImageType > TokenReaderType;
  typedef itk::TensorVoting3D< TensorImageType > TensorVotingFilterType;
  typedef itk::ImageRegionIterator< TensorImageType > TensorIteratorType;

  typedef itk::Vector< double, Dimension > VectorType;
  typedef itk::Image< VectorType, Dimension > VectorImageType;

  // Declare the type of multiscale vesselness filter
  typedef itk::MultiscaleStructMeasureImageFilter< FeatureImageType, InputImageType> MultiscalePlateFilterType;

  typedef itk::TensorToSaliencyImageFilter< TensorImageType, VectorImageType > SaliencyFilterType;
  typedef itk::VectorIndexSelectionCastImageFilter< VectorImageType, InputImageType > IndexFilterType;
  typedef itk::BinaryThresholdImageFilter< InputImageType, SegmentImageType > ThresholdType;

  typedef itk::GradientWeightedDistanceImageFilter< FeatureImageType, InputImageType, SegmentImageType >
    DistanceFilterType;

  typedef itk::InvertIntensityImageFilter< InputImageType, InputImageType > RInvertType;

  typedef itk::MorphologicalWatershedImageFilter2< InputImageType, SegmentImageType > WatershedFilterType;

  typedef itk::RescaleIntensityImageFilter< InputImageType, FeatureImageType > RescaleFilterType;

  typedef itk::ShapeRelabelImageFilter< SegmentImageType > RelabelFilterType;
  typedef itk::ImageRegionIterator< FeatureImageType > FeatureIteratorType;
  typedef itk::ImageRegionIterator< InputImageType > InputIteratorType;
  typedef itk::ImageRegionIterator< SegmentImageType > SegmentIteratorType;

  int opt = atof( argv[4] );
  double m_Alpha = atof( argv[5] );
  double sigma = atof( argv[6] );
  double omega = atof( argv[7] );
  double m_ThresholdMin = atof( argv[8] );

  // Read in the raw image
  FeatureReaderType::Pointer raw = FeatureReaderType::New();
  raw->SetFileName ( argv[1] );
  raw->Update();
  FeatureImageType::Pointer input = raw->GetOutput();
  input->DisconnectPipeline();
  std::cout << "Read input image..." << std::endl;

  if (opt >= 2)
  {
    // Create a vesselness Filter
    MultiscalePlateFilterType::Pointer MultiscalePlateFilter = MultiscalePlateFilterType::New();
    MultiscalePlateFilter->SetInput( raw->GetOutput() );
    MultiscalePlateFilter->SetObjectType( 0 );
    MultiscalePlateFilter->SetSigmaMin( sigma );
    MultiscalePlateFilter->SetSigmaMax( sigma );
    MultiscalePlateFilter->SetNumberOfSigmaSteps( 1 );
    MultiscalePlateFilter->Update();
    std::cout << "Planarity filtering complete... " << std::endl;

  // Planarity filter
  InputImageType::Pointer planarity = MultiscalePlateFilter->GetOutput();
  TensorImageType::Pointer eigen =  MultiscalePlateFilter->GetEigenMatrix();

    if (opt == 3)
    {
      // Fill the sparse token image
      InputIteratorType iIt( planarity, planarity->GetLargestPossibleRegion() );
      TensorIteratorType It( eigen, eigen->GetLargestPossibleRegion() );
      iIt.GoToBegin();
      It.GoToBegin();
      MatrixPixelType p;
      VectorType u;
      double q;
      while( !It.IsAtEnd() )
      {
        q = iIt.Get();
        u = It.Get()[Dimension-1];

        for( unsigned int i = 0; i < Dimension; i++ )
          for( unsigned int j = 0; j < Dimension; j++ )
          {
            p[i][j] = p[j][i] = u[i]*u[j]*q;
          }
        It.Set( p );
        ++iIt;
        ++It;
      }
      std::cout << "Filled input image..." << std::endl;

      // Do tensor voting
      TensorVotingFilterType::Pointer tensorVote = TensorVotingFilterType::New();
      tensorVote->SetInput( eigen );
      tensorVote->SetSigma( omega );
      tensorVote->SetUseSparseVoting( false );
      tensorVote->SetNumberOfThreads( 12 );
      tensorVote->Update();
      std::cout << "Voting complete..." << std::endl;

      SaliencyFilterType::Pointer saliency = SaliencyFilterType::New();
      saliency->SetInput( tensorVote->GetOutput() );
      saliency->Update();

      IndexFilterType::Pointer componentExtractor = IndexFilterType::New();
      componentExtractor->SetInput( saliency->GetOutput() );
      componentExtractor->SetIndex( 2 );
      componentExtractor->Update();

      RescaleFilterType::Pointer rescale = RescaleFilterType::New();
      rescale->SetInput( componentExtractor->GetOutput() );
      rescale->SetOutputMinimum( 0 );
      rescale->SetOutputMaximum( 255 );
      rescale->Update();
      input = rescale->GetOutput();
      input->DisconnectPipeline();
      std::cout << "Extract planes complete..." << std::endl;
    }
    else
    {
      RescaleFilterType::Pointer rescale = RescaleFilterType::New();
      rescale->SetInput( MultiscalePlateFilter->GetOutput() );
      rescale->SetOutputMinimum( 0 );
      rescale->SetOutputMaximum( 255 );
      rescale->Update();
      input = rescale->GetOutput();
      input->DisconnectPipeline();
    }
  }

  SegmentReaderType::Pointer reader = SegmentReaderType::New();
  reader->SetFileName ( argv[2] );
  reader->Update();
  SegmentImageType::Pointer m_FgImg = reader->GetOutput();
  m_FgImg->DisconnectPipeline();

  SegmentImageType::Pointer m_ForegroundImg = SegmentImageType::New();
  m_ForegroundImg->CopyInformation( input );
  m_ForegroundImg->SetRegions( input->GetLargestPossibleRegion() );
  m_ForegroundImg->Allocate();
  m_ForegroundImg->FillBuffer( 0 );

  SegmentIteratorType fIt ( m_ForegroundImg, m_ForegroundImg->GetLargestPossibleRegion() );
  SegmentIteratorType tIt ( m_FgImg, m_FgImg->GetLargestPossibleRegion() );
  FeatureIteratorType sIt ( input, input->GetLargestPossibleRegion() );

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
  std::cout << "Computed foreground" << std::endl;

  DistanceFilterType::Pointer distFilter = DistanceFilterType::New();
  distFilter->SetInput ( raw->GetOutput() );
  distFilter->SetForeground ( m_ForegroundImg );
  distFilter->SetLargestCellRadius( 8.0 );
  distFilter->SetNucleiSigma ( 0.4 );
  distFilter->SetAlpha( m_Alpha );
  distFilter->Update();
  std::cout << "Computed distance map" << std::endl;

  RInvertType::Pointer idistance = RInvertType::New();
  idistance->SetInput( distFilter->GetOutput() );
  idistance->SetMaximum( 10.0 );
  idistance->Update();
  std::cout << "Inverted distance map" << std::endl;

  WatershedFilterType::Pointer wshed = WatershedFilterType::New();
  wshed->SetInput( idistance->GetOutput() );
  wshed->SetMarkWatershedLine( false );
  wshed->SetLevel( 1.0 );
  wshed->FullyConnectedOn();
  wshed->SetNumberOfThreads( 12 );
  wshed->SetForegroundImage( m_FgImg );
  wshed->Update();
  SegmentImageType::Pointer output = wshed->GetOutput();
  output->DisconnectPipeline();
  std::cout << "Computed watershed segmentation" << std::endl;

  RelabelFilterType::Pointer relabel = RelabelFilterType::New();
  relabel->SetInput ( output );
  relabel->SetBackgroundValue ( 0 );
  relabel->SetReverseOrdering ( false );
  relabel->Update();

  typedef itk::ImageFileWriter< SegmentImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( relabel->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->SetUseCompression( true );
  writer->Update();
  std::cout << "Writing complete" << std::endl;

  return EXIT_SUCCESS;
}
