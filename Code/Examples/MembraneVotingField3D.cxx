/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 585 $  // Revision of last commit
  Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
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
#include "itkSymmetricSecondRankTensor.h"
#include "itkTensorVoting3D.h"
#include "itkTensorToSaliencyImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImage eigenImage outputSaliencyImage sigma";
    std::cerr << "" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned char, Dimension >   FeatureImageType;
  typedef itk::Image< double, Dimension > InputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageRegionIterator< InputImageType > InputIteratorType;

  typedef itk::Matrix< double, Dimension, Dimension> PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > TokenReaderType;
  typedef itk::TensorVoting3D< ImageType > TensorVotingFilterType;
  typedef itk::ImageRegionIterator< ImageType > IteratorType;

  typedef itk::Vector< double, Dimension > VectorType;
  typedef itk::Image< VectorType, Dimension > VectorImageType;
  typedef itk::TensorToSaliencyImageFilter< ImageType, VectorImageType > SaliencyFilterType;
  typedef itk::ImageFileWriter< FeatureImageType > WriterType;
  typedef itk::VectorIndexSelectionCastImageFilter< VectorImageType, InputImageType > IndexFilterType;
  typedef itk::RescaleIntensityImageFilter< InputImageType, FeatureImageType > RescaleFilterType;

  InputImageType::Pointer planarity;
  {
    ReaderType::Pointer sReader = ReaderType::New();
    sReader->SetFileName( argv[1] );
    sReader->Update();
    planarity = sReader->GetOutput();
    planarity->DisconnectPipeline();
  }

  ImageType::Pointer input;
  {
    TokenReaderType::Pointer tReader = TokenReaderType::New();
    tReader->SetFileName( argv[2] );
    tReader->Update();
    input = tReader->GetOutput();
    input->DisconnectPipeline();
  }

  // Rotate around a circle and fill the pixels
  InputIteratorType iIt( planarity, planarity->GetLargestPossibleRegion() );
  IteratorType It( input, input->GetLargestPossibleRegion() );
  iIt.GoToBegin();
  It.GoToBegin();
  PixelType p;
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
  tensorVote->SetInput( input );
  tensorVote->SetSigma( atof( argv[4] ) );
  tensorVote->SetUseSparseVoting( false );
  tensorVote->SetNumberOfThreads( 12 );

  try
    {
    tensorVote->Update();
    }
  catch ( itk::ExceptionObject e )
    {
    std::cerr << "Error: " << e << std::endl;
    }
  std::cout << "Voting complete..." << std::endl;

  SaliencyFilterType::Pointer saliency = SaliencyFilterType::New();
  saliency->SetInput( tensorVote->GetOutput() );
  saliency->Update();

  IndexFilterType::Pointer componentExtractor = IndexFilterType::New();
  componentExtractor->SetInput( saliency->GetOutput() );
  componentExtractor->SetIndex( 2 );
  std::cout << "Saliency complete..." << std::endl;

  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( componentExtractor->GetOutput() );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();

  // Write the output
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput ( rescale->GetOutput() );
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
