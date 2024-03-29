/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianSmoothed3DToMembranenessMeasureImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/04/01 21:19:46 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkMultiscaleStructMeasureImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTensorVoting3D.h"
#include "itkTensorToSaliencyImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"

#include "itkTimeProbe.h"

int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iInputImage oSaliencyImage iSigmaFiltering iSigmaVoting <NumberOfThreads> ";
    std::cerr << "<alpha> <beta> <gamma> <saliencyid>" << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const unsigned int    Dimension = 3;
  typedef unsigned short InputPixelType;
  typedef double        InternalPixelType;

  // Declare the types of the images
  typedef itk::Vector< double, Dimension > VectorType;
  typedef itk::Matrix< double, Dimension, Dimension> TokenPixelType;

  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< InternalPixelType, Dimension> InternalImageType;
  typedef itk::Image< TokenPixelType, Dimension > TokenImageType;
  typedef itk::Image< VectorType, Dimension > VectorImageType;

  typedef itk::ImageRegionIterator< TokenImageType > TokenIteratorType;
  typedef itk::ImageRegionIterator< InternalImageType > OutputIteratorType;

  typedef itk::ImageFileReader< InputImageType  > ImageReaderType;
  typedef itk::ImageFileReader< TokenImageType > TokenReaderType;
  typedef itk::ImageFileWriter< InternalImageType > ImageWriterType;

  typedef itk::MultiscaleStructMeasureImageFilter< InputImageType, InternalImageType> MultiscalePlateFilterType;
  typedef itk::TensorVoting3D< TokenImageType > TensorVotingFilterType;
  typedef itk::TensorToSaliencyImageFilter< TokenImageType, VectorImageType > SaliencyFilterType;
  typedef itk::VectorIndexSelectionCastImageFilter< VectorImageType, InternalImageType > IndexFilterType;
  typedef itk::RescaleIntensityImageFilter< InternalImageType, InputImageType > RescaleFilterType;
  typedef itk::ImageFileWriter< InputImageType > WriterType;


  double m_Alpha = 0.5;
  double m_Beta  = 0.5;
  double m_Gamma = 8.0;
  unsigned int saliencyId = 2;

  unsigned int NumberOfThreads = 24;
  if ( argc > 5 )
  {
    NumberOfThreads = atoi( argv[5] );
  }

  if ( argc > 6 )
  {
    m_Alpha = atof( argv[6] );
    m_Beta = atof( argv[7] );
    m_Gamma = atof( argv[8] );
  }

  if ( argc > 9 )
  {
    saliencyId = atoi( argv[9] );
  }


  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();
  std::cout << "Read input image..." << std::endl;

  // Declare the type of multiscale vesselness filter
  MultiscalePlateFilterType::Pointer MultiscalePlateFilter =
    MultiscalePlateFilterType::New();
  MultiscalePlateFilter->SetInput( reader->GetOutput() );
  MultiscalePlateFilter->SetObjectType( 0 );
  MultiscalePlateFilter->SetSigmaMin( atof(argv[3])  );
  MultiscalePlateFilter->SetSigmaMax( atof(argv[3]) );
  MultiscalePlateFilter->SetNumberOfSigmaSteps( 1 );
  MultiscalePlateFilter->SetAlpha( m_Alpha );
  MultiscalePlateFilter->SetBeta( m_Beta );
  MultiscalePlateFilter->SetGamma( m_Gamma );
  try
    {
    MultiscalePlateFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Planarity complete..." << std::endl;

  InternalImageType::Pointer planarity = MultiscalePlateFilter->GetOutput();
  TokenImageType::Pointer eigen = MultiscalePlateFilter->GetEigenMatrix();

  TokenImageType::Pointer token = TokenImageType::New();
  token->SetRegions( eigen->GetLargestPossibleRegion() );
  token->CopyInformation( eigen );
  token->Allocate();

  // Create tokens of all kinds
  OutputIteratorType iIt( planarity, planarity->GetLargestPossibleRegion() );
  TokenIteratorType It( eigen, eigen->GetLargestPossibleRegion() );
  TokenIteratorType tokIt( token, eigen->GetLargestPossibleRegion() );
  iIt.GoToBegin();
  It.GoToBegin();
  tokIt.GoToBegin();
  TokenPixelType p;
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
    tokIt.Set( p );
    ++iIt;
    ++It;
    ++tokIt;
  }
  std::cout << "Filled token image..." << std::endl;



  // Measure time taken
  itk::TimeProbe cputimer;
  cputimer.Start();

  // Do tensor voting
  TensorVotingFilterType::Pointer tensorVote = TensorVotingFilterType::New();
  tensorVote->SetInput( token );
  tensorVote->SetSigma( atof( argv[4] ) );
//  tensorVote->SetSaliencyImage( planarity );
//  tensorVote->SetEigenMatrixImage( eigen );
  tensorVote->SetUseSparseVoting( false );
  tensorVote->SetNumberOfThreads( NumberOfThreads );

  try
    {
    tensorVote->Update();
    }
  catch ( itk::ExceptionObject e )
    {
    std::cerr << "Error: " << e << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Voting complete..." << std::endl;

  cputimer.Stop();
  std::cout << "Tensor voting filter took " << cputimer.GetMean() << " seconds" << std::endl;


  SaliencyFilterType::Pointer saliency = SaliencyFilterType::New();
  saliency->SetInput( tensorVote->GetOutput() );
  saliency->SetComputeEigenMatrix( false );
  saliency->SetNumberOfThreads( 24 );
  saliency->Update();

  IndexFilterType::Pointer componentExtractor = IndexFilterType::New();
  componentExtractor->SetInput( saliency->GetOutput() );
  componentExtractor->SetIndex( saliencyId );
  std::cout << "Saliency complete..." << std::endl;

  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( componentExtractor->GetOutput() );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();

  // Write the output
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput ( rescale->GetOutput() );
  writer->SetFileName ( argv[2] );

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

