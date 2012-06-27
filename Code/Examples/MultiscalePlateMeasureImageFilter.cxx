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

int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " InputImage PlanarityOutputImage eigenMatrixOutput sigma" << std::endl;
    return EXIT_FAILURE;
    }


  // Define the dimension of the images
  const unsigned int    Dimension = 3;
  typedef unsigned char InputPixelType;
  typedef double        OutputPixelType;

  // Declare the types of the images
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< OutputPixelType, Dimension> OutputImageType;
  typedef itk::ImageFileReader< InputImageType  > ImageReaderType;
  typedef itk::ImageFileWriter< OutputImageType > ImageWriterType;

  ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName ( argv[1] );

  std::cout << "Reading input image : " << argv[1] << std::endl;
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err )
    {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Declare the type of multiscale vesselness filter
  typedef itk::MultiscaleStructMeasureImageFilter<
    InputImageType, OutputImageType> MultiscalePlateFilterType;

  // Create a vesselness Filter
  MultiscalePlateFilterType::Pointer MultiscalePlateFilter =
    MultiscalePlateFilterType::New();
  MultiscalePlateFilter->SetInput( reader->GetOutput() );
  MultiscalePlateFilter->SetObjectType( 0 );
  MultiscalePlateFilter->SetSigmaMin( atof(argv[4])  );
  MultiscalePlateFilter->SetSigmaMax( atof(argv[4]) );
  MultiscalePlateFilter->SetNumberOfSigmaSteps( 1 );


  try
    {
    MultiscalePlateFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Writing out the enhanced image to " <<  argv[2]
    << std::endl;

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( MultiscalePlateFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Matrix< double, Dimension, Dimension> MatrixType;
  typedef itk::Image< MatrixType, Dimension > EigenMatrixImageType;
  typedef itk::ImageFileWriter< EigenMatrixImageType > EigenMatrixWriterType;
  EigenMatrixWriterType::Pointer eigMatrixWriter = EigenMatrixWriterType::New();
  eigMatrixWriter->SetFileName( argv[3] );
  eigMatrixWriter->SetInput ( MultiscalePlateFilter->GetEigenMatrix() );

  try
    {
    eigMatrixWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

