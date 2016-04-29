/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1445 $  // Revision of last commit
  Date: $Date: 2010-05-11 11:50:11 -0400 (Tue, 11 May 2010) $  // Date of last commit
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
#include "itkMatrix.h"
#include "itkMembraneScaleSelectionImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkDiscreteHessianGaussianImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkStructMeasureImageFilter.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " InputMembraneImage EigenVectorImage PlanarityFunctionImage" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef unsigned char                             FeaturePixelType;
  typedef itk::Image< FeaturePixelType, Dimension > FeatureImageType;

  typedef double                                    InputPixelType;
  typedef itk::Image< double, Dimension >           InputImageType;

  typedef itk::HessianRecursiveGaussianImageFilter< FeatureImageType > HessianFilterType;

  // Define image of matrix pixel type
  typedef itk::SymmetricSecondRankTensor< double, Dimension > TensorType;
  typedef itk::Image< TensorType, Dimension>  TensorImageType;

  typedef itk::MembraneScaleSelectionImageFilter< TensorImageType, InputImageType > ScaleSelectionFilterType;
  typedef itk::DiscreteHessianGaussianImageFunction< FeatureImageType, InputPixelType > HessianGaussianImageFunctionType;

  typedef itk::StructMeasureImageFilter< double > StructFilterType;
  typedef StructFilterType::EigenMatrixImageType EigenMatrixImageType;

  typedef itk::ImageFileReader< FeatureImageType  > ImageReaderType;
  typedef itk::ImageFileWriter< EigenMatrixImageType >   MatrixWriterType;
  typedef itk::ImageFileWriter< InputImageType >    InputWriterType;
  typedef itk::ImageRegionIteratorWithIndex< TensorImageType > TensorIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > InputIteratorType;

  // Read membrane image
  ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();
  std::cout << "Reading input image : " << argv[1] << std::endl;

  // Use a blind sigma estimate of 0.45 to compute the Hessian
  HessianFilterType::Pointer blindHessianFilter = HessianFilterType::New();
  blindHessianFilter->SetInput( reader->GetOutput() );
  blindHessianFilter->SetSigma( 1.0 );
  blindHessianFilter->SetNormalizeAcrossScale( true );
  blindHessianFilter->Update();
  TensorImageType::Pointer hessian = blindHessianFilter->GetOutput();
  std::cout << "Blind scale Hessian computed..." << std::endl;

  // Using the Hessian eigen system, compute a new sigma at each pixel
  ScaleSelectionFilterType::Pointer scaleFilter = ScaleSelectionFilterType::New();
  scaleFilter->SetInput ( hessian );
  scaleFilter->SetSigmaMin( 0.4 );
  scaleFilter->SetSigmaMax( 1.0 );
  scaleFilter->SetNumberOfThreads(12);
  scaleFilter->Update();
  InputImageType::Pointer sigmaImage = scaleFilter->GetOutput();
  std::cout << "Sigma scale computed..." << std::endl;

  HessianGaussianImageFunctionType::Pointer hessianFunction = HessianGaussianImageFunctionType::New();
  hessianFunction->SetInputImage( reader->GetOutput() );
  hessianFunction->SetNormalizeAcrossScale( true );
  hessianFunction->SetUseImageSpacing( true );
  hessianFunction->SetInterpolationMode( HessianGaussianImageFunctionType::NearestNeighbourInterpolation );
  hessianFunction->Initialize( );

  TensorIteratorType hIt( hessian, hessian->GetLargestPossibleRegion() );
  InputIteratorType sIt( sigmaImage, sigmaImage->GetLargestPossibleRegion() );
  hIt.GoToBegin();
  sIt.GoToBegin();
  {
    hessianFunction->SetSigma( sIt.Get() );//0.7
    hIt.Set( hessianFunction->EvaluateAtIndex( sIt.GetIndex() ) );
    ++hIt;
    ++sIt;
  }
  std::cout << "Scale-dependent Hessian computed..." << std::endl;

  StructFilterType::Pointer m_StructFilter = StructFilterType::New();
  m_StructFilter->SetInput ( hessian );
  m_StructFilter->SetObjectType( 0 );
  m_StructFilter->Update();
  std::cout << "Planarity function computed..." << std::endl;

  // Write scaled hessian image
  MatrixWriterType::Pointer hessianWriter = MatrixWriterType::New();
  hessianWriter->SetFileName( argv[2] );
  hessianWriter->SetInput ( m_StructFilter->GetEigenMatrix() );

  // Write planarity output
  InputWriterType::Pointer planarityFunction = InputWriterType::New();
  planarityFunction->SetFileName( argv[3] );
  planarityFunction->SetInput (  m_StructFilter->GetOutput() );

  try
    {
    hessianWriter->Update();
    planarityFunction->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }
