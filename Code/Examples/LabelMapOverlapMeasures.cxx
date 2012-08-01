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
#include "itkLabelMapOverlapMeasures.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " groundtruth segmentation" << std::endl;
    return EXIT_FAILURE;
    }

  double m_ConfusionThreshold = 0.75;
  double m_ThresholdForUndersegmentation = 0.75;
  if (argc > 3)
  {
    m_ConfusionThreshold = atof( argv[3] );
  }
  if (argc > 4)
  {
    m_ThresholdForUndersegmentation = atof( argv[4] );
  }

  const unsigned int Dimension = 3;
  typedef unsigned int PixelType;
  typedef itk::Image< PixelType, Dimension > SegmentImageType;
  typedef itk::ImageFileReader< SegmentImageType > ReaderType;
  typedef itk::LabelMapOverlapMeasures< SegmentImageType > LabelMapOverlapMeasuresFilterType;

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName ( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName ( argv[2] );
  reader2->Update();

  LabelMapOverlapMeasuresFilterType::Pointer filter = LabelMapOverlapMeasuresFilterType::New();
  filter->SetInput( 0, reader1->GetOutput() );
  filter->SetInput( 1, reader2->GetOutput() );
  filter->SetConfusionThreshold( m_ConfusionThreshold );
  filter->SetThresholdForUndersegmentation( m_ThresholdForUndersegmentation );
  filter->Update();

  std::cout << double(filter->GetNumberOfMatches())/filter->GetNumberOfGroundTruthLabels() << ' ';
  std::cout << double(filter->GetNumberOfMatches())/filter->GetNumberOfSegLabels() << std::endl;

//  std::cout << filter->GetNumberOfOverSegmentations() << ' ';
//  std::cout << filter->GetNumberOfUnderSegmentations() << ' ';
//  std::cout << filter->GetAverageValueOfDiceMetric() << ' ';
//  std::cout << filter->GetAverageValueOfL2Metric() << std::endl;

  return EXIT_SUCCESS;
  }
