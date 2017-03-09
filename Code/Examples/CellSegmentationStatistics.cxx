/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 870 $  // Revision of last commit
  Date: $Date: 2009-11-23 13:44:26 -0500 (Mon, 23 Nov 2009) $  // Date of last commit
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
#include "itkCellSegmentationStatistics.h"
#include "itkCastImageFilter.h"
#include <map>
#include <fstream>
#include <stdio.h>

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " segmentImage rawImage statFile"
    << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned char, Dimension >  FeatureImageType;
  typedef itk::Image< float, Dimension >          InputImageType;
  typedef itk::Image< unsigned int, Dimension >   SegmentImageType;

  typedef itk::ImageFileReader< FeatureImageType > FeatureReaderType;
  typedef itk::ImageFileReader< SegmentImageType > SegmentReaderType;
  typedef itk::ImageFileWriter< SegmentImageType > SegmentWriterType;

  typedef itk::CellSegmentationStatistics< FeatureImageType, InputImageType, SegmentImageType > FilterType;

  SegmentImageType::Pointer voronoi;
  {
    SegmentReaderType::Pointer vreader = SegmentReaderType::New();
    vreader->SetFileName ( argv[1] );
    vreader->Update();
    voronoi = vreader->GetOutput();
    voronoi->DisconnectPipeline();
  }

  FeatureImageType::Pointer rawImg;
  {
    FeatureReaderType::Pointer reader_raw = FeatureReaderType::New();
    reader_raw->SetFileName ( argv[2] );
    reader_raw->Update();
    rawImg = reader_raw->GetOutput();
    rawImg->DisconnectPipeline();
  }

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput ( voronoi );
  filter->SetRawImage ( rawImg );
  filter->SetLargestCellRadius ( 5.0 );
  filter->GenerateData();

#ifndef NDEBUG
  std::cout << "Writing out cell statistics...." << std::endl;
#endif


  std::fstream statFile;
  statFile.open ( argv[3], std::ios::out | std::ios::app);
  if ( statFile.is_open() )
  {
    unsigned int labels = filter->GetNumOfLabels();
    //statFile << labels << std::endl;
    for ( unsigned int i = 0; i < labels; i++ )
    {
// std::cout << i << ' ' << filter->m_LabelLookup[i] << std::endl;
        statFile << filter->m_EllipsoidDiameter[i] << std::endl;
        //statFile << filter->m_Moments[i] << std::endl;

      // statFile << filter->m_Size[i] << std::endl;
      //statFile << filter->m_PhysicalSize[i] << ' ' << filter->m_Mean[i] << ' ' << filter->m_Mean[i]*filter->m_Size[i] << std::endl;

      //statFile << filter->m_OrientedBoundingBoxSize[i][0] << ' '
      //    << filter->m_OrientedBoundingBoxSize[i][1] << ' ' << filter->m_OrientedBoundingBoxSize[i][2] << std::endl;

      //std::cout << filter->m_Size[i] << ' ' << std::endl;
//        statFile << i << ' ' << std::endl;
//      statFile << filter->m_LabelLookup[i] << ' ';
//       statFile << "PhysicalSize:                 " << filter->m_PhysicalSize[i] << std::endl;
//        statFile << "Mean:                         " << filter->m_Mean[i] << std::endl;
//        statFile << "Size in voxels:               " << filter->m_Size[i] << std::endl;
//        statFile << "Centroid:                     " << filter->m_Centroid[i] << std::endl;
//        statFile << "CenterOfGravity:              " << filter->m_CenterOfGravity[i] << std::endl;
//        statFile << "Moments about principal axis: " << filter->m_Moments[i] << std::endl;
//        statFile << "Binary moments about axis     " << filter->m_BinaryMoments[i] << std::endl;
//        statFile << "Binary Principal axis:        " << std::endl << filter->m_BinaryPrincipalAxis[i] << std::endl;
    }
  }
  statFile.close();

  return EXIT_SUCCESS;
  }
