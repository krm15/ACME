/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1367 $  // Revision of last commit
  Date: $Date: 2010-04-30 03:01:43 -0400 (Fri, 30 Apr 2010) $  // Date of last commit
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

#ifndef __itkLabelMapOverlapMeasures_txx
#define __itkLabelMapOverlapMeasures_txx

#include "itkLabelMapOverlapMeasures.h"

namespace itk
{
template < class TInputImage >
LabelMapOverlapMeasures< TInputImage >::LabelMapOverlapMeasures()
{
  m_NumberOfMatches = 0;
  m_NumberOfOverSegmentations = 0;
  m_NumberOfUnderSegmentations = 0;
  m_AverageValueOfDiceMetric = 0.0;
  m_AverageValueOfL2Metric = 0.0;

  m_ConfusionThreshold = 0.75;
  m_ThresholdForUndersegmentation = 0.75;

  this->Superclass::SetNumberOfRequiredInputs( 2 );
}


template< class TInputImage >
void
LabelMapOverlapMeasures< TInputImage >::
IdentifyNumberOfCells()
{
  ShapeConverterPointer gShapeConverter = ShapeConverterType::New();
  gShapeConverter->SetInput( this->GetInput( 0 ) );
  gShapeConverter->SetBackgroundValue( 0 );
  gShapeConverter->Update();
  ShapeLabelMapPointer gShapeLabelMap = gShapeConverter->GetOutput();
  m_NumberOfGroundTruthLabels = gShapeLabelMap->GetNumberOfLabelObjects();
  std::cout << "Ng = " << m_NumberOfGroundTruthLabels << std::endl;

  ShapeConverterPointer sShapeConverter = ShapeConverterType::New();
  sShapeConverter->SetInput( this->GetInput( 1 ) );
  sShapeConverter->SetBackgroundValue( 0 );
  sShapeConverter->Update();
  ShapeLabelMapPointer sShapeLabelMap = sShapeConverter->GetOutput();
  m_NumberOfSegLabels = sShapeLabelMap->GetNumberOfLabelObjects();
  std::cout << "Ns = " << m_NumberOfSegLabels << std::endl;
}

template< class TInputImage >
double
LabelMapOverlapMeasures< TInputImage >::
EvaluateHaussdorf( unsigned int iR, unsigned int iC )
{
  BinaryThresholdFilterPointer thresh1 = BinaryThresholdFilterType::New();
  thresh1->SetLowerThreshold ( iR );
  thresh1->SetUpperThreshold ( iR );
  thresh1->SetInsideValue ( 1 );
  thresh1->SetOutsideValue ( 0 );
  thresh1->SetInput ( this->GetInput( 0 ) );
  thresh1->Update();

  BinaryThresholdFilterPointer thresh2 = BinaryThresholdFilterType::New();
  thresh2->SetLowerThreshold ( iC );
  thresh2->SetUpperThreshold ( iC );
  thresh2->SetInsideValue ( 1 );
  thresh2->SetOutsideValue ( 0 );
  thresh2->SetInput ( this->GetInput( 1 ) );
  thresh2->Update();

  HausdorffDistanceFilterPointer hdorff = HausdorffDistanceFilterType::New();
  hdorff->SetInput1( thresh1->GetOutput() );
  hdorff->SetInput2( thresh2->GetOutput() );
  hdorff->SetUseImageSpacing( true );
  hdorff->Update();
  double dist = hdorff->GetAverageHausdorffDistance();
  return dist;
}


template< class TInputImage >
void
LabelMapOverlapMeasures< TInputImage >::
CreateConfusionMatrix()
{
  m_OverlapArea.set_size( m_NumberOfGroundTruthLabels+1, m_NumberOfSegLabels+1 );
  m_OverlapArea.fill( 0 );

  m_GroundTruthArea.set_size( m_NumberOfGroundTruthLabels+1 );
  m_GroundTruthArea.fill( 0 );

  m_SegmentedArea.set_size( m_NumberOfSegLabels+1 );
  m_SegmentedArea.fill( 0 );

  unsigned int r, c;

  ImageConstIteratorType gIt( this->GetInput(0), this->GetInput(0)->GetLargestPossibleRegion() );
  ImageConstIteratorType sIt( this->GetInput(1), this->GetInput(1)->GetLargestPossibleRegion() );
  gIt.GoToBegin();
  sIt.GoToBegin();

  while( !sIt.IsAtEnd() )
  {
    r = static_cast<unsigned int>( gIt.Get() );
    m_GroundTruthArea[r]++;

    c = static_cast<unsigned int>( sIt.Get() );
    m_SegmentedArea[c]++;

    m_OverlapArea[r][c]++;
    ++gIt;
    ++sIt;
  }

  double s;
  double den = 0;
  PairVectorType mv;
  for( unsigned int i = 1; i <= m_NumberOfGroundTruthLabels; i++ )
    for( unsigned int j = 1; j <= m_NumberOfSegLabels; j++ )
    {
      den = static_cast<double>( m_GroundTruthArea[i] + m_SegmentedArea[j] - m_OverlapArea[i][j] );
      s = static_cast<double>(m_OverlapArea[i][j])/den;
      mv[0] = i;
      mv[1] = j;
      if ( s > m_ConfusionThreshold )
      {
        while ( m_Dice.find(s) != m_Dice.end() )
          s = s-0.00001;
        m_Dice[s] = mv;
      }
    }
}

template< class TInputImage >
void
LabelMapOverlapMeasures< TInputImage >::
AssignMatchingLabels()
{
  // Finding matching and store in m_TrackLabels
  // Some divided cells are also matched up on the parent id
  if ( m_Dice.size() > 0 )
  {
    PairVectorType mv;
    DiceMapType::iterator loc;
    unsigned int r, c;
    loc = m_Dice.end();
    do
    {
      --loc;
      mv = loc->second;
      r = mv[0];
      c = mv[1];

      if ( ( m_GroundTruthLabels.find(r) == m_GroundTruthLabels.end() ) &&
        ( m_SegLabels.find(c) == m_SegLabels.end() ) )
      {
        //m_AverageValueOfDiceMetric += loc->first;
        //m_AverageValueOfL2Metric += EvaluateHaussdorf(r,c);
        m_GroundTruthLabels[r] = c;
        m_SegLabels[c] = r;
      }
    }
    while(loc != m_Dice.begin() );
  }
  m_NumberOfMatches = m_SegLabels.size();
  //m_AverageValueOfDiceMetric /= m_NumberOfMatches;
  //m_AverageValueOfL2Metric /= m_NumberOfMatches;
}


template< class TInputImage >
void
LabelMapOverlapMeasures< TInputImage >::
IdentifyMismatch()
{
  double groundTruthArea, segArea;
  double overlapArea;
  for( unsigned int i = 1; i <= m_NumberOfGroundTruthLabels; i++ )
  {
    // Check that the label is not matched
    if ( m_GroundTruthLabels.find(i) == m_GroundTruthLabels.end() )
    {
      groundTruthArea = m_GroundTruthArea[i];
      for( unsigned int j = 1; j <= m_NumberOfSegLabels; j++ )
      {
        segArea = m_SegmentedArea[j];
        overlapArea = m_OverlapArea[i][j];
      }

      if ( overlapArea < m_ThresholdForUndersegmentation * groundTruthArea )
      {
        // oversegmented
        m_NumberOfOverSegmentations++;
      }
      else
      {
        // undersegmented
        m_NumberOfUnderSegmentations++;
      }
    }
  }
}

template < class TInputImage >
void
LabelMapOverlapMeasures< TInputImage >::
GenerateData()
{
  // Number of cells
  IdentifyNumberOfCells();
  std::cout << "Identified number of cells..." << std::endl;

  // Create a map of confusion matrix entries above the threshold
  CreateConfusionMatrix();
  std::cout << "Creating confusion matrix complete..." << std::endl;

  // 1:1 mapping assignment
  AssignMatchingLabels();
  std::cout << "Matching complete..." << std::endl;

  // Identify over and undersegmentations
  IdentifyMismatch();
  std::cout << "Mismatch id complete..." << std::endl;
}


template < class TInputImage >
void
LabelMapOverlapMeasures< TInputImage >::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

} /* end namespace itk */

#endif
