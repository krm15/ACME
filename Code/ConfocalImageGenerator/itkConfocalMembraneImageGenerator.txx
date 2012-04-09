/*=========================================================================
  Author: $Author: agouaillard $  // Author of last commit
  Version: $Rev: 868 $  // Revision of last commit
  Date: $Date: 2009-11-20 18:42:47 -0500 (Fri, 20 Nov 2009) $  // Date of last commit
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

#ifndef __itkConfocalMembraneImageGenerator_txx
#define __itkConfocalMembraneImageGenerator_txx

#include "itkConfocalMembraneImageGenerator.h"

namespace itk
{
template < class TInputImage >
ConfocalMembraneImageGenerator< TInputImage >
::ConfocalMembraneImageGenerator()
{
  m_NumberOfClasses = 5;
  m_Membrane = 0;
}


template < class TInputImage >
void
ConfocalMembraneImageGenerator< TInputImage >::
ComputeNoiseModel()
{
  SPFilterPointer sp = SPFilterType::New();
  sp->SetInput( m_Membrane );
  sp->SetProbability( 0.07 );
  sp->Update();

  // Implement an exact convolution with PSF
  PSFFilterPointer psf2 = PSFFilterType::New();
  psf2->SetInput( sp->GetOutput() );
  psf2->SetVariance( 0.5 );
  psf2->SetMaximumKernelWidth( 3 );
  psf2->Update();

  PoissonFilterPointer poisson = PoissonFilterType::New();
  poisson->SetInput( psf2->GetOutput() );
  poisson->SetScale( 1.0 );
  poisson->Update();

//   SpeckleFilterPointer speckle = SpeckleFilterType::New();
//   speckle->SetInput( poisson->GetOutput() );
//   speckle->SetStandardDeviation( 0.3 );
//   speckle->Update();

  NormalFilterPointer gauss = NormalFilterType::New();
  gauss->SetInput( poisson->GetOutput() );
  gauss->SetStandardDeviation( 7.0 );
  gauss->SetMean( 0 );
  gauss->Update();

  m_Membrane = gauss->GetOutput();
  m_Membrane->DisconnectPipeline();
}


template < class TInputImage >
void
ConfocalMembraneImageGenerator< TInputImage >::
GenerateVoronoiTesselation()
{
  if ( !m_Membrane )
  {
    ImageRegionType region;
    ImageIndexType index;
    index[0] = index[1] = index[2] = 0;

    region.SetIndex( index );
    region.SetSize( m_Size );

    SamplePointer sample = SampleType::New();
    MeasurementVectorType mv;

    m_Membrane = ImageType::New();
    m_Membrane->SetRegions( region );
    m_Membrane->SetSpacing( m_Spacing );
    m_Membrane->SetOrigin( m_Origin );
    m_Membrane->Allocate();

    ImagePointType p;
    IteratorType It( m_Membrane, m_Membrane->GetLargestPossibleRegion() );
    It.GoToBegin();
    while( !It.IsAtEnd() )
    {
      index = It.GetIndex();
      m_Membrane->TransformIndexToPhysicalPoint( index, p );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        mv[i] = static_cast< float > ( p[i] );
      sample->PushBack(mv);
      ++It;
    }

    TreeGeneratorPointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample(sample);
    treeGenerator->SetBucketSize(16);
    treeGenerator->Update();

    EstimatorPointer estimator = EstimatorType::New();

    vnl_random r;
    typename EstimatorType::ParametersType initialMeans( m_NumberOfClasses * ImageDimension );
    for( unsigned int i = 0; i < m_NumberOfClasses; i++ )
      for( unsigned int j = 0; j < ImageDimension; j++ )
        initialMeans[ i*ImageDimension+j ] = r.drand32(0, 1) * m_Size[j]* m_Spacing[j];

    estimator->SetParameters(initialMeans);
    estimator->SetKdTree(treeGenerator->GetOutput());
    estimator->SetMaximumIteration( 20 );
    estimator->SetCentroidPositionChangesThreshold( 2.5 );
    estimator->StartOptimization();
    typename EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();

    DecisionRulePointer decisionRule = DecisionRuleType::New();
    ClassifierPointer classifier = ClassifierType::New();
    classifier->SetDecisionRule( (itk::DecisionRuleBase::Pointer) decisionRule);
    classifier->SetSample( sample );
    classifier->SetNumberOfClasses( m_NumberOfClasses );

    std::vector< unsigned int > classLabels;
    classLabels.resize( m_NumberOfClasses );
    for( unsigned int i = 0; i < m_NumberOfClasses; i++ )
      classLabels[i] = i;
    classifier->SetMembershipFunctionClassLabels( classLabels );

    std::vector< MembershipFunctionPointer > membershipFunctions;
    typename MembershipFunctionType::OriginType origin( sample->GetMeasurementVectorSize() );
    int t = 0;
    for( unsigned int i = 0; i < m_NumberOfClasses ; i++ )
      {
      membershipFunctions.push_back( MembershipFunctionType::New() );
      for( unsigned int j = 0;j < sample->GetMeasurementVectorSize(); j++ )
        origin[j] = estimatedMeans[t++];
      membershipFunctions[i]->SetOrigin( origin );
      classifier->AddMembershipFunction( membershipFunctions[i].GetPointer() );
      }
    classifier->Update();

    typename ClassifierType::OutputType* membershipSample = classifier->GetOutput();
    typename ClassifierType::OutputType::ConstIterator iter = membershipSample->Begin();

    It.GoToBegin();
    while( iter != membershipSample->End() )
    {
      It.Set( static_cast<ImagePixelType>( iter.GetClassLabel() ) );
      ++iter;
      ++It;
    }

    PSFFilterPointer psf1 = PSFFilterType::New();
    psf1->SetInput( m_Membrane );
    psf1->SetVariance( 0.5 );
    psf1->SetMaximumKernelWidth( 3 );
    psf1->Update();

    ImagePixelType m,n;
    IteratorType sIt( psf1->GetOutput(), m_Membrane->GetLargestPossibleRegion() );
    sIt.GoToBegin();
    It.GoToBegin();
    while( !It.IsAtEnd() )
    {
      n = It.Get();
      m = sIt.Get();
      if ( m == static_cast<float>( n ) )
        It.Set( 0 );
      else
        {
        if ( r.drand32(0, 1) > 0.85 )
          It.Set( 100 );
        }
      ++It;
      ++sIt;
    }
  }
}

template < class TInputImage >
void
ConfocalMembraneImageGenerator< TInputImage >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name:        " << GetNameOfClass() << std::endl;
  os << indent << "Number of classes: " << this->m_NumberOfClasses << std::endl;
}

} /* end namespace itk */

#endif
