/*=========================================================================
  Author: $Author: arnaudgelas $  // Author of last commit
  Version: $Rev: 567 $  // Revision of last commit
  Date: $Date: 2009-08-17 11:47:32 -0400 (Mon, 17 Aug 2009) $  // Date of last commit
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

#ifndef __itkSizeThresholdImageFilter_txx
#define __itkSizeThresholdImageFilter_txx

#include "itkSizeThresholdImageFilter.h"

namespace itk
{
//  Software Guide : BeginCodeSnippet
template < class TInput, class TOutput >
SizeThresholdImageFilter< TInput, TOutput >::
SizeThresholdImageFilter()
{
  m_InputCast = InputCastType::New();
  m_OutputCast = OutputCastType::New();
  m_Map2Image = MapToImageType::New();

  m_MinComponentSize = 200;
  m_MaxComponentSize = 20000;
  m_MaxComponent = false;
  m_Label = 1;

  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->Superclass::SetNumberOfRequiredOutputs( 1 );

  this->Superclass::SetNthOutput( 0, TOutput::New() );
}


template < class TInput, class TOutput >
void
SizeThresholdImageFilter< TInput, TOutput >::
GenerateData()
{
  m_InputCast->SetInput( this->GetInput() );
  m_InputCast->Update();
  ImagePointer input = m_InputCast->GetOutput();

  ImageSpacingType spacing = input->GetSpacing();
  ImagePointType origin = input->GetOrigin();

  // convert the labeled image in a collection of label objects
  ShapeConverterPointer shapeConverter = ShapeConverterType::New();
  shapeConverter->SetInput ( input );
  shapeConverter->SetBackgroundValue ( 0 );
  shapeConverter->Update();

  ShapeLabelMapPointer labelMap = shapeConverter->GetOutput();

  m_Map2Image->SetInput( labelMap );
  m_Map2Image->Update();

  ImagePointer comp = m_Map2Image->GetOutput();

  //std::cout << labelMap->GetNumberOfLabelObjects() << std::endl;

  unsigned int maxSize = 0, label = 0, maxSizeLabel = 0;

  if ( m_MaxComponent )
  {
    // Identify the maximum size
    for( unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); i++ )
    {
      label = labelMap->GetLabels()[i];
      const ShapeLabelObjectType *labelObject = labelMap->GetLabelObject ( label );

      if ( labelObject->Size() > maxSize )
      {
        maxSize = labelObject->Size();
        maxSizeLabel = label;
      }
    }

    // Make all other labels 0
    IteratorType cIt( comp, comp->GetLargestPossibleRegion() );
    for ( cIt.GoToBegin(); !cIt.IsAtEnd(); ++cIt )
    {
      if ( cIt.Get() != maxSizeLabel )
      {
        cIt.Set( 0 );
      }
      else
      {
        cIt.Set( m_Label );
      }
    }
  }
  else
  {
    for( unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); i++ )
    {
      label = labelMap->GetLabels()[i];
      const ShapeLabelObjectType *labelObject = labelMap->GetLabelObject ( label );

      if ( ( labelObject->Size() < m_MinComponentSize ) ||
        ( labelObject->Size() > m_MaxComponentSize ) )
      {
        // Update the voronoi assignment
        IteratorType cIt( comp, labelObject->GetBoundingBox() );
        for ( cIt.GoToBegin(); !cIt.IsAtEnd(); ++cIt )
        {
          if ( cIt.Get() == label )
          {
            cIt.Set( 0 );
          }
        }
      }
    }
  }

  RelabelFilterPointer relabel = RelabelFilterType::New();
  relabel->SetInput ( comp );
  relabel->SetBackgroundValue ( 0 );
  relabel->SetReverseOrdering ( false );
  relabel->SetAttribute ( 0 );
  relabel->Update();

  m_OutputCast->SetInput( relabel->GetOutput() );
  m_OutputCast->GraftOutput( this->GetOutput() );
  m_OutputCast->Update();

  this->GraftOutput( m_OutputCast->GetOutput() );
}

template < class TInput, class TOutput >
void
SizeThresholdImageFilter< TInput, TOutput >::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Minimum ComponentSize: " << GetMinComponentSize() <<
    std::endl;
  os << indent << "Maximum ComponentSize: " << GetMaxComponentSize() <<
    std::endl;
}

} /* end namespace itk */

#endif
