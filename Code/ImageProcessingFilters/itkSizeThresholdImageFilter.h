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

#ifndef __itkSizeThresholdImageFilter_h
#define __itkSizeThresholdImageFilter_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageToImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkShapeRelabelImageFilter.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkVector.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

#include <vnl/vnl_vector.h>
#include "itkNumericTraits.h"

namespace itk {

template < class TInput, class TOutput >
class ITK_EXPORT SizeThresholdImageFilter : public ImageToImageFilter<
TInput, TOutput >
{
public:

  typedef SizeThresholdImageFilter             Self;
  typedef ImageToImageFilter< TInput, TOutput > Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int,
    TInput::ImageDimension);

  /** Method for creation through object factory */
  itkNewMacro(Self);

  /** Run-time type information */
  itkTypeMacro( SizeThresholdImageFilter, ImageToImageFilter );

  /** Display */
  void PrintSelf( std::ostream& os, Indent indent ) const;

  typedef Image< unsigned long,ImageDimension >  ImageType;
  typedef typename ImageType::Pointer             ImagePointer;
  typedef typename ImageType::ConstPointer        ImageConstPointer;
  typedef typename ImageType::PixelType           ImagePixelType;
  typedef typename ImageType::RegionType          ImageRegionType;
  typedef typename ImageType::SizeType            ImageSizeType;
  typedef typename ImageSizeType::SizeValueType   ImageSizeValueType;
  typedef typename ImageType::IndexType           ImageIndexType;
  typedef typename ImageIndexType::IndexValueType ImageIndexValueType;
  typedef typename ImageType::SpacingType         ImageSpacingType;
  typedef typename ImageType::PointType           ImagePointType;

  typedef ImageRegionIterator< ImageType > IteratorType;
  typedef CastImageFilter< TInput, ImageType > InputCastType;
  typedef CastImageFilter< ImageType, TOutput > OutputCastType;
  typedef BinaryThresholdImageFilter< ImageType, ImageType >
    ThresholdFilterType;
  typedef typename ThresholdFilterType::Pointer ThresholdFilterPointer;

  typedef unsigned int                                  LabelType;
  typedef ShapeLabelObject< LabelType, ImageDimension > ShapeLabelObjectType;
  typedef LabelMap< ShapeLabelObjectType >              ShapeLabelMapType;
  typedef typename ShapeLabelMapType::Pointer ShapeLabelMapPointer;
  typedef LabelImageToShapeLabelMapFilter< ImageType, ShapeLabelMapType >  ShapeConverterType;
  typedef typename ShapeConverterType::Pointer ShapeConverterPointer;
  typedef ShapeLabelMapFilter< ShapeLabelMapType > ShapeFilterType;
  typedef ShapeRelabelImageFilter< ImageType >     RelabelFilterType;
  typedef typename RelabelFilterType::Pointer RelabelFilterPointer;
  typedef LabelMapToLabelImageFilter< ShapeLabelMapType, ImageType > MapToImageType;

  itkGetConstMacro( MinComponentSize, unsigned int );
  itkSetMacro( MinComponentSize, unsigned int );
  itkGetConstMacro( MaxComponentSize, unsigned int );
  itkSetMacro( MaxComponentSize, unsigned int );
  itkGetConstMacro( MaxComponent, bool );
  itkSetMacro( MaxComponent, bool );
  itkGetConstMacro( Label, unsigned int );
  itkSetMacro( Label, unsigned int );

protected:

  SizeThresholdImageFilter();
  ~SizeThresholdImageFilter() {}
  void GenerateData();

  typename InputCastType::Pointer m_InputCast;
  typename OutputCastType::Pointer m_OutputCast;
  typename MapToImageType::Pointer  m_Map2Image;

  unsigned int m_MinComponentSize;
  unsigned int m_MaxComponentSize;
  bool m_MaxComponent;
  unsigned int m_Label;

private:

  SizeThresholdImageFilter(Self&);   // intentionally not implemented
  void operator=(const Self&);          // intentionally not implemented
};

} /* namespace itk */

#include "itkSizeThresholdImageFilter.txx"
#endif
