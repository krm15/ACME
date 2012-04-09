/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 756 $  // Revision of last commit
  Date: $Date: 2009-10-20 11:50:49 -0400 (Tue, 20 Oct 2009) $  // Date of last commit
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

#ifndef __itkConfocalMembraneImageGenerator_h
#define __itkConfocalMembraneImageGenerator_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageToListAdaptor.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkEuclideanDistance.h"
#include "itkSampleClassifier.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "vnl/vnl_random.h"
#include "itkShotNoiseImageFilter.h"
#include "itkSpeckleNoiseImageFilter.h"
#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSaltAndPepperNoiseImageFilter.h"

namespace itk
{
template < class TInputImage >
class ITK_EXPORT ConfocalMembraneImageGenerator : public Object
{
  public:
    typedef ConfocalMembraneImageGenerator  Self;
    typedef Object                          Superclass;
    typedef SmartPointer< Self >            Pointer;
    typedef SmartPointer< const Self >      ConstPointer;

    itkStaticConstMacro ( ImageDimension, unsigned int,
                          TInputImage::ImageDimension );

    /** Method for creation through object factory */
    itkNewMacro ( Self );

    /** Run-time type information */
    itkTypeMacro ( ConfocalMembraneImageGenerator, Object );

    /** Display */
    void PrintSelf ( std::ostream& os, Indent indent ) const;

    typedef TInputImage                      ImageType;
    typedef typename ImageType::Pointer      ImagePointer;
    typedef typename ImageType::ConstPointer ImageConstPointer;
    typedef typename ImageType::PixelType    ImagePixelType;
    typedef typename ImageType::RegionType   ImageRegionType;
    typedef typename ImageType::SizeType     ImageSizeType;
    typedef typename ImageSizeType::SizeValueType ImageSizeValueType;
    typedef typename ImageType::SpacingType  ImageSpacingType;
    typedef typename ImageType::IndexType    ImageIndexType;
    typedef typename ImageType::PointType    ImagePointType;

    typedef Vector< float, ImageDimension > MeasurementVectorType;
    typedef Statistics::ListSample< MeasurementVectorType > SampleType;
    typedef typename SampleType::Pointer SamplePointer;
    typedef Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
    typedef typename TreeGeneratorType::Pointer TreeGeneratorPointer;
    typedef typename TreeGeneratorType::KdTreeType TreeType;
    typedef Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
    typedef typename EstimatorType::Pointer EstimatorPointer;
    typedef Statistics::EuclideanDistance< MeasurementVectorType > MembershipFunctionType;
    typedef typename MembershipFunctionType::Pointer MembershipFunctionPointer;
    typedef MinimumDecisionRule DecisionRuleType;
    typedef typename DecisionRuleType::Pointer DecisionRulePointer;
    typedef Statistics::SampleClassifier< SampleType > ClassifierType;
    typedef typename ClassifierType::Pointer ClassifierPointer;

    typedef DiscreteGaussianImageFilter< ImageType, ImageType >  PSFFilterType;
    typedef typename PSFFilterType::Pointer PSFFilterPointer;
    typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;
    typedef CastImageFilter< ImageType, ImageType > CastFilterType;
    typedef typename CastFilterType::Pointer CastFilterPointer;
    typedef ShotNoiseImageFilter< ImageType, ImageType > PoissonFilterType;
    typedef typename PoissonFilterType::Pointer PoissonFilterPointer;
    typedef SpeckleNoiseImageFilter< ImageType, ImageType > SpeckleFilterType;
    typedef typename SpeckleFilterType::Pointer SpeckleFilterPointer;
    typedef AdditiveGaussianNoiseImageFilter< ImageType, ImageType > NormalFilterType;
    typedef typename NormalFilterType::Pointer NormalFilterPointer;
    typedef SaltAndPepperNoiseImageFilter< ImageType, ImageType > SPFilterType;
    typedef typename SPFilterType::Pointer SPFilterPointer;

    itkGetConstMacro ( NumberOfClasses, unsigned int );
    itkSetMacro ( NumberOfClasses, unsigned int );

    void SetMembraneImage ( ImagePointer p )
    {
      p = m_Membrane;
    }

    ImagePointer GetMembraneImage ( void )
    {
      GenerateVoronoiTesselation();
      ComputeNoiseModel();
      return m_Membrane;
    }

    void SetOrigin( ImagePointType& origin )
    {
      m_Origin = origin;
    }

    void SetSize( ImageSizeType& size )
    {
      m_Size = size;
    }

    void SetSpacing( ImageSpacingType& spacing )
    {
      m_Spacing = spacing;
    }

  protected:
    ConfocalMembraneImageGenerator();
    ~ConfocalMembraneImageGenerator() {}

    void GenerateVoronoiTesselation();
    void ComputeNoiseModel();

    unsigned int m_NumberOfClasses;
    ImageSpacingType m_Spacing;
    ImagePointType m_Origin;
    ImageSizeType m_Size;
    ImagePointer m_Membrane;

  private:
    ConfocalMembraneImageGenerator ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkConfocalMembraneImageGenerator.txx"
#endif
