#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkSampleClassifierFilter.h"
#include <itkDiscreteGaussianImageFilter.h>
#include "itkBinaryThinningImageFilter3D.h"
#include "itkCastImageFilter.h"

#include "vnl/vnl_random.h"
#include "itkShotNoiseImageFilter.h"
#include "itkSpeckleNoiseImageFilter.h"
#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include <itkNeighborhoodIterator.h>
#include "itkZeroFluxNeumannBoundaryCondition.h"

int main(int argc, char *argv[])
{
  if ( argc < 5 )
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " input outputSegmentation outputWithNoise sigma s" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;

  typedef itk::Image< unsigned char, Dimension > FeatureImageType;
  typedef itk::Image< float, Dimension > InputImageType;
  typedef itk::Image< unsigned int, Dimension >  SegmentImageType;
  typedef itk::ImageFileReader< InputImageType > InputReaderType;
  typedef itk::ImageFileWriter< SegmentImageType > SegmentWriterType;
  typedef itk::ImageFileWriter< FeatureImageType > FeatureWriterType;

  typedef itk::ImageRegionIteratorWithIndex< InputImageType > InputIteratorType;
  typedef itk::Vector<float, Dimension > MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
  typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
  typedef TreeGeneratorType::KdTreeType TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
  typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
  typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
  typedef itk::Statistics::SampleClassifierFilter<SampleType> ClassifierType;
  typedef ClassifierType::ClassLabelVectorObjectType ClassLabelVectorObjectType;
  typedef ClassifierType::ClassLabelVectorType ClassLabelVectorType;
  typedef ClassifierType::ClassLabelType ClassLabelType;
  typedef ClassifierType::MembershipFunctionVectorObjectType MembershipFunctionVectorObjectType;
  typedef ClassifierType::MembershipFunctionVectorType MembershipFunctionVectorType;

  typedef itk::CastImageFilter< InputImageType, SegmentImageType > CastFilterType;
  typedef itk::DiscreteGaussianImageFilter< InputImageType, InputImageType > SmoothingFilterType;
  typedef itk::RecursiveGaussianImageFilter< InputImageType, InputImageType > PSFFilterType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
  typedef itk::RescaleIntensityImageFilter< InputImageType, FeatureImageType > RescaleFilterType;
  typedef itk::ShotNoiseImageFilter< FeatureImageType, FeatureImageType > PoissonFilterType;
  typedef itk::AdditiveGaussianNoiseImageFilter< InputImageType, InputImageType > NormalFilterType;
  typedef itk::ResampleImageFilter< SegmentImageType, SegmentImageType > ResampleFilterType;
  typedef itk::LinearInterpolateImageFunction< SegmentImageType > InterpolatorType;
  typedef itk::AffineTransform< double, Dimension > TransformType;
  typedef itk::LinearInterpolateImageFunction< InputImageType > InputInterpolatorType;
  typedef itk::ResampleImageFilter< InputImageType, InputImageType > InputResampleFilterType;

  InputImageType::SpacingType spacing;
  spacing[0] = spacing[1] = 0.2;
  spacing[2] = 0.2;

  InputImageType::RegionType region;
  InputImageType::SizeType size;
  size[0] = size[1] = 256;
  size[2] = 256;
  InputImageType::PointType origin1;
  origin1[0] = origin1[1] = origin1[2] = 0.0;
  InputImageType::IndexType index;
  index[0] = index[1] = index[2] = 0;

  region.SetIndex( index );
  region.SetSize( size );

  SampleType::Pointer sample = SampleType::New();
  MeasurementVectorType mv;

  // Allocate the output image
  InputImageType::Pointer outputImage = InputImageType::New();
  outputImage->SetRegions( region );
  outputImage->SetSpacing( spacing );
  outputImage->SetOrigin( origin1 );
  outputImage->Allocate();

  InputImageType::PointType p;
  InputIteratorType It( outputImage, outputImage->GetLargestPossibleRegion() );
  It.GoToBegin();
  while( !It.IsAtEnd() )
  {
    index = It.GetIndex();
    outputImage->TransformIndexToPhysicalPoint( index, p );
    mv[0] = static_cast< float > ( p[0] );
    mv[1] = static_cast< float > ( p[1] );
    mv[2] = static_cast< float > ( p[2] );
    sample->PushBack(mv);
    ++It;
  }

  TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample(sample);
  treeGenerator->SetBucketSize(16);
  treeGenerator->Update();

  EstimatorType::Pointer estimator = EstimatorType::New();

  vnl_random r;

  unsigned int numOfClasses = atoi( argv[1] );
  EstimatorType::ParametersType initialMeans( numOfClasses * Dimension );
  for( unsigned int i = 0; i < numOfClasses; i++ )
  {
    for( unsigned int j = 0; j < Dimension; j++ )
    {
      initialMeans[ i*Dimension+j ] = r.drand32(0, 1) * size[j] * spacing[j];
    }
  }

  estimator->SetParameters(initialMeans);
  estimator->SetKdTree(treeGenerator->GetOutput());
  estimator->SetMaximumIteration( 20 );
  estimator->SetCentroidPositionChangesThreshold( 2.5 );
  estimator->StartOptimization();

  EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();

  DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
  ClassifierType::Pointer classifier = ClassifierType::New();
  classifier->SetDecisionRule( decisionRule);
  classifier->SetInput( sample );
  classifier->SetNumberOfClasses( numOfClasses );

  ClassLabelVectorObjectType::Pointer classLabelsObject = ClassLabelVectorObjectType::New();
  ClassLabelVectorType& classLabelsVector = classLabelsObject->Get();

  for( unsigned int i = 1; i <= numOfClasses; i++ )
  {
    ClassLabelType class1 = i;
    classLabelsVector.push_back( class1 );
  }
  classifier->SetClassLabels( classLabelsObject );

  MembershipFunctionVectorObjectType::Pointer membershipFunctionVectorObject =
    MembershipFunctionVectorObjectType::New();
  MembershipFunctionVectorType& membershipFunctionVector =
    membershipFunctionVectorObject->Get();

  int indexCount = 0;
  for ( unsigned int i = 0 ; i < numOfClasses; i++ )
    {
    MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
    MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
    for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ )
      {
      centroid[j] = estimatedMeans[indexCount++];
      }
    membershipFunction->SetCentroid( centroid );
    membershipFunctionVector.push_back( membershipFunction.GetPointer() );
    }
  classifier->SetMembershipFunctions( membershipFunctionVectorObject );

  try
  {
    classifier->Update();
  }
  catch ( itk::ExceptionObject e )
  {
    std::cerr << "Exception: " << e << std::endl;
  }
  std::cout << "Classification complete..." << std::endl;

  const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
  ClassifierType::MembershipSampleType::ConstIterator iter = membershipSample->Begin();

  It.GoToBegin();
  while ( iter != membershipSample->End() )
    {
    It.Set( iter.GetClassLabel() );
    ++iter;
    ++It;
    }

  // Z sampling
  size[2] = size[2]/5;
  spacing[2] = spacing[2]*5;

  CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput( outputImage );
  caster->Update();
  std::cout << "Casting complete..." << std::endl;

  // Resample and write out here
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  InterpolatorType::Pointer interp = InterpolatorType::New();

  ResampleFilterType::Pointer resample1 = ResampleFilterType::New();
  resample1->SetTransform ( transform );
  resample1->SetInterpolator ( interp );
  resample1->SetInput ( caster->GetOutput() );
  resample1->SetSize ( size );
  resample1->SetOutputOrigin ( origin1 );
  resample1->SetOutputSpacing ( spacing );
  resample1->SetDefaultPixelValue ( 0 );
  resample1->Update();
  std::cout << "Resampling complete..." << std::endl;

  SegmentWriterType::Pointer writer1 = SegmentWriterType::New();
  writer1->SetFileName(argv[2]);
  writer1->SetInput( resample1->GetOutput() );
  writer1->Update();
  std::cout << "Manual segmentation complete..." << std::endl;

  itk::ZeroFluxNeumannBoundaryCondition<InputImageType> neumann;
  typedef itk::NeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType::RadiusType radius = {1, 1, 1};
  NeighborhoodIteratorType it( radius, outputImage, outputImage->GetLargestPossibleRegion() );
  it.SetBoundaryCondition( neumann );

  InputImageType::Pointer output = InputImageType::New();
  output->CopyInformation( outputImage );
  output->SetRegions( outputImage->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer(0);

  IteratorType oIt( output, output->GetLargestPossibleRegion() );
  oIt.GoToBegin();
  float m,n;
  for ( it.GoToBegin(), oIt.GoToBegin(); !it.IsAtEnd(); ++it, ++oIt)
  {
    m = it.GetCenterPixel();
    for (unsigned int i = 0; i < 27; i++)//27 is the size of neighborhood
    {
      n = it.GetPixel(i);
      if ( m > n )
        oIt.Set( 1 );
    }
  }
  outputImage = 0;

  PSFFilterType::Pointer psfx = PSFFilterType::New();
  psfx->SetInput( output );
  psfx->SetDirection( 0 );
  psfx->SetOrder( PSFFilterType::ZeroOrder );
  psfx->SetNormalizeAcrossScale( false );
  psfx->SetSigma( 0.2 );
  psfx->Update();

  PSFFilterType::Pointer psfy = PSFFilterType::New();
  psfy->SetInput( psfx->GetOutput() );
  psfy->SetDirection( 1 );
  psfy->SetOrder( PSFFilterType::ZeroOrder );
  psfy->SetNormalizeAcrossScale( false );
  psfy->SetSigma( 0.2 );
  psfy->Update();

  PSFFilterType::Pointer psfz = PSFFilterType::New();
  psfz->SetInput( psfy->GetOutput() );
  psfz->SetDirection( 2 );
  psfz->SetOrder( PSFFilterType::ZeroOrder );
  psfz->SetNormalizeAcrossScale( false );
  psfz->SetSigma( 0.4 );
  psfz->Update();

  NormalFilterType::Pointer gauss = NormalFilterType::New();
  gauss->SetInput( psfz->GetOutput() );
  gauss->SetStandardDeviation( atof(argv[4]) );
  gauss->SetMean( 0 );
  gauss->Update();

  // Resample and write out here
  InputInterpolatorType::Pointer interp2 = InputInterpolatorType::New();
  InputResampleFilterType::Pointer resample2 = InputResampleFilterType::New();
  resample2->SetTransform ( transform );
  resample2->SetInterpolator ( interp2 );
  resample2->SetInput ( gauss->GetOutput() );
  resample2->SetSize ( size );
  resample2->SetOutputOrigin ( origin1 );
  resample2->SetOutputSpacing ( spacing );
  resample2->SetDefaultPixelValue ( 0 );
  resample2->Update();

  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( resample2->GetOutput() );//
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();

  PoissonFilterType::Pointer poisson = PoissonFilterType::New();
  poisson->SetInput( rescale->GetOutput() );
  poisson->SetScale( atof(argv[5]) );
  poisson->Update();

  FeatureWriterType::Pointer writer2 = FeatureWriterType::New();
  writer2->SetFileName(argv[3]);
  writer2->SetInput( poisson->GetOutput() );
  writer2->Write();
  writer2->Update();

  return 0;
}

