/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStickFieldGenerator2D.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:32 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorVoting3D_h
#define __itkTensorVoting3D_h

#include "itkTensorVoting.h"
#include "itkThreadSafeMersenneTwisterRandomVariateGenerator.h"
#include "itkBallFieldGenerator3D.h"
#include "itkPlateFieldGenerator3D.h"
#include "itkStickFieldGenerator3D.h"
#include "itkTensorToSaliencyImageFilter.h"
#include "itkComposeVotesFromLookupImageFilter3D.h"
#include <itkAzimuthElevationToCartesianTransform.h>

// #include "itkImageFileWriter.h"

namespace itk
{

/** \class TensorVoting3D
 * This calculator performs tensor voting in 3D.
 * It is templated over the input image type.
 *
 * \ingroup Operators
 */
template <class TInputImage >
class ITK_EXPORT TensorVoting3D :
 public TensorVoting< TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef TensorVoting3D              Self;
  typedef TensorVoting<TInputImage >  Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TensorVoting3D, TensorVoting );

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   ImageConstPointer;
  typedef typename InputImageType::PixelType      PixelType;
  typedef typename InputImageType::IndexType      IndexType;
  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType        SizeValueType;
  typedef typename InputImageType::PointType      PointType;
  typedef typename InputImageType::SpacingType    SpacingType;

  typedef StickFieldGenerator3D< InputImageType >   StickGeneratorType;
  typedef typename StickGeneratorType::Pointer      StickGeneratorPointer;
  typedef PlateFieldGenerator3D< InputImageType >   PlateGeneratorType;
  typedef typename PlateGeneratorType::Pointer      PlateGeneratorPointer;
  typedef BallFieldGenerator3D< InputImageType >    BallGeneratorType;
  typedef typename BallGeneratorType::Pointer       BallGeneratorPointer;

  typedef GenerateRotationMatrixHelper< InputImageType > RotationMatrixHelperType;
  typedef AzimuthElevationToCartesianTransform< double, ImageDimension > CoordinateTransformType;
  typedef typename CoordinateTransformType::Pointer CoordinateTransformPointer;

  /** Type definition for the output image. */
  typedef Image< bool, ImageDimension >      TokenImageType;
  typedef typename TokenImageType::Pointer   TokenImagePointer;
  typedef ImageRegionIteratorWithIndex< TokenImageType > TokenIteratorType;

  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionConstIteratorWithIndex< InputImageType > ConstIteratorType;
  typedef ImageRegionIteratorWithIndex< InputImageType > IteratorType;

  typedef Image< VectorType, ImageDimension > VectorImageType;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  typedef ImageRegionIterator< VectorImageType > VectorIteratorType;

  typedef std::list< IndexType > IdType;
  typedef Vector< IdType, ImageDimension > VectorInternalType;
  typedef Image<VectorInternalType, ImageDimension> InternalImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;
  typedef ImageRegionIteratorWithIndex< InternalImageType > InternalIteratorType;
  
  typedef Image< double, ImageDimension > DoubleImageType;
  typedef typename DoubleImageType::Pointer DoubleImagePointer;
  typedef TensorToSaliencyImageFilter< InputImageType, VectorImageType > SaliencyFilterType;
  typedef VectorIndexSelectionCastImageFilter< VectorImageType, DoubleImageType > IndexFilterType;
  typedef ImageRegionIteratorWithIndex< DoubleImageType > DoubleIteratorType;

  typedef ComposeVotesFromLookupImageFilter3D< InternalImageType >
    ComposeVotesFilterType;
  typedef typename ComposeVotesFilterType::Pointer ComposeVotesFilterPointer;

  typedef Statistics::ThreadSafeMersenneTwisterRandomVariateGenerator
    RandomGeneratorType;
  typedef typename RandomGeneratorType::Pointer RandomGeneratorPointer;

protected:
  TensorVoting3D();
  virtual ~TensorVoting3D() {}

  void InitializeVotingFields();
  void ComputeLookup();
  void IntegrateVotes();

private:
  TensorVoting3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorVoting3D.txx"
#endif

#endif /* __itkTensorVoting3D_h */
