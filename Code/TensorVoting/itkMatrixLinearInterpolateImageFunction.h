/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMatrixLinearInterpolateImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2009-10-29 11:18:45 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMatrixLinearInterpolateImageFunction_h
#define __itkMatrixLinearInterpolateImageFunction_h

#include "itkMatrixInterpolateImageFunction.h"

namespace itk
{

/**
 * \class MatrixLinearInterpolateImageFunction
 * \brief Linearly interpolate a matrix image at specified positions.
 *
 * MatrixLinearInterpolateImageFunction linearly interpolates a vector
 * image intensity non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type.
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for Matrix images. For
 * scalar images use LinearInterpolateImageFunction.
 *
 * \ingroup ImageFunctions ImageInterpolators
 *
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT MatrixLinearInterpolateImageFunction :
  public MatrixInterpolateImageFunction<TInputImage,TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef MatrixLinearInterpolateImageFunction                  Self;
  typedef MatrixInterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MatrixLinearInterpolateImageFunction,
    MatrixInterpolateImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::ValueType      ValueType;
  typedef typename Superclass::RealType       RealType;

  /** Grab the vector dimension from the superclass. */
  itkStaticConstMacro(Dimension, unsigned int,
                       Superclass::Dimension);

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType       IndexType;
  typedef typename IndexType::IndexValueType  IndexValueType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is Matrix<double,Dimension,Dimension> */
  typedef typename Superclass::OutputType OutputType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & index ) const;

protected:
  MatrixLinearInterpolateImageFunction();
  ~MatrixLinearInterpolateImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  MatrixLinearInterpolateImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;

};

} // end namespace itk

# include "itkMatrixLinearInterpolateImageFunction.txx"

#endif
