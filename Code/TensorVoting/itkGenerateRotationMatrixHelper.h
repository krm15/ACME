/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGenerateRotationMatrixHelper.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:32 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGenerateRotationMatrixHelper_h
#define __itkGenerateRotationMatrixHelper_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{

/** \class GenerateRotationMatrixHelper
 * This calculator generates the rotation matrix.
 * It is templated over the output image type.
 *
 * \ingroup Operators
 */
template <class TOutputImage>
class ITK_EXPORT GenerateRotationMatrixHelper : public Object
{
public:
  /** Standard class typedefs. */
  typedef GenerateRotationMatrixHelper  Self;
  typedef Object                        Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GenerateRotationMatrixHelper, Object);

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  /** Type definition for the output image. */
  typedef TOutputImage                       ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::PixelType      PixelType;
  typedef typename ImageType::IndexType      IndexType;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename ImageType::PointType      PointType;
  typedef typename ImageType::SpacingType    SpacingType;

  typedef vnl_vector< double > vnlVectorType;
  typedef Vector< double, ImageDimension > VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension> MatrixType;
  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;

  void ComputeRotationMatrixWithAxis( double ux, double uy, double uz,
                                      double theta, MatrixType& R )
  {
    R[0][0] = ux*ux + ( 1 - ux*ux )*vcl_cos( theta );
    R[0][1] = ux*uy*( 1 - vcl_cos(theta) ) - uz*vcl_sin( theta );
    R[0][2] = ux*uz*( 1 - vcl_cos(theta) ) + uy*vcl_sin( theta );

    R[1][0] = ux*uy*( 1 - vcl_cos(theta) ) + uz*vcl_sin( theta );
    R[1][1] = uy*uy + ( 1 - uy*uy )*vcl_cos( theta );
    R[1][2] = uy*uz*( 1 - vcl_cos(theta) ) - ux*vcl_sin( theta );


    R[2][0] = ux*uz*( 1 - vcl_cos(theta) ) - uy*vcl_sin( theta );
    R[2][1] = uy*uz*( 1 - vcl_cos(theta) ) + ux*vcl_sin( theta );
    R[2][2] = uz*uz + ( 1 - uz*uz )*vcl_cos( theta );
  }

  void ComputeRotationMatrix( VectorType& u, MatrixType& R )
  {
    VectorType zAxis, rAxis;
    zAxis.Fill( 0 );
    zAxis[ImageDimension-1] = 1;
    rAxis = CrossProduct( zAxis, u );

    if ( rAxis.GetVnlVector().is_zero() )
    {
      R.SetIdentity();
      return;
    }

    rAxis.Normalize();

    double c = zAxis * u;
    double theta = vcl_acos(  c );

    ComputeRotationMatrixWithAxis( rAxis[0], rAxis[1], rAxis[2], theta, R );

    vnlVectorType p, q;
    p = R.GetVnlMatrix() * u.GetVnlVector();
    q = p - zAxis.GetVnlVector();
    if ( q.magnitude() > 0.01 )
    {
      ComputeRotationMatrixWithAxis( rAxis[0], rAxis[1], rAxis[2], -theta, R );
    }
  }

  GenerateRotationMatrixHelper(){}
  virtual ~GenerateRotationMatrixHelper() {};

protected:

private:
  GenerateRotationMatrixHelper(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#endif /* __itkGenerateRotationMatrixHelper_h */
