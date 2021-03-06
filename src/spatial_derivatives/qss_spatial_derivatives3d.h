#ifndef INCLUDED_QSS_SPATIAL_DERIVATIVES_3D_H
#define INCLUDED_QSS_SPATIAL_DERIVATIVES_3D_H

#include "QSSLIB_config.h"

/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define QSS3D_HJ_ENO2_TEST_NEW            qss3dhjeno2testnew_
#define QSS3D_CENTRAL_GRAD_ORDER2         qss3dcentralgradorder2_
#define QSS3D_GET_RHS_SECOND              qss3dgetrhs2nd_
#define QSS3D_SET_VAR_CURV_ADV            qss3dsetvarcurvadv_
#define QSS3D_SET_VAR_NORM                qss3dsetvarnorm_
#define QSS3D_SIGNED_LINEAR_EXTRAPOLATION  qss3dsignedlinearextrapolation_


 void QSS3D_HJ_ENO2_TEST_NEW(
   QSSLIB_REAL *phi,
   QSSLIB_REAL *lse_rhs,
   QSSLIB_REAL *vel_n,
   const int *ilo_phi_gb,
   const int *ihi_phi_gb, 
   const int *jlo_phi_gb,
   const int *jhi_phi_gb,
   const int *klo_phi_gb,
   const int *khi_phi_gb,
   const int *ilo_fb, 
   const int *ihi_fb,
   const int *jlo_fb,
   const int *jhi_fb,
   const int *klo_fb,
   const int *khi_fb,
   const QSSLIB_REAL  *dx,
   const QSSLIB_REAL  *dy,
   const QSSLIB_REAL  *dz);
     
/*! 
 * LSM3D_CENTRAL_GRAD_ORDER2() computes the second-order, central, 
 * finite difference approximation to the gradient of \f$ \phi \f$ 
 * using the formula:
 * 
 *    \f[
 *
 *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
 *        \frac{ \phi_{i+1} - \phi_{i-1} }{ 2 dx }
 *
 *    \f]
 *
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 */
void QSS3D_CENTRAL_GRAD_ORDER2( 
  QSSLIB_REAL *phi_x,
  QSSLIB_REAL *phi_y,
  QSSLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const QSSLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const QSSLIB_REAL *dx,
  const QSSLIB_REAL *dy,
  const QSSLIB_REAL *dz);


void QSS3D_GET_RHS_SECOND(
   QSSLIB_REAL *phi,
   QSSLIB_REAL *vel_n,
   QSSLIB_REAL *b,
   QSSLIB_REAL *vel_x, 
   QSSLIB_REAL *vel_y, 
   QSSLIB_REAL *vel_z,
   QSSLIB_REAL *max_H_over_dX,
   QSSLIB_REAL *lse_rhs,
   const int *ilo_phi_gb,
   const int *ihi_phi_gb, 
   const int *jlo_phi_gb,
   const int *jhi_phi_gb,
   const int *klo_phi_gb,
   const int *khi_phi_gb, 
   const int *ilo_fb, 
   const int *ihi_fb,
   const int *jlo_fb,
   const int *jhi_fb,
   const int *klo_fb,
   const int *khi_fb,
   const QSSLIB_REAL  *dx,
   const QSSLIB_REAL  *dy,
   const QSSLIB_REAL  *dz);  


void QSS3D_SET_VAR_CURV_ADV(
    QSSLIB_REAL      *mask,
    QSSLIB_REAL      *mask_x,
    QSSLIB_REAL      *mask_y,
    QSSLIB_REAL      *mask_z,
    QSSLIB_REAL      *var_b,
    QSSLIB_REAL      *vel_x,
    QSSLIB_REAL      *vel_y,
    QSSLIB_REAL      *vel_z,
    QSSLIB_REAL      *b_max_over_dx,
    QSSLIB_REAL      *max_U_over_dx,
    QSSLIB_REAL      *const_b,
    const int *ilo_gb, 
    const int *ihi_gb,
    const int *jlo_gb, 
    const int *jhi_gb,
    const int *klo_gb, 
    const int *khi_gb,
    const int *ilo_fb,
    const int *ihi_fb,
    const int *jlo_fb,
    const int *jhi_fb,
    const int *klo_fb,
    const int *khi_fb,
    const QSSLIB_REAL *dx);
    
/* The variational normal velocity must be calculated each step of the da loop.
    So it is more efficient to do it separately, and not calculate the mask_x etc.
    derivatives at each step
*/
void QSS3D_SET_VAR_NORM(
    QSSLIB_REAL      *mask,
    QSSLIB_REAL      *mask_x,
    QSSLIB_REAL      *mask_y,
    QSSLIB_REAL      *mask_z,
    QSSLIB_REAL      *var_a,
    QSSLIB_REAL      *a0,
    QSSLIB_REAL      *ca,
    const int *ilo_gb, 
    const int *ihi_gb,
    const int *jlo_gb, 
    const int *jhi_gb,
    const int *klo_gb, 
    const int *khi_gb,
    const int *ilo_fb,
    const int *ihi_fb,
    const int *jlo_fb,
    const int *jhi_fb,
    const int *klo_fb,
    const int *khi_fb,
    const QSSLIB_REAL *dx);
    
void QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
  QSSLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *jlo_gb, const int *jhi_gb,
  const int *klo_gb, const int *khi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *jlo_fb, const int *jhi_fb,
  const int *klo_fb, const int *khi_fb,
  const int *bdry_location_idx);
  
 
   
  #endif
