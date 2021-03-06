#ifndef INCLUDED_QSS_TVD_RUNGE_KUTTA_3D_H
#define INCLUDED_QSS_TVD_RUNGE_KUTTA_3D_H

#include "QSSLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_tvd_runge_kutta3d.h
 *
 * \brief
 * @ref lsm_tvd_runge_kutta3d.h provides support for time integration of
 * partial differential equations in three space dimensions via 
 * total-variation diminishing Runge-Kutta methods.  Support is provided 
 * for first-, second-, and third-order time integration.
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
#define QSS3D_RK1_STEP                      qss3drk1step_


/*!
 * LSM3D_RK1_STEP() takes a single first-order Runge-Kutta (i.e. Forward
 * Euler) step.
 *
 * Arguments:
 *  - u_next (out):  u(t_cur+dt)
 *  - u_cur (in):    u(t_cur)
 *  - rhs (in):      right-hand side of time evolution equation
 *  - dt (in):       step size
 *  - *_gb (in):     index range for ghostbox
 *  - *_fb (in):     index range for fillbox
 *
 * Return value:     none
 */
void QSS3D_RK1_STEP(
  QSSLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const QSSLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const QSSLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const QSSLIB_REAL *dt);

#endif
