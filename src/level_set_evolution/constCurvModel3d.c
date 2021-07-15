 /* Include constant curvature drivers for all dimensions here */
 
/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Local headers */
#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_reinitialization3d.h"
#include "qss_grid.h"
#include "qss_macros.h"
#include "constCurvModel3d.h"
#include "qss_general_util.h"
#include "connectivity.h"

#define DSZ  sizeof(QSSLIB_REAL)

/* Main driver for constant curvature level set method model */

QSSLIB_REAL constCurvModel3d(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out) 
{

    QSSLIB_REAL     zero = 0.0;
    char fname[256];
  
    int flag = 0, OUTER_STEP = 0, INNER_STEP, idx;     
  
    QSSLIB_REAL dt = 0, dt_sub, t = 0;
    QSSLIB_REAL   mask_sign = -1;
                                                                                                    
    QSSLIB_REAL eps, cur_max_H_over_dX = -1, cfl_number = 0.9;
    QSSLIB_REAL max_abs_err = 100, vol_phi = 100, vol_very_small;
        
    vol_very_small = (g->dx)[0]*(g->dx)[1]*(g->dx)[2];    
    eps = (options->eps_coefficient)*(g->dx[0]);

    QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]), &(g->dx[2]), &eps);
		        
    while( (t < options->tmax)  && (max_abs_err > options->eps_stop) && (vol_phi > vol_very_small))
    { 
        /* outer loop */
        OUTER_STEP++;
        dt_sub = 0;
        COPY_DATA(p->phi_prev,p->phi,g)
        
        /* Begin parallel region */
        #pragma omp parallel default(none) shared(p, g, flag, cfl_number,\
            cur_max_H_over_dX, zero, dt_sub, vol_phi, max_abs_err, eps, options) \
            private(INNER_STEP, dt)
        {

            INNER_STEP = 0;
            QSSLIB_REAL max_H_over_dX;
            int    bdry_location_idx = 9; /* all boundaries */
            QSSLIB_REAL disconn_overlap = 0;
            
            /* Set up variables for multi-threading */
            int cur_thread, cur_klo_fb, cur_khi_fb, num_threads, nslices, i;
            int cur_klo_gb, cur_khi_gb;
                       
            cur_thread = omp_get_thread_num();
            num_threads = omp_get_num_threads();
    
            if (num_threads > 1) flag = 1;

            nslices = g->khi_fb - g->klo_fb + 1;

            cur_klo_fb = g->klo_fb + nslices*cur_thread/num_threads;
            cur_khi_fb = g->klo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            double t1 = omp_get_wtime();
            /* Keeping track of thread-local ghost boundaries, mainly for imposing mask */
            if (cur_khi_fb > (g->khi_fb))
	            cur_khi_fb = (g->khi_fb);

            if (cur_thread == 0)
                cur_klo_gb = cur_klo_fb - 3;
            else
                cur_klo_gb = cur_klo_fb;
            
            if (cur_thread == (num_threads - 1))
                cur_khi_gb = cur_khi_fb + 3;
            else
                cur_khi_gb = cur_khi_fb;
            

            while( dt_sub < options->tplot )
            { 
                /* inner loop */
                INNER_STEP++;
                      
                QSS3D_GET_RHS_SECOND(p->phi, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y, p->external_velocity_z,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	 
	                
	            #pragma omp barrier
	            
	            #pragma omp critical
	            {
	                if(max_H_over_dX > cur_max_H_over_dX)
	                    cur_max_H_over_dX = max_H_over_dX;
	            }
	
	            /* 
	               Barrier to ensure same cur_max_H_over_dX across all threads.
	            */
	            #pragma omp barrier	            
	            
	            /* get final correct dt due to parabolic (curvature) term */
	            dt =  cfl_number / (cur_max_H_over_dX + options->b_max_over_dx + options->max_U_over_dx);
		

	            QSS3D_RK1_STEP(p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);
		
		        #pragma omp barrier
	            IMPOSE_MASK_PAR(p->phi, p->mask, p->phi_next, &(options->overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            
	            #pragma omp barrier
	            IMPOSE_MASK_PAR(p->phi, p->mask_w, p->phi, &(disconn_overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            
	            #pragma omp barrier
	            IMPOSE_MASK_PAR(p->phi, p->mask_nw, p->phi, &(disconn_overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	              
            	#pragma omp barrier

                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
	            /* 
	               There should be an implicit barrier after omp single,
	               so all threads should sync here
	            */
	             #pragma omp single
	            {
	                SET_DATA_TO_CONSTANT(p->lse_rhs,g,zero);
	                dt_sub += dt;
	            }
	            //printf("%f\n", dt_sub);
	            	    
            }   /* End inner loop */
            
            double t2 = omp_get_wtime();

            if (cur_thread == 0)
	            printf("%d threads: Level set time = %lf\n", num_threads, t2 - t1);
            
         }   /* End parallel region */
            
    /* 
        Reinitialization of the level set function 
        - may want to parallelize the following functions later
    */
        
        t += dt_sub;
        
        double t3 = omp_get_wtime();
        
        if(options->check_connectivity)
            trapComponents_mask(p, g, options);
        
        printf("Reinitializing....");
        reinitialize3d_subcell_fix_qss(p,g,options);
        printf("Reinitialized\n");

            
        /* compute stopping criteria */
        /* max abs error */
        QSS3D_MAX_NORM_DIFF_LOCAL(&max_abs_err,p->phi,GB_DIMS, p->phi_prev,
		        GB_DIMS, FB_DIMS, &(options->err_check_zone));

        QSSLIB_REAL vol_phi_old = vol_phi;
        
        QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]),&(g->dx[2]), &eps);
		
        if(options->checkpoint)
        {
            sprintf(fname,"checkpoint_phi");
            writeDataArrayQSS(p->phi,g,fname,GZIP);
            sprintf(fname,"checkpoint_phi_prev");
            writeDataArrayQSS(p->phi_prev,g,fname,GZIP);
        }
        
        double t4 = omp_get_wtime();
        printf("connectivity time = %lf\n", t4 - t3);

        printf("t = %4.1f\t",t);
        printf("max_abs_err = %4.3f,\t",max_abs_err);
        printf("vol_nwp = %4.3f\n", vol_phi);

        if (fabsf(vol_phi - vol_phi_old) < 0.1*options->eps_stop) break;
    } /* End outer loop */
    
    MERGE_SETS(p->phi, p->mask_nw, g);
    
    return t;
}


