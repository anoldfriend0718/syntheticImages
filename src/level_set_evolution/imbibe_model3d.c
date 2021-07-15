 
/* System headers */
#include <stdio.h>
#include <stdlib.h>

/* Local headers */
#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_general_util.h"
#include "qss_reinitialization3d.h"
#include "qss_grid.h"
#include "qss_macros.h"
#include "constCurvModel3d.h"
#include "connectivity.h"

/* Main driver for constant curvature level set method model */

void imbibe_model3d(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out) 
{  
    char fname[256];
    QSSLIB_REAL a0 = options->a;
    QSSLIB_REAL da = options->b * options->dc;
    int step = options->init_step + 1;
    QSSLIB_REAL t_step = 0;
    
    /* Initialize disconnected component masks */
    initializeDisconnectedMasks(p->mask_w, g);
    initializeDisconnectedMasks(p->mask_nw, g);
    
    COPY_DATA(p->mask_disconn_init, p->mask_nw, g);
    
    /* Create reservoir inlet */
    createReservoirInlet3d(p, g);
    
    /* Set amin */
    QSSLIB_REAL amin = options->amin;

    /* Compute gradients of mask */
    QSS3D_CENTRAL_GRAD_ORDER2(p->mask_x, p->mask_y, p->mask_z,
        GB_DIMS, p->mask, GB_DIMS, FB_DIMS,
        &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
    
    /* Impose reservoir inlet, if in options */
    if (options->reservoir_inlet) {
        IMPOSE_MIN(p->phi, p->phi, p->phi_extra, g);
        IMPOSE_MASK(p->phi, p->mask, p->phi, g);
    }
    
    /* Get inlet index and outlet index*/
    getMainInd_new(options, p, g);
    
    while (a0 > amin) 
    { 
        printf("\nStep %d: a = %f\n", step, a0);
        
        /* Compute variational curvature term b and advective term vel */      
        QSS3D_SET_VAR_CURV_ADV(p->mask, p->mask_x, p->mask_y, p->mask_z,
     	    p->curvature_coeff, p->external_velocity_x, p->external_velocity_y,
     	    p->external_velocity_z, &(options->b_max_over_dx), &(options->max_U_over_dx),
            &(options->b), GB_DIMS, FB_DIMS, &(g->dx[0]));
            
        /* compute variational a. variational b and vel are already computed. */
        QSS3D_SET_VAR_NORM(p->mask, p->mask_x, p->mask_y, p->mask_z,
     	    p->normal_velocity, &(a0), &(options->theta), GB_DIMS, FB_DIMS, &(g->dx[0])); 
     	
     	if(options->check_connectivity) trapComponents_mask(p, g, options);
     	
        t_step = constCurvModel3d(options, p, g, fp_out);
        
        a0 -= da;
        
        /* Impose reservoir inlet, if in options */
        if (options->reservoir_inlet) {
            IMPOSE_MIN(p->phi, p->phi, p->phi_extra, g);
            IMPOSE_MASK(p->phi, p->mask, p->phi, g);
        }
        
        /* Get inlet index and outlet index, again. Things might have moved*/
        getMainInd_new(options, p, g);
        
        if (t_step < options->tmax)
            sprintf(fname,"data_step_%d",step);
	    else
		    sprintf(fname,"data_step_%d_no_eq",step);

	    step += 1;
	    writeDataArrayQSS(p->phi,g,fname,GZIP);
        
    }
      
}


