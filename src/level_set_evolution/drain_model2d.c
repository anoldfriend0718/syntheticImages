 
/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/* Local headers */
#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_spatial_derivatives2d.h"
#include "qss_tvd_runge_kutta2d.h"
#include "qss_data_arrays.h"
#include "qss_util2d.h"
#include "qss_reinitialization2d.h"
#include "qss_general_util.h"
#include "qss_grid.h"
#include "constCurvModel2d.h"
#include "qss_macros.h"
#include "connectivity.h"

/* Main driver for constant curvature level set method model */

void drain_model2d(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out) 
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

    /* Set amax so that at least 80 steps are possible */
    QSSLIB_REAL amax = options->amax;

    if (options->mask_reinit)
        qss_reinitialize_mask2d(p, g, options);
        
    /* Compute gradients of mask */
    QSS2D_CENTRAL_GRAD_ORDER2(p->mask_x, p->mask_y,
        GB_DIMS_2D, p->mask, GB_DIMS_2D, FB_DIMS_2D,
        &((g->dx)[0]),&((g->dx)[1]));
    
    /* Get inlet index and outlet index */
    getMainInd_new(options, p, g);
    
    while (a0 < amax) 
    { 
        printf("\nStep %d: a = %f\n", step, a0);
        
        /* Compute variational curvature term b and advective term vel */      
        QSS2D_SET_VAR_CURV_ADV(p->mask, p->mask_x, p->mask_y,
     	    p->curvature_coeff, p->external_velocity_x, p->external_velocity_y,
     	    &(options->b_max_over_dx), &(options->max_U_over_dx), &(options->b),
            GB_DIMS_2D, FB_DIMS_2D, &(g->dx[0]));
            
        /* compute variational a. variational b and vel are already computed. */
        QSS2D_SET_VAR_NORM(p->mask, p->mask_x, p->mask_y,
     	    p->normal_velocity, &(a0), &(options->theta), GB_DIMS_2D, FB_DIMS_2D, &(g->dx[0])); 
     	  
     	if(options->check_connectivity) trapComponents_mask(p, g, options);

        t_step = constCurvModel2d(options, p, g, fp_out, a0);
        
        a0 += da;
        
	    if (t_step < options->tmax)
        	sprintf(fname,"data_step_%d",step);
	    else
		    sprintf(fname,"data_step_%d_no_eq",step);

	    step += 1;
        writeDataArrayQSS(p->phi,g,fname,GZIP);
        
    }

      
}


