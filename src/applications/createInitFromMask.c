
/* System headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "QSSLIB_config.h"
#include "qss_initialization3d.h"
#include "qss_initialization2d.h"
#include "qss_data_arrays.h"
#include "qss_macros.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_reinitialization3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_spatial_derivatives3d.h"


int main(int argc, char *argv[])
{
    
    QSSLIB_REAL  x_lo[3], x_hi[3];
    int  n[3], n_local[3], i, j, k, idx, idx1, dim;
    char    fname_mask[256], fname_grid[256], fname[256];	

    sprintf(fname_mask,"%s",argv[1]);
    sprintf(fname_grid,"%s",argv[2]);
    
    Grid *g = readGridFromBinaryFile(fname_grid);
    
    QSSLIB_REAL    *data_init, *mask, *data_reservoir;
    QSSLIB_REAL    center_x, center_y, center_z, radius;

    data_init = (QSSLIB_REAL *)malloc((g->num_gridpts)*DSZ);
    data_reservoir = (QSSLIB_REAL *)malloc((g->num_gridpts)*DSZ);
    
    mask = (QSSLIB_REAL *)malloc((g->num_gridpts)*DSZ);
    
    mask = readDataArrayQSS(g->grid_dims_ghostbox,fname_mask);
    
    if (g->num_dims == 3){
        radius = 0.7*(g->x_hi[0] - g->x_lo[0]);    
        center_y = 0.5*(g->x_lo[1] + g->x_hi[1]);
        center_x = (g->x_lo[0]) - radius*0.7; //-0.1
        center_z = 0.5*(g->x_lo[2] + g->x_hi[2]);;//  - radius*0.7;

        
        createSphere(data_init,center_x,center_y,center_z,radius,-1,g);
        
        /* Impose reservoir inlet */
        
        QSSLIB_REAL corner_x, corner_y, corner_z;
        QSSLIB_REAL side_length_x, side_length_y, side_length_z;
  
        corner_x = g->x_lo_ghostbox[0];
        corner_y = g->x_lo_ghostbox[1];
        corner_z = g->x_lo_ghostbox[2];
  
        side_length_x = 3*g->dx[0];
        side_length_y = g->x_hi_ghostbox[1] - g->x_lo_ghostbox[1];
        side_length_z = g->x_hi_ghostbox[2] - g->x_lo_ghostbox[2];
  
        createBox(data_reservoir, corner_x, corner_y, corner_z,
            side_length_x, side_length_y, side_length_z,
            -1, g);
        IMPOSE_MIN(data_init, data_init, data_reservoir, g);
        
    } else if (g->num_dims == 2) {
        center_y = 0.5*(g->x_lo[1] + g->x_hi[1]);
        center_x = 0.5*(g->x_lo[0] + g->x_hi[0]);
        radius = 90*g->dx[0];

        createCircle(data_init,center_x,center_y,radius,-1,g);
    }  	
    
        
    IMPOSE_MASK(data_init,mask,data_init,g);  
    
    

    sprintf(fname,"data_init");
    writeDataArrayQSS(data_init,g,fname,GZIP);
    
    destroyGrid(g);     
    free(data_init);
    free(data_reservoir);
    free(mask);


    return 0;
}

