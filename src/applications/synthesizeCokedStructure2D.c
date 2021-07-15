#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "QSSLIB_config.h"
#include "qss_util3d.h"
#include "qss_grid.h"
#include "drain_top.h"
#include "qss_data_arrays.h"
#include "drain_model3d.h"
#include "drain_model2d.h"
#include "qss_macros.h"
#include "qss_options.h"
#include "qss_general_util.h"

// /* number of grid cells for "film" */
// #define n_cells 1

// /* fraction of indices close to mask turned back to pore space */
// #define frac_close 0.5

void grow(unsigned char *data, int ind, Grid *g, QSSLIB_REAL *mask) {
    int w = g->grid_dims_ghostbox[1];
    int h = g->grid_dims_ghostbox[0];
    int d = g->grid_dims_ghostbox[2];
    
    int i, neighbor_idx;
    int nxy = h*w;
    /* Get neighbors of idx */
    /* ind -1, ind + 1, ind + h, ind - h, ind + nxy, ind - nxy */
    int neighbors[6] = {ind - 1, ind + 1, ind + h, ind - h, ind + nxy, ind - nxy};

    for (i = 0; i < 6; i++) {
        neighbor_idx = neighbors[i];
        if ((neighbor_idx >= 0) && (neighbor_idx < g->num_gridpts)) {
            if (mask[neighbor_idx] < 0) /* if neighbor is in pore space */
                data[neighbor_idx] = 0; /* change neighbor to mask */
        }    
    }
}


int main(int argc, char* argv[])
{
  //input parameter 
  int n_cells=atoi(argv[1]);
  float frac_close=atof(argv[2]);
  printf("growing parameter: layer cell: %d, fraction: %f \n",n_cells,frac_close);
  int grow_num=atoi(argv[3]);
  printf("growing number: %d\n",grow_num);

    /* structure containing all arrays */
  QSS_DataArrays *p; 
  
  int *indices;
  /* grid structure */
  Grid *g;   
  Options *o;
  
  o = createOptionsDefault();
  
  int     n1[3], n2[3], i, temp, rand_idx, num_part_pixels;
  char    fname[256];
  FILE    *fp_out; 
  unsigned char *data;  
  p = allocateQSSDataArrays();
  
  sprintf(fname, "grid.gz");
  g = readGridFromBinaryFile(fname);
  
  allocateMemoryForQSSDataArrays(p,g);
  data = (unsigned char *)malloc(g->num_gridpts*sizeof(unsigned char));
  indices = (int *)malloc(g->num_gridpts*sizeof(int));
  
  /* Initialize indices array to -1 */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  sprintf(fname, "mask.gz");
  p->mask = readDataArrayQSS(n2,fname);
  
  /* Add a uniform film next to the mask everywhere, and store in p->phi */
  COPY_DATA(p->phi, p->mask, g);
  for (i = 0; i < g->num_gridpts; i++)
        p->phi[i] += n_cells*(g->dx[0]); //set the N layer of the pore-space near the grain surface as grain (positive values)
  
  reinitialize2d_subcell_fix_qss(p, g, o);


  /* convert cells in pore space, without the film to 1, rest is 0 */
  /* This is uniform grain-attaching, which is not typically observed */
  IMPOSE_UCHAR(data, p->phi, g);
  sprintf(fname,"uniformFilmUchar");
  writeDataArrayUchar(data, g, fname, 0);
  printf("output synthetic image: %s \n",fname);
  
  int idx = 0;
  /* convert cells in pore space (total) to 1, rest is 0 */
  IMPOSE_UCHAR(data, p->mask, g);
  sprintf(fname,"originalPoreSpace");
  writeDataArrayUchar(data, g, fname, 0);
  
  /* Get indices close to mask */
  //p->mask[i] < 0 is the voxel in the pore space
  //p->mask[i] > -n_cells*(g->dx[0]) is the voxel within the nearest-grain-space N layer pore space
  for (i = 0; i < g->num_gridpts; i++){
    if((p->mask[i] < 0) && (p->mask[i] >= -n_cells*(g->dx[0]))){
        indices[idx] = i;
        idx++;
    }
  }
  /*  number of indices near to mask (rest of indices array is -1)*/
  int n_indices_near = idx;
  int num_turn_back = frac_close * n_indices_near;
  for (i = n_indices_near-1; i >= 0; --i){
    //generate a random number [0, n-1]
    int j = rand() % (i+1);
    //swap the last element with element at random index
    temp = indices[i];
    indices[i] = indices[j];
    indices[j] = temp;
  }
  /* Now the indices array should be jumbled, so just pick the first num_turn_back indices
    to turn to mask */
    /* This is grain-attaching mode, initial grain-attaching, only convert a portion*/ 
  for (i = 0; i < num_turn_back; i++) {
    idx = indices[i];
    data[idx] = 0; //set the solid 
  }
  sprintf(fname,"randFilmTurnBack");
  writeDataArrayUchar(data, g, fname, 0);
  printf("output synthetic image: %s \n",fname);
  
  /* Grow particle, first pass */
  /* This is grain-attaching with coarsening processing*/
  for (i = 0; i < num_turn_back; i++) {
    idx = indices[i];
    //grow all the 4(2D)/6(3D) neighboring voxels of the previous randomly picked voxels into the solid
    grow(data, idx, g, p->mask); 
  }
  sprintf(fname,"randFilmTurnBack_grow_1");
  writeDataArrayUchar(data, g, fname, 0);
  printf("output synthetic image: %s \n",fname);
  

  /* grow again */ 
  for(int grow_index=2; grow_index<=grow_num; grow_index++)
  {
     //reset the indices 
    SET_DATA_TO_CONSTANT(indices, g, -1);
    /* get indices of particles to grow */
     /* Get indices close to mask */
    idx = 0;
    for (i = 0; i < g->num_gridpts; i++){
      //p->mask[i] < 0 is the voxel in the original pore space (mask is fixed during simulation)
      //data[i] == 0 is the solid voxel 
      // "and" means the new growed solid voxel in the original pore space 
      if( (p->mask[i] < 0) && (data[i] == 0) ){
          indices[idx] = i;
          idx++;
      }
    }
    num_part_pixels = idx;
    
    /* Grow particle */
    /* This is grain-attaching with coarsening processing*/

    for (i = 0; i < num_part_pixels; i++) {
      idx = indices[i];
      //grow all the 4(2D)/6(3D) neighboring voxels of all the new solid voxels into the solid
      grow(data, idx, g, p->mask);
    }
    
    char fname[30]="randFilmTurnBack_grow_";
    char index[5];
    sprintf(index, "%d", grow_index);
    strcat(fname, index);
    writeDataArrayUchar(data, g, fname, 0);

    printf("output synthetic image: %s \n",fname);
  }

  destroyQSSDataArrays(p);
  destroyGrid(g); 
  free(data);
  free(indices);
}
