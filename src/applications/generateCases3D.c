#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

/* number of grid cells for "film" */
#define N_CELLS 10

/* fraction of indices close to mask turned back to pore space */
#define FRAC_CLOSE 0.5

/* fraction of indices far from mask turned back to pore space */
#define FRAC_FAR 0.0001

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


int main()
{

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
        p->phi[i] += N_CELLS*(g->dx[0]);
  
  reinitialize3d_subcell_fix_qss(p, g, o);
  /*
  sprintf(fname,"uniformFilm");
  writeDataArrayQSS(p->phi,g,fname,GZIP);
  */
  
  /* convert cells in pore space, without the film to 1, rest is 0 */
  /* This is uniform grain-attaching, which is not typically observed */
  IMPOSE_UCHAR(data, p->phi, g);
  
  sprintf(fname,"uniformFilmUchar");
  writeDataArrayUchar(data, g, fname, 0);
  
  int idx = 0;
  /* convert cells in pore space (total) to 1, rest is 0 */
  IMPOSE_UCHAR(data, p->mask, g);
  sprintf(fname,"originalPoreSpace");
  writeDataArrayUchar(data, g, fname, 0);
  
  /* Get indices close to mask */
  for (i = 0; i < g->num_gridpts; i++){
    if((p->mask[i] < 0) && (p->mask[i] > -N_CELLS*(g->dx[0]))){
        indices[idx] = i;
        idx++;
    }
  }
  
  /*  number of indices near to mask (rest of indices array is -1)*/
  int n_indices_near = idx;
  
  int num_turn_back = FRAC_CLOSE * n_indices_near;
  
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
    data[idx] = 0;
  }
  
  sprintf(fname,"randFilmTurnBack");
  writeDataArrayUchar(data, g, fname, 0);
  
  /* Grow particle, first pass */
  /* This is grain-attaching with coarsening processing*/
  
  for (i = 0; i < num_turn_back; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_1");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy starts here */  
/* tertiary growth */

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
    /* This is grain-attaching with coarsening processing*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_2");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  


/* copy starts here */  
/* 4th growth */

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, 3rd pass */
    /* This is grain-attaching with coarsening processing*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_3");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 5th growth */

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, 4th pass */
    /* This is grain-attaching with coarsening processing*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_4");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 6th growth */
  /* This is grain-attaching with coarsening processing*/


  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_5");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 7th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_6");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  


/* copy starts here */  
/* 8th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_7");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 9th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_8");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 10th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_9");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 11th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_10");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 12th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_11");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 13th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_12");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 14th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_13");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  

/* copy starts here */  
/* 15th growth */
  /* This is grain-attaching with coarsening processing*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices close to mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if( (p->mask[i] < 0) && (data[i] == 0) ){
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randFilmTurnBack_grow_14");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here */  


  /* Do same thing, but for indices far from mask */
  SET_DATA_TO_CONSTANT(indices, g, -1);

  IMPOSE_UCHAR(data, p->mask, g);
  
  /* get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -N_CELLS*(g->dx[0])) {
        indices[idx] = i;
        idx++;
    }
  }
  
  /* number of indices far (rest of indices array is -1)*/
  int n_indices_far = idx;
  
  int num_turn_back_far = FRAC_FAR * n_indices_far;
  
 for (i = n_indices_far-1; i >= 0; --i){
    //generate a random number [0, n-1]
    int j = rand() % (i+1);

    //swap the last element with element at random index
    temp = indices[i];
    indices[i] = indices[j];
    indices[j] = temp;
}
  
  /* Now the indices array should be jumbled, so just pick the first num_turn_back indices
    to turn to mask */
  /* This is the dispersed pore-filling habit*/
    
  for (i = 0; i < num_turn_back_far; i++) {
    idx = indices[i];
    data[idx] = 0;
  }
  
  sprintf(fname,"randPoreTurnBack");
  writeDataArrayUchar(data, g, fname, 0);
  
  
  /* Grow particle, first pass */
  /* This is the coarse pore filling pore habit*/
  
  for (i = 0; i < num_turn_back_far; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_1");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
    /* This is the coarse pore filling pore habit*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_2");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/  

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, 3rd pass */
    /* This is the coarse pore filling pore habit*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_3");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, 4th pass */
    /* This is the coarse pore filling pore habit*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_4");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, 5th pass */
    /* This is the coarse pore filling pore habit*/

  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_5");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
    /* This is the coarse pore filling pore habit*/

  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_6");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_7");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_8");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_9");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_10");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_11");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_12");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_13");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

/* copy starts here*/

  /* grow again */
  SET_DATA_TO_CONSTANT(indices, g, -1);
  
  /* get indices of particles to grow */
    /* Get indices far from mask */
  idx = 0;
  for (i = 0; i < g->num_gridpts; i++){
    if(p->mask[i] < -(g->dx[0]) && (data[i] == 0)) {
        indices[idx] = i;
        idx++;
    }
  }
  num_part_pixels = idx;
  
  /* Grow particle, second pass */
  for (i = 0; i < num_part_pixels; i++) {
    idx = indices[i];
    grow(data, idx, g, p->mask);
  }
  
  sprintf(fname,"randPoreTurnBack_grow_14");
  writeDataArrayUchar(data, g, fname, 0);
  
/* copy ends here*/ 

  destroyQSSDataArrays(p);
  destroyGrid(g); 
  free(data);
  free(indices);
}
