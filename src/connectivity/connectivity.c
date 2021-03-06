#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "qss_data_arrays.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "connectivity.h"
#include "qss_macros.h"

typedef struct _union_find{
    int* parent;
    int size;
} union_find;

union_find make_sets(int size) {
    union_find result;
    result.parent = malloc(sizeof(int) * size);
    result.size = size;
    int i;
    for (i = 0; i < size; ++i) {
        result.parent[i] = size;
    }

    return result;
}

int find(union_find uf, int i) {
    if (uf.parent[i] < uf.size)
        return uf.parent[i] = find(uf, uf.parent[i]);
    return i;
}

void do_union(union_find uf, int i, int j) {
    int pi = find(uf, i);
    int pj = find(uf, j);
    if (pi == pj) {
        return;
    }
    if (pi < pj) {
        // link the smaller group to the larger one
        uf.parent[pi] = pj;
    } else if (pi > pj) {
        // link the smaller group to the larger one
        uf.parent[pj] = pi;
    } else {
        // equal rank: link arbitrarily and increase rank
        uf.parent[pj] = pi;
        ++uf.parent[pi];  
    }
}

void unionCoords2d(int x, int y, int x2, int y2, union_find component, unsigned char *input, int h, int w)
{
    int ind1 = x*h + y;
    int ind2 = x2*h + y2;
    if (y2 < h && x2 < w && input[ind1] && input[ind2] && y2 >= 0 && x2 >= 0)
        do_union(component, ind1, ind2);
}

void findConnectivity2d(QSS_DataArrays *p, Grid *g)
{
    int i, j, x, y;
    int w, h, c, c1;
    
    w = g->grid_dims_ghostbox[1];
    h = g->grid_dims_ghostbox[0];
               
    union_find component = make_sets(w*h);
        
    for (x = 0; x < w; x++)
        for (y = 0; y < h; y++){
            unionCoords2d(x, y, x+1, y, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x, y+1, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x-1, y, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x, y-1, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x+1, y+1, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x-1, y+1, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x+1, y-1, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x-1, y-1, component, p->phi_bin, h, w);
        }

    /* update connectivity array */
    for (x = 0; x < w; x++)
    {
        for (y = 0; y < h; y++)
        {
            c = x*h + y;
            if (p->phi_bin[c] == 0)
            {
                p->connectivity[c] = 0;
                continue;
            }
            p->connectivity[c] = find(component, c);

        }
    }
    
    free(component.parent);
}

void unionCoords3d(int x, int y, int z, int x2, int y2, int z2, union_find component, \
        unsigned char *input, int h, int w, int d)
{
    int ind1, ind2, nxy;
    
    nxy = h*w;
    ind1 = y + x*h + z*nxy;
    ind2 = y2 + x2*h + z2*nxy;
    
    if (y2 < h && x2 < w && z2 < d && input[ind1] && input[ind2] \
            && y2 >= 0 && x2 >= 0 && z2 >= 0)
        do_union(component, ind1, ind2);
}

void findConnectivity3d(QSS_DataArrays *p, Grid *g)
{
    int i, j, k, x, y, z;
    int w, h, d, c, c1;
    
    w = g->grid_dims_ghostbox[1];
    h = g->grid_dims_ghostbox[0];
    d = g->grid_dims_ghostbox[2];

    union_find component = make_sets(w*h*d);
    
    for (z = 0; z < d; z++){
        for (x = 0; x < w; x++){
            for (y = 0; y < h; y++){
                /* 6-connectivity */
                unionCoords3d(x, y, z, x+1, y, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y+1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y, z+1, component, p->phi_bin, h, w, d);
                
                unionCoords3d(x, y, z, x-1, y, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y-1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y, z-1, component, p->phi_bin, h, w, d);
                
                /* 18-connectivity */
                unionCoords3d(x, y, z, x-1, y-1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y-1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y+1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y+1, z, component, p->phi_bin, h, w, d);

                unionCoords3d(x, y, z, x-1, y, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y-1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y-1, z+1, component, p->phi_bin, h, w, d);
                        
                unionCoords3d(x, y, z, x+1, y, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y+1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y+1, z-1, component, p->phi_bin, h, w, d);
                
                /* 26-connectivity */
                unionCoords3d(x, y, z, x-1, y-1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y-1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y+1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y+1, z-1, component, p->phi_bin, h, w, d);

                unionCoords3d(x, y, z, x-1, y-1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y-1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y+1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y+1, z+1, component, p->phi_bin, h, w, d);
            }
        }
    }
    
    for (z = 0; z < d; z++){
        for (x = 0; x < w; x++){
            for (y = 0; y < h; y++){
                c = y + x*h + z*w*h;
                if (p->phi_bin[c] == 0)
                {
                    p->connectivity[c] = 0;
                    continue;
                }
            
                p->connectivity[c] = find(component,c);
            }
        }
    }
}

void getMainInd_new(Options *o, QSS_DataArrays *p, Grid *g)
{
    int i, j, k, idx, nw_x_ind, nw_y_ind, nw_z_ind, w_x_ind, w_y_ind, w_z_ind, ind_init;
    
    if (g->num_dims == 2) {
    /* Juanes2d case: inlet is in the center of domain */
    
    nw_x_ind = 0.5*(g->ilo_fb + g->ihi_fb);
    nw_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
    
    o->phi_nw_ind = nw_x_ind + nw_y_ind * (g->grid_dims_ghostbox[0]);
    
    w_x_ind = g->ilo_fb;
    w_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
    
    /* Copy wetting phase into scratch1, to generate a masked copy of the wetting phase */
    COPY_DATA(p->scratch1, p->phi, g);
    NEGATE_DATA(p->scratch1, g);
    IMPOSE_MASK(p->scratch1, p->mask, p->scratch1, g);
    
    /* search for a wetting point where no mask is present in the first x-slice */
    i = g->ilo_fb;
    for( j = g->jlo_fb; j <= g->jhi_fb; j++) {
        idx = i + j*(g->grid_dims_ghostbox[0]);
        if( p->scratch1[idx] < 0) {
            o->phi_w_ind = idx;
            break;    
        }
    }  
    
    } else if (g->num_dims == 3) {
    
    /* Search for inlet (or non-wetting index point) in initial x-ghost boundary */
    nw_x_ind = g->ilo_gb;
    nw_y_ind = 0.5*(g->jlo_gb + g->jhi_gb);
    nw_z_ind = 0.5*(g->klo_gb + g->khi_gb);
    
    /* Start searching near the center */
    ind_init = nw_x_ind + nw_y_ind * (g->grid_dims_ghostbox[0]) 
        + nw_z_ind * (g->grid_dims_ghostbox[0]) * (g->grid_dims_ghostbox[1]);
    
    /* Search for a point near the initial point, as that may be masked */
    for ( i = ind_init; i < g->num_gridpts; i++) {
        if (p->phi[i] < 0) {
            o->phi_nw_ind = i;
            break;
        }
    }
    
    /* Copy wetting phase into scratch1, to generate a masked copy of the wetting phase */
    COPY_DATA(p->scratch1, p->phi, g);
    NEGATE_DATA(p->scratch1, g);
    IMPOSE_MASK(p->scratch1, p->mask, p->scratch1, g);
    
    /* Search for outlet (or wetting index point) in last x-ghost boundary */
    w_x_ind = g->ihi_gb;
    w_y_ind = 0.5*(g->jlo_gb + g->jhi_gb);
    w_z_ind = 0.5*(g->klo_gb + g->khi_gb);
    
    /* Start searching near the center */
    ind_init = w_x_ind + w_y_ind * (g->grid_dims_ghostbox[0]) 
        + w_z_ind * (g->grid_dims_ghostbox[0]) * (g->grid_dims_ghostbox[1]);
        
    /* Search for a point near the initial point, as that may be masked */    
    for ( i = ind_init; i < g->num_gridpts; i++) {
        if (p->scratch1[i] < 0) {
            o->phi_w_ind = i;
            break;
        }
    }   
    }
        
}

void getMainInd(Options *o, QSS_DataArrays *p, Grid *g)
{
    int i;
    for (i = 0; i < g->num_gridpts; i++)
        if (p->phi[i] < 0) {
            o->phi_nw_ind = i;
            break;
        }
        
        
    /* Copy wetting phase into scratch1, to generate a masked copy of the wetting phase */
    COPY_DATA(p->scratch1, p->phi, g);
    NEGATE_DATA(p->scratch1, g);
    IMPOSE_MASK(p->scratch1, p->mask, p->scratch1, g);
    
    for (i = 0; i < g->num_gridpts; i++)
        if (p->scratch1[i] < 0) {
            o->phi_w_ind = i;
            break;
        }

}