/*
  Copyright (c) 2016 Lara Raad <lara.raad@cmla.ens-cachan.fr>,
                     Bruno Galerne <bruno.galerne@parisdescartes.fr>

  Quilting is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "blending.h"
#include "io_png.h"

#define min(a,b) (a<=b)?(a):(b)
#define max(a,b) (a>b)?(a):(b)
#define abs(a) (a<0)?-(a):(a)

/**
*   This method blends new block in the output image with the corresponding 
*   "new edges" defined by the boudary mask. The new block will be composed of
*   the synthesized part in the top and/or left part of the minimal error 
*   cut and of the new block from the source image in the rest, with a mix of
*   both in the transition zone if the smoothed option is set to 1.
*   @author Lara Raad
*   @param boundary_mask    boundary_mask is the mask used to redefine 
*                           the ragged edges of the patches
*   @param out_im           The synthesized texture
*   @param src_im           The input texture sample
*   @param patch_src        Corner position of the patch to copy from 
*                           the input texture sample
*   @param temp             Corner position of the location were to add the 
*                           new patch in the output image
*   @param ps               The size of the patches
*   @param ps               The overlap size
*   @date 24/06/2014
*/
int blend_patch(Image boundary_mask, Image *out_im, Image src_im, 
                Corner patch_src, Corner temp, int ps, int os)
{
    int pos_o, pos_s, pos_m, ii, jj, kk, inc_o, inc_s, x, y;

    int smooth = 1;

    x = temp.x;
    y = temp.y;
    inc_s = src_im.w*src_im.h;
    inc_o = out_im->w*out_im->h;
    
    if (smooth && os!= 1)
        smooth_mask(&boundary_mask,ps,os,temp);
    
    for(ii=0; ii<ps; ii++)
    {
        for(jj=0; jj<ps; jj++)
        {
            pos_m = ii*boundary_mask.w + jj;
            for(kk=0; kk<3; kk++)
            {
                pos_o = (ii+x)*out_im->w + (jj+y) + kk*inc_o;
                pos_s = (ii+patch_src.x)*src_im.w + (jj+patch_src.y) + kk*inc_s;
                out_im->img[pos_o] = 
                        out_im->img[pos_o]*(boundary_mask.img[pos_m]) +
                        src_im.img[pos_s]*(1-(boundary_mask.img[pos_m]));
            }
        }
    }

    return 0;
}

/**
*   This method smoothes the boudary mask using a Gaussian kernel. 
*   The transition between the old part and the new part is smoother.
*   @author Lara Raad
*   @param boundary_mask     boundary_mask is the mask used to redefine 
*                            the ragged edges of the patches
*   @param ps                The patch size
*   @param os                The overlap size
*   @param temp              Corner position of the location were to add 
*                            the new patch in the output image
*   @date 07/07/2014
*/
int smooth_mask(Image *boundary_mask, int ps, int os, Corner temp)
{
    int ii;
    int kernel_size = 3;
    float kernel[kernel_size];
    int center = kernel_size - 1;
    float var = 0.8;//gaussian kernel's variance
    
    //* Create Gaussian kernel *//
    for(ii=0; ii<kernel_size; ii++)
    {
        kernel[ii] = expf(-(ii*ii)/(2*var));
    }

    //* Filtering the boundary mask*//
    int X1[2]; int Y1[2]; int X2[2]; int Y2[2];
    if (temp.x>0 && temp.y>0)//L-shape overlap
    {
        X1[0] = 0; X1[1] = os; Y1[0] = 0; Y1[1] = ps;
        X2[0] = os; X2[1] = ps; Y2[0] = 0; Y2[1] = os;
        filter_mask(X1,Y1,X2,Y2,ps,boundary_mask,kernel,center);
    }
    else if (temp.x==0 && temp.y>0)//vertival overlap
    {
        X1[0] = 0; X1[1] = os; Y1[0] = 0; Y1[1] = os;
        X2[0] = os; X2[1] = ps; Y2[0] = 0; Y2[1] = os;
        filter_mask(X1,Y1,X2,Y2,ps,boundary_mask,kernel,center);
    }
    else if (temp.x>0 && temp.y==0)//horizontal overlap
    {
        X1[0] = 0; X1[1] = os; Y1[0] = 0; Y1[1] = os;
        X2[0] = 0; X2[1] = os; Y2[0] = os; Y2[1] = ps;
        filter_mask(X1,Y1,X2,Y2,ps,boundary_mask,kernel,center);
    }
    
    return 0;
}


/**
*   This method partially filters an image. We pass X1, X2, X3 and X4
*   to mark out the area we wish to filter. It is done in two steps to 
*   to contemplate the L-shaped case. With X1 and Y1 we filter a first area.
*   With X2 ans Y2 we filter a second area. X1 and X2 have the coordinates
*   in x that mark out the first and second area vertically. And Y1 and Y2
*   the coordinates in y that mark out the first and second area horizontally.
*   @author Lara Raad
*   @param X1                2-d vector with the x-coordinates that mark out 
*                            horizontally the first area
*   @param Y1                2-d vector with the y-coordinates that mark out 
*                            horizontally the first area
*   @param X2                2-d vector with the x-coordinates that mark out 
*                            horizontally the second area
*   @param Y2                2-d vector with the y-coordinates that mark out 
*                            horizontally the second area
*   @param ps                The patch size
*   @param boundary_mask     Mask to filter
*   @param kernel            The kernel used to filter the mask
*   @param center            The kernel's center
*   @date 07/07/2014
*/
int filter_mask(int X1[],int Y1[], int X2[], int Y2[], int ps, 
                Image *boundary_mask, float kernel[], int center)
{
    int left, right, top, bottom, ii, jj, kk, ll, pos, pos_aux, pos_blend;
    float tot_weight;
    Image aux_mask;
    
    //* Initialize auxiliar mask *//
    initialize_image(&aux_mask, ps, ps, 1);
    
    for(ii=0; ii<X1[1]; ii++)
    {
        for(jj=0; jj<Y1[1]; jj++)
        {
            left = max(0, jj-center); right = min(Y1[1]-1, jj+center);
            top = max(0, ii-center); bottom = min(X1[1]-1, ii+center);
            pos_aux = ii*aux_mask.w + jj;
            tot_weight = 0;
            for(kk=top; kk<=bottom; kk++)
            {
                for(ll=left; ll<=right; ll++)
                {
                    pos_blend = kk*boundary_mask->w + ll;
                    aux_mask.img[pos_aux] +=  
                        boundary_mask->img[pos_blend]*
                        kernel[abs(ii-kk)]*kernel[abs(jj-ll)];
                    tot_weight += kernel[abs(ii-kk)]*kernel[abs(jj-ll)];
                }
            }
            aux_mask.img[pos_aux] /= tot_weight;
        }
    }

    for(ii=X2[0]; ii<X2[1]; ii++)
    {
        for(jj=Y2[0]; jj<Y2[1]; jj++)
        {
            left = max(Y2[0], jj-center); right = min(Y2[1]-1, jj+center);
            top = max(X2[0], ii-center); bottom = min(X2[1]-1, ii+center);
            pos_aux = ii*aux_mask.w + jj;
            tot_weight = 0;
            for(kk=top; kk<=bottom; kk++)
            {
                for(ll=left; ll<=right; ll++)
                {
                    pos_blend = kk*boundary_mask->w + ll;
                    aux_mask.img[pos_aux] +=  
                        boundary_mask->img[pos_blend]*
                        kernel[abs(ii-kk)]*kernel[abs(jj-ll)];
                    tot_weight += kernel[abs(ii-kk)]*kernel[abs(jj-ll)];
                }
            }
            aux_mask.img[pos_aux] /= tot_weight;
        }
    }

    //* Updating values of smoothed mask *//
    for(ii=0; ii<X1[1]; ii++)
    {
        for(jj=0; jj<Y1[1]; jj++)
        {
            pos = ii*aux_mask.w + jj;
            boundary_mask->img[pos] = aux_mask.img[pos];
        }
    }
    for(ii=X2[0]; ii<X2[1]; ii++)
    {
        for(jj=Y2[0]; jj<Y2[1]; jj++)
        {
            pos = ii*aux_mask.w + jj;
            boundary_mask->img[pos] = aux_mask.img[pos];
        }
    }

    //* free memory *//
    free(aux_mask.img);
    
    return 0;
}
