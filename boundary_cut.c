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

#include "boundary_cut.h"
#include "random_number.h"

#define min(a,b) (a<=b)?(a):(b)
#define max(a,b) (a>b)?(a):(b)
#define abs(a) (a<0)?-(a):(a)
#define true 1 
#define false 0

/**
*   This method computes the cumulative veritcal error in the overlap region 
*   and the corresponding patches to these errors.
*   @author Lara Raad
*   @param src_im       The input texture sample
*   @param out_im       The synthesized texture
*   @param E            E is the cumulative error in the overlap
*   @param paths        paths saves the minimal path to each pixel
*                       of the overlap zone
*   @param os           Overlap size 
*   @param ps           Patch size
*   @param x            X position of the location were to add
*                       the new patch in the output image
*   @param y            Y position of the location were to add 
*                       the new patch in the output image
*   @param patch_src    Corner position of the patch to copy 
*                       from the input texture sample
*   @date 22/06/2014
*/
int cumulative_vertical_error(Image *out_im, Image src_im, Image *E,
                 Image *paths, int ps, int os, int x, int y, Corner patch_src)
{
    int ii,jj, kk, ll, mm, pos_o, pos_e, pos_s, posE_x, pos_P,
         aux1, aux2, aux3, m1, m2, inc_o, inc_s;
    Image e;// e is the error along the overlap region
    
    initialize_image(&e, ps, os,1);
    inc_o = out_im->h*out_im->w; // number of pixels in 1 channel of output
    inc_s = src_im.h*src_im.w; // number of pixels in 1 channel of source
    
    kk=0; ll=0; mm=0;

    for (ii=0; ii<ps; ii++)
    {
        kk = ps-1-ii;
        ll = ps-1-(ii-1);
        mm = ps-2-(ii-1);
        
        for (jj=0; jj<os; jj++)
        {      
            pos_o = (kk + x)*out_im->w + (y + jj);
            pos_e = kk*e.w + jj;
            pos_s = (patch_src.x + kk)*src_im.w +(patch_src.y +jj);
            posE_x = ll*e.w;
            pos_P = (mm)*paths->w + jj;

            assert(pos_o<(int)(out_im->h*out_im->w));
            assert(pos_e<(int)(e.h*e.w));


            e.img[pos_e] =  (out_im->img[pos_o] - src_im.img[pos_s]) * 
                (out_im->img[pos_o] - src_im.img[pos_s]) 
                +
                (out_im->img[pos_o+inc_o] - src_im.img[pos_s+inc_s]) *
                (out_im->img[pos_o+inc_o] - src_im.img[pos_s+inc_s]) 
                +
                (out_im->img[pos_o+2*inc_o] - src_im.img[pos_s+2*inc_s]) *
                (out_im->img[pos_o+2*inc_o] - src_im.img[pos_s+2*inc_s]);

            //build cumulative error image E 
            //in the first row of the vertical overlap zone the error 
            // image e and the cumulative error match
            if (ii==0)
            {
                E->img[pos_e] = e.img[pos_e];                
            }
            if (ii > 0)
            {
                assert(pos_s<(int)(src_im.h*src_im.w));
                assert(pos_P<(int)(paths->w*paths->h));
                assert(posE_x<=(int)((ps-1)*e.w));
                assert(posE_x>=0);
                assert(pos_P>=0);
                aux2 = E->img[posE_x + jj];
                if (jj==0)
                {
                    //left edge of the vertical overlap zone
                    aux3 = E->img[posE_x + (jj+1)];
                    m1 = min(aux2,aux3);                  
                    E->img[pos_e] = e.img[pos_e] + m1;
                    
                    if ( m1 == aux2 )
                        paths->img[pos_P] = 0;
                    else if ( m1 == aux3 )
                        paths->img[pos_P] = 1;   
                }
                else if (jj==(os-1))
                {
                    //right edge of the vertical overlap zone
                    aux1 = E->img[posE_x + (jj-1)];
                    m1 = min(aux1,aux2); 
                    E->img[pos_e] = e.img[pos_e] + m1;

                    if ( m1 == aux1 )
                        paths->img[pos_P] = -1;
                    else if ( m1 == aux2 )
                        paths->img[pos_P] = 0;
                }
                else
                {
                    //rest of pixels of the vertical overlap zone
                    aux1 = E->img[posE_x+ (jj-1)];
                    aux3 = E->img[posE_x + (jj+1)];
                    m1 = min(aux1,aux2);
                    m2 = min(m1,aux3);
                    E->img[pos_e] = e.img[pos_e] + m2;
                    
                    if ( m2 == aux1 )
                        paths->img[pos_P] = -1;
                    else if ( m2 == aux2 )
                        paths->img[pos_P] = 0;
                    else if ( m2 == aux3 )
                        paths->img[pos_P] = 1;
                }
                
            }            
        }
    }
    
    // release memory
    free(e.img);
    
    return 0;
}

/**
*   This method computes the cumulative horizontal error in the overlap region
*   and the corresponding patches to these errors.
*   @author Lara Raad
*   @param src_im       The input texture sample
*   @param out_im       The synthesized texture
*   @param E            E is the cumulative error in the overlap
*   @param paths        paths saves the minimal path to each pixel of the
*                       overlap zone
*   @param os           Overlap size 
*   @param ps           Patch size
*   @param x            X position of the location were to add
*                       the new patch in the output image
*   @param y            Y position of the location were to add
*                       the new patch in the output image
*   @param patch_src    Corner position of the patch to copy from
*                       the input texture sample
*   @date 22/06/2014
*/

int cumulative_horizontal_error(Image *out_im, Image src_im, Image *E, 
                Image *paths, int ps, int os, int x, int y, Corner patch_src)
{
    int ii, jj, kk, ll, mm, pos_o, pos_e, pos_s, 
        posE_y, pos_P, aux1, aux2, aux3, m1, m2;
    Image e;

    initialize_image(&e, os, ps,1);
    int inc_o = out_im->h*out_im->w;// number of pixels in 1 channel of output
    int inc_s = src_im.h*src_im.w;// number of pixels in 1 channel of source

    kk=0; ll=0; mm=0;

    for (jj=0; jj<ps; jj++)
    {
        kk = ps-1-jj;
        ll = ps - 1 - (jj-1);
        mm = ps - 2 - (jj-1);
        
        for (ii=0; ii<os; ii++)
        {
            pos_o = (x + ii)*out_im->w + (y + kk);
            pos_e = ii*e.w + kk;
            pos_s = (patch_src.x + ii)*src_im.w +(patch_src.y +kk);
            posE_y = ll;
            pos_P = ii*paths->w + mm;

            assert(pos_o<(int)(out_im->h*out_im->w));
            assert(pos_e<(int)(e.h*e.w));
            
            
            e.img[pos_e] = 
                (out_im->img[pos_o] - src_im.img[pos_s]) *
                (out_im->img[pos_o] - src_im.img[pos_s]) 
                +
                (out_im->img[pos_o+inc_o] - src_im.img[pos_s+inc_s]) *
                (out_im->img[pos_o+inc_o] - src_im.img[pos_s+inc_s]) 
                +
                (out_im->img[pos_o+2*inc_o] - src_im.img[pos_s+2*inc_s]) *
                (out_im->img[pos_o+2*inc_o] - src_im.img[pos_s+2*inc_s]);

            // build cumulative error image E 
            // in the first colmun of the horizontal overlap zone the error
            // image e and the cumulative error match
            if (jj==0)
                E->img[pos_e] = e.img[pos_e];
            if (jj > 0)
            {
                assert(pos_s<(int)(src_im.h*src_im.w));
                assert(posE_y<(int)ps);
                assert(pos_P<(int)(paths->w*paths->h));
                assert(posE_y>=0);
                assert(pos_P>=0);
                aux2 = E->img[(ii)*E->w + posE_y];
                if (ii==0)
                {
                    //considering the top edge of the horizontal overlap zone
                    aux3 = E->img[(ii+1)*E->w + posE_y];
                    m1 = min(aux2,aux3);                  
                    E->img[pos_e] = e.img[pos_e] + m1;
                    
                    if ( m1 == aux2 )
                        paths->img[pos_P] = 0;
                    else if ( m1 == aux3 )
                        paths->img[pos_P] = 1;
                }
                else if (ii==(os-1))
                {
                    //considering the bottom edge of the horizontal overlap zone
                    aux1 = E->img[(ii-1)*E->w + posE_y];
                    m1 = min(aux1,aux2); 
                    E->img[pos_e] = e.img[pos_e] + m1;

                    if ( m1 == aux1 )
                        paths->img[pos_P] = -1;
                    else if ( m1 == aux2 )
                        paths->img[pos_P] = 0;
                }
                else
                {
                    //considering the rest of pixels of the horizontal overlap
                    aux1 = E->img[(ii-1)*E->w + posE_y];
                    aux3 = E->img[(ii+1)*E->w + posE_y];
                    m1 = min(aux1,aux2);
                    m2 = min(m1,aux3);
                    E->img[pos_e] = e.img[pos_e] + m2;
                    
                    if ( m2 == aux1 )
                        paths->img[pos_P] = -1;
                    else if ( m2 == aux2 )
                        paths->img[pos_P] = 0;
                    else if ( m2 == aux3 )
                        paths->img[pos_P] = 1;
                }
            }
        }
    }
    
    // Release memory
    free(e.img);
    
    return 0;
}

/**
*   This method creates the boundary mask for the vertical overlap. 
*   @author Lara Raad
*   @param E_V             Cumulative error across the vertical overlap area
*   @param paths_V         Optimal paths to reach the positions of the last row
*   @param ps              Patch size
*   @param os              Overlap size
*   @param boundary_mask   Mask used to redefine the ragged edges of the patches
*                          in the blending step
*   @date 24/06/2014
*/
int create_vertical_mask(Image *E_V, int ps, int os, Image *boundary_mask, 
                         Image paths_V)
{
    int ymin, ii, jj, counter, num_rand;
    int ys[os];
    float minE;
    ymin = 0;
    assert((0)*E_V->w + ymin<E_V->w*E_V->h);
    minE = E_V->img[(0)*E_V->w + ymin];
    for (jj=1; jj<os; jj++)
    {
        assert((0)*E_V->w + jj < E_V->w*E_V->h);
        if (E_V->img[(0)*E_V->w + jj] < minE)
        {
            minE = E_V->img[(0)*E_V->w + jj];
        }
    }
    
    counter = 0;
    for (jj=0; jj<os; jj++)
    {
        if (E_V->img[(0)*E_V->w + jj]==minE)
        {
            ys[counter] = jj;
            counter++; 
        }
    }

    assert(counter>0);
    assert(counter<=os);
    num_rand = random_number(NULL) % counter;
    assert(num_rand < counter);
    assert(num_rand < os);
    ymin = ys[num_rand];

    //* initialize vertical mask *//
    for (jj=0; jj<=ymin; jj++)
    {
        assert((0)*boundary_mask->w + jj<boundary_mask->w*boundary_mask->h);
        boundary_mask->img[(0)*boundary_mask->w + jj] = 1;
    }
    
    //* complete vertical mask *//
    for (ii=1; ii<ps; ii++)
    {
        assert((ii-1)*(int)paths_V.w + ymin < (int)(paths_V.w*paths_V.h));
        ymin = ymin + paths_V.img[(ii-1)*(int)paths_V.w + ymin];
        for (jj=0; jj<=ymin; jj++)
        {
            assert(ii*boundary_mask->w + jj < 
                   boundary_mask->h*boundary_mask->w);
            boundary_mask->img[ii*boundary_mask->w + jj] = 1;
        }
    }
    
    return 0;
}

/**
*   This method creates the boundary mask for the horizontal overlap.
*   @author Lara Raad
*   @param E_H              The cumulative error across the horizontal 
*                           overlap arear
*   @param paths_H          The optimal paths to reach the positions 
*                           of the last column
*   @param ps               Patch size
*   @param os               Overlap size
*   @param boundary_mask    Mask used to redefine the ragged edges of 
*                           the patches in the blending step
*   @date 24/06/2014
*/
int create_horizontal_mask(Image *E_H, int ps, int os, Image *boundary_mask, 
                           Image paths_H)
{
    int xmin, ii, jj, counter, num_rand;
    int xs[os];
    float minE;
    xmin = 0;
    assert(xmin*E_H->w + 0 < E_H->w*E_H->h);
    minE = E_H->img[xmin*E_H->w + 0];
    for (ii=1; ii<os; ii++)
    {
        assert(ii*E_H->w + 0 < E_H->w*E_H->h);
        if (E_H->img[ii*E_H->w + 0] < minE)
        {
            minE = E_H->img[ii*E_H->w + 0];

        }
    }
    
    counter = 0;
    for (ii=0; ii<os; ii++)
    {
        if (E_H->img[ii*E_H->w + 0] == minE)
        {
            xs[counter] = ii;
            counter++; 
        }
    }
    
    assert(counter>0);
    assert(counter<=os);
    num_rand = random_number(NULL) % counter;
    assert(num_rand<counter);
    assert(num_rand<os);
    xmin = xs[num_rand];

    //* initialize horizontal mask *//
    for (ii=0; ii<=xmin; ii++)
    {
        assert(ii*boundary_mask->w + 0 < boundary_mask->w*boundary_mask->h);
        boundary_mask->img[ii*boundary_mask->w + 0] = 1;
    }
    
    //* complete horizontal boundary mask *//    
    for (jj=1; jj<=ps-1; jj++)
    {
        assert(xmin*(int)paths_H.w + (jj-1) < (int)(paths_H.w*paths_H.h));
        xmin = xmin + paths_H.img[xmin*(int)paths_H.w + (jj-1)];
        for (ii=0; ii<=xmin; ii++)
        {
            assert(ii*boundary_mask->w + jj < 
                   boundary_mask->h*boundary_mask->w);
            boundary_mask->img[ii*boundary_mask->w + jj] = 1;
        }
    }
    
    return 0;
}

/**
* This method creates the boundary mask for the L-shape overlap. 
* @author Lara Raad
* @param E_V            Cumulative error across the vertical overlap area
* @param paths_V        Optimal paths to reach the positions of the first row
* @param E_H            Cumulative error across the horizontal overlap area
* @param paths_H        Optimal paths to reach the positions of the first column
* @param ps             Patch size
* @param os             Overlap size
* @param boundary_mask  Mask used to redefine the ragged edges of 
*                       the patches in the blending step
*   @date 24/06/2014
*/
int create_Lshape_mask(Image *E_V, Image *E_H, int ps, int os, 
                       Image *boundary_mask, Image paths_V, Image paths_H)
{
    int xmin, ymin, x0, y0, ii, jj, counter, num_rand;
    int diag[os];
    float minE;
    xmin = 0;
    ymin = 0;
    minE = E_V->img[0] + E_H->img[0];
    for (ii=1; ii<os; ii++)
    {
        if (E_V->img[ii*E_V->w + ii] + E_H->img[ii*E_H->w + ii] < minE)
        {
            minE = E_V->img[ii*E_V->w + ii] + E_H->img[ii*E_H->w + ii];
        }
    }
    counter = 0;
    for (ii=0; ii<os; ii++)
    {
        if (E_V->img[ii*E_V->w + ii] + E_H->img[ii*E_H->w + ii] == minE)
        {
            diag[counter] = ii;
            counter++; 
        }
    }
    assert(counter>0);
    assert(counter<=os);
    num_rand = random_number(NULL) % counter;
    assert(num_rand<counter);
    assert(num_rand<os);
    xmin = diag[num_rand];
    ymin = xmin;
    x0 = xmin;
    y0 = ymin;
    //* initialize L-shape mask *//
    for (ii=0; ii<=xmin; ii++)
        for(jj=0; jj<=ymin; jj++)
    {
        boundary_mask->img[ii*boundary_mask->w + jj] = 1;
    }
    
    //* complete L-shape boundary mask*//
    for (ii=x0+1; ii<ps; ii++)
    {
        assert((ii-1)*(int)paths_V.w + ymin < (int)(paths_V.w*paths_V.h));
        ymin = ymin + paths_V.img[(ii-1)*(int)paths_V.w + ymin];
        for (jj=0; jj<=ymin; jj++)
        {
            assert(ii*boundary_mask->w + jj < 
                   boundary_mask->h*boundary_mask->w);
            boundary_mask->img[ii*boundary_mask->w + jj] = 1;
        }
    }
    for (jj=y0+1; jj<ps; jj++)
    {
        assert(xmin*(int)paths_H.w + (jj-1) < (int)(paths_H.w*paths_H.h));
        xmin = xmin + paths_H.img[xmin*(int)paths_H.w + (jj-1)];
        for (ii=0; ii<=xmin; ii++)
        {
            assert(ii*boundary_mask->w + jj < 
                   boundary_mask->h*boundary_mask->w);
            boundary_mask->img[ii*boundary_mask->w + jj] = 1;
        }
    }
    return 0;
}

/**
*   This method computes the boundary mask needed 
*   to then blend correctly the patch in the
*   output texture.
*   @author Lara Raad
*   @param src_im       The input texture sample
*   @param out_im       The synthesized texture
*   @param os           The size of the overlap area betwenn patches
*   @param ps           The size of the patches
*   @param temp         Left top corner position of the patch under 
*                       construction in the output image
*   @param patch_src    Corner position of the patch to copy from the 
*                       input texture sample
*   @date 22/06/2014
*/
int compute_boundary_mask(Image src_im, Image out_im, Image *boundary_mask,
                          int os, int ps, Corner * temp, Corner patch_src)
{
    int x, y;
    x = temp->x;
    y = temp->y;
    
    if(os==1)
    {
        int i,j;
        if (x > 0 && y > 0)                    //* L-shape overlap *//
        {
            for(i=0;i<ps;i++)
            {
                boundary_mask->img[i*boundary_mask->w] = 1;
            }
        }
        else if (x == 0 && y > 0)            //* vertical overlap *//
        {
            for(j=0;j<ps;j++)
            {
                boundary_mask->img[j] = 1;
            }
        }
        else if (y == 0 && x > 0)            //* horizontal overlap *//
        {
            boundary_mask->img[0] = 1;
            for(i=1;i<ps;i++)
            {
                boundary_mask->img[i*boundary_mask->w] = 1;
            }
            for(j=1;j<ps;j++)
            {
                boundary_mask->img[j] = 1;
            }
        }
    }
    else
    {
        if (x > 0 && y > 0)                    //* L-shape overlap *//
        {
            Image E_V, E_H, paths_V, paths_H;
            initialize_image(&paths_V, ps-1, os, 1);
            initialize_image(&paths_H, os, ps-1, 1);
            initialize_image(&E_V, ps, os, 1);
            initialize_image(&E_H, os, ps, 1);
            cumulative_vertical_error(&out_im, src_im, &E_V, &paths_V, 
                                      ps, os, x, y, patch_src);
            cumulative_horizontal_error(&out_im, src_im, &E_H, &paths_H, 
                                        ps, os, x, y, patch_src);

            //* create L-shape blending mask *//
            create_Lshape_mask(&E_V, &E_H, ps, os, boundary_mask, 
                               paths_V, paths_H);

            //* release memory *//
            free(E_V.img);
            free(E_H.img);
            free(paths_V.img);
            free(paths_H.img);
        } 
        else if (x == 0 && y > 0)            //* vertical overlap *//
        {
            Image E_V, paths_V;
            initialize_image(&paths_V, ps-1, os, 1);
            initialize_image(&E_V, ps, os, 1);
            cumulative_vertical_error(&out_im, src_im, &E_V, &paths_V, ps, 
                                      os, x, y, patch_src);

            //* create vertical blending mask *//
            create_vertical_mask(&E_V, ps, os, boundary_mask, paths_V);

            //* release memory *//
            free(E_V.img);
            free(paths_V.img);
        }
        else if (y == 0 && x > 0)            //* horizontal overlap *//
        {
            Image E_H, paths_H;
            initialize_image(&paths_H, os, ps-1, 1);
            initialize_image(&E_H, os, ps, 1);
            cumulative_horizontal_error(&out_im, src_im, &E_H, 
                                        &paths_H, ps, os, x, y, patch_src);

            //* create horizontal blending mask *//
            create_horizontal_mask(&E_H, ps, os, boundary_mask, paths_H);

            //* release memory *//
            free(E_H.img);
            free(paths_H.img);
        }
    }

    return 0;
}
