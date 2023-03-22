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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <fftw3.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

#include "quilting.h"
#include "patch_search.h"
#include "boundary_cut.h"
#include "blending.h"
#include "random_number.h"

#define MAX_DIST     FLT_MAX

/**
*   This method synthesizes a texture image quilting patches taken from an
*    input sample texture. 
*   @author Lara Raad
*   @param src_im     The input texture sample
*   @param wp         The size of the patches
*   @param wo         The size of the overlap area betwenn patches
*   @param p_row      Number of patches per row to paste in the synthesized img
*   @param p_col      Number of patches per col to paste in the synthesized img
*   @param tol        The threshold used to select patches from the input img
*   @date 24/05/2012
*/
Image quilting_synthesis(Image src_im, int wp, int wo, int p_row, int p_col, 
                         float tol, Image * position_map, Image * synth_map,
                         long * seed)
{
    Image out_im;// synthesized image 
    long double *patch_normsH, *patch_normsV, *patch_normsL;
    Corner patch_src;// position (x,y) of the patch to paste in the output img
    Corner temp;// temporary position (x,y) of the patch under construction
    int step_size;// step size between each patch to paste in the output img
    int ii, jj, inc_s, h_s, w_s, w_o;
    int fft_size; // size of the r2c ffts for the source image

    int Nc = p_col;
    int Nr = p_row;
    
    /*************************************************************************/
    /*** INITIALIZATION                                                    ***/
    /*************************************************************************/
    
    //* Initialize ouput image : out_im *//
    out_im  = initialize(src_im, wp, wo, p_row, p_col,&patch_src);    
    
    //* Initiialize synthesis map *//
    initialize_image(synth_map, out_im.h, out_im.w, 3);
    
    h_s = src_im.h;
    w_s = src_im.w;
    w_o = out_im.w;
    
    inc_s = h_s*w_s;
    step_size = wp - wo;
    fft_size = h_s*((int) (w_s/2) + 1);
    
    int h_dict, w_dict;
    h_dict = h_s-wp+1;
    w_dict = w_s-wp+1;
    
    // create the norm of each valid patch in src_im
    // initialize the image of the patches norms from src_im (horizontal ov)
    patch_normsH = (long double *) malloc(inc_s*sizeof(long double));
    // initialize the image of the patches norms from src_im (vertical ov)
    patch_normsV = (long double *) malloc(inc_s*sizeof(long double));
    // initialize the image of the patches norms from src_im (L-shape ov)
    patch_normsL = (long double *) malloc(inc_s*sizeof(long double));

    compute_input_patch_norms(src_im, patch_normsH, 
                              patch_normsV, patch_normsL, wp, wo);
    
    //create the fft of the src_im for each channel
    fftwl_complex *fft_srcR, *fft_srcG, *fft_srcB;
    fft_srcR = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*(fft_size));
    fft_srcG = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*(fft_size));
    fft_srcB = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*(fft_size));
    im_fft(src_im.img, fft_srcR, h_s, w_s);
    im_fft(src_im.img + inc_s, fft_srcG, h_s, w_s);
    im_fft(src_im.img + 2*inc_s, fft_srcB, h_s, w_s);
    
    Image aux_mask;// Mask of all zeros to initialize the synthesis map
    initialize_image(&aux_mask, wp, wp, 1);

    temp.x = 0;
    temp.y = 0;
    update_synth_map(synth_map, position_map, &patch_src, &temp, wp, aux_mask);
    
    /*************************************************************************/
    /*** MAIN SYNTHESIS LOOP (RUN IN PARALLEL)                             ***/
    /*************************************************************************/

    for(jj=0; jj<(int) (Nr + 2*(Nc-1) -1)*step_size + 1; jj+=step_size)
    {
        #pragma omp parallel
        {
            /* initialize random generator seed for each thread */
            long s;
            if ( seed != NULL ) s = *seed ^ omp_get_thread_num();
            else                s = time(NULL) ^ omp_get_thread_num();
            random_number(&s);

            #pragma omp for
            for(ii=0; ii<(int) (Nc-1)*step_size+1; ii+=step_size)
            {
                //not for the seed patch
                if ((ii>0 || jj>0) && (jj-2*ii >= 0) && (jj-2*ii < w_o-wp+1))
                {
                    Corner temp;
                    Corner patch_pos;
                    Image boundary_mask;
                
                    //* Initialize boundary mask *//
                    initialize_image(&boundary_mask, wp, wp, 1);
                
                    temp.x = ii;
                    temp.y = jj-2*ii;
                
                    // Patch search: Randomly pick a patch that
                    // is close to the minimal distance
                    patch_search(src_im, out_im, temp, wp, patch_normsH, 
                                 patch_normsV, patch_normsL, fft_srcR, fft_srcG,
                                 fft_srcB, h_dict, w_dict, &patch_pos, tol);

                    // compute boundary mask
                    compute_boundary_mask(src_im, out_im, &boundary_mask, 
                                          wo, wp, &temp, patch_pos);

                    // Updating synthesis map
                    update_synth_map(synth_map, position_map, &patch_pos, 
                                     &temp, wp, boundary_mask);
                
                    // blending new patch into out_im
                    blend_patch(boundary_mask, &out_im, src_im, patch_pos, 
                                temp, wp, wo);

                    free(boundary_mask.img);
                }
            }
        }
    }
    
    //free memory
    free(patch_normsH);
    free(patch_normsV);
    free(patch_normsL);
    free(aux_mask.img);
    fftwl_free(fft_srcR);
    fftwl_free(fft_srcG);
    fftwl_free(fft_srcB);

    return out_im;
}
