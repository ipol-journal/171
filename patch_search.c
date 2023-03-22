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
#include <fftw3.h>
#include <assert.h>
#include <float.h>

#include "patch_search.h"
#include "random_number.h"

#define MAX_DIST        FLT_MAX
#define MIN_DIST        1.0
#define MAX_INT         10000000
#define SMOOTH 1

/**
 * @file distances.c 
 * This file contains all the functions for the patch search procedure that 
 * computes the distance of a patch to all input patches and then randomly
 * pick on of the patches that are close to the minimal distance.
 */

/**
 * Convention for image sizes
 * w_img = nx = width of image : Indexes x or j
 * h_img = ny = height of image : Indexes y or i
 * 
 */

/**
 * @brief This function computes the cross-correlation between img1 and img2 
 * using FFT, given the fft of img2 as input.
 * The result is a long double image having the same size (already allocated).
 * @param img1, img2, two pointers to float images
 * @param cross_corr image of the cross-correlation
 * @param nx, ny image sizes
 */
static void cross_correlation_fft2_known(float *img1, fftwl_complex *fft2, 
                                long double *cross_corr, size_t nx, size_t ny)
{    
    long double *ldimg;
    fftwl_complex *fft1;
    fftwl_plan plan_forward, plan_backward;
    int i, j, si, sj;
    size_t size = nx*ny;
    int fftsize = (int) (ny*((int) (nx/2)+1));
    long double re, im;
    int imsize = (int) (nx*ny);
    long double invimsize = 1./(long double) imsize;
    
    /* Memory allocation */
    ldimg = (long double *) malloc(size*sizeof(long double));
    #pragma omp critical (make_plan)
    fft1 = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*(fftsize));
    
    /* Copy symmetric version of img1 into a long double array
       and compute FFT */
    for(i=0; i < (int) ny; i++)
    {
         // Define (si, sj) = (-i,-j) modulo (nx,ny)
        si = (i==0) ? 0 : (ny-i);
        for(j=0; j < (int) nx; j++)
        {
            sj = (j==0) ? 0 : (nx-j);
            assert(nx*si+sj<size);
            ldimg[nx*si+sj] = (long double) img1[nx*i+j];
        }
    }
    #pragma omp critical (make_plan)
    plan_forward = fftwl_plan_dft_r2c_2d(ny,nx,ldimg,fft1,FFTW_ESTIMATE);
    fftwl_execute(plan_forward);

    /* Pointwisze multiplication of Fourier transforms */
    for(i=0; i<fftsize; i++)
    {
        re = fft1[i][0];
        im = fft1[i][1];
        fft1[i][0] = re*fft2[i][0] - im*fft2[i][1];
        fft1[i][1] = re*fft2[i][1] + im*fft2[i][0];
    }
    /* Backward FFT */
    #pragma omp critical (make_plan)
    plan_backward = fftwl_plan_dft_c2r_2d(ny,nx,fft1,cross_corr,FFTW_ESTIMATE);
    fftwl_execute(plan_backward);

    /* normalization of the convolution */
    for(i=0; i<(int) size; i++)
    {
        cross_corr[i] *= invimsize;
    }
    
    /* free memory */
    free(ldimg);
    #pragma omp critical (make_plan)
    fftwl_destroy_plan(plan_forward);
    #pragma omp critical (make_plan)
    fftwl_destroy_plan(plan_backward);
    #pragma omp critical (make_plan)
    fftwl_free(fft1);
}

/**
 * @brief Computes the squared norms of each RGB patches of the input according
 * to the 3 possible overlap regions.The computation is done in applying the 
 * cross-correlation between the extended mask and the sum over the 3 channels 
 * of the squares of the pixel values.
 * @param src_im input RGB Image
 * @param patch_normsH long double image for the horizontal overlap
 *                     (must be allocated)
 * @param patch_normsV long double image for the vertical overlap 
 *                       (must be allocated)
 * @param patch_normsL long double image for the L-shape overlap 
 *                        (must be allocated)
 * @param patch_sz patch size 
 * @param overlap_sz size of the overlap at the edge of the patch
 * @author Lara Raad, Bruno Galerne
 */
void compute_input_patch_norms(Image src_im, 
                    long double *patch_normsH, long double *patch_normsV, 
                    long double *patch_normsL, int patch_sz, int overlap_sz)
{
    float *mask_ext;
    long double *sq_rgb_src_im;
    fftwl_complex *fft_sq_rgb_src_im;
    fftwl_plan plan_forward;
    int i,j;
    int size = (int) (src_im.h*src_im.w);
    int fftsize = (int) src_im.h*((int) (src_im.w/2)+1);
    // allocate memory
    mask_ext = (float *) calloc(size, sizeof(float));
    for(i=0; i< size; i++) mask_ext[i] = 0.;
    sq_rgb_src_im = (long double *) malloc(size*sizeof(long double));
    #pragma omp critical (make_plan)
    fft_sq_rgb_src_im = 
        (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*(fftsize));
    
    // Compute sq_rgb_src_im : each pixel is the sum of the squares of
    // the pixels of the 3 RGB channels
    for(i=0; i<size; i++)
    {
        sq_rgb_src_im[i] = (long double) src_im.img[i]
                         * (long double) src_im.img[i]
                       +   (long double) src_im.img[i+size]
                         * (long double) src_im.img[i+size]
                       +   (long double) src_im.img[i+2*size]
                         * (long double) src_im.img[i+2*size];
    }
    // Compute the fft of this image
    #pragma omp critical (make_plan)
    plan_forward = fftwl_plan_dft_r2c_2d(src_im.h,src_im.w,sq_rgb_src_im,
                                         fft_sq_rgb_src_im,FFTW_ESTIMATE);
    fftwl_execute(plan_forward);

    // Horizontal overlap 
    for(i=0;i<patch_sz;i++)
    {
        for(j=0;j<patch_sz;j++)
        {
            if(i<overlap_sz)
                mask_ext[src_im.w*i+j] = 1.; 
            else 
                mask_ext[src_im.w*i+j] = 0.;
        }
    }
    // compute the norm of all patches with cross-correlation
    cross_correlation_fft2_known(mask_ext, fft_sq_rgb_src_im, 
                                 patch_normsH, src_im.w, src_im.h);
    
    // Vertical overlap
    for(i=0;i<patch_sz;i++)
    {
        for(j=0;j<patch_sz;j++)
        {
            if(j<overlap_sz) 
                mask_ext[src_im.w*i+j] = 1.; 
            else 
                mask_ext[src_im.w*i+j] = 0.;
        }
    }    
    // compute the norm of all patches with cross-correlation
    cross_correlation_fft2_known(mask_ext, fft_sq_rgb_src_im, 
                                 patch_normsV, src_im.w, src_im.h);
    
    // L-shape overlap
    for(i=0;i<patch_sz;i++)
    {
        for(j=0;j<patch_sz;j++)
        {
            if(i<overlap_sz) 
                mask_ext[src_im.w*i+j] = 1.; 
            else if (j<overlap_sz) 
                mask_ext[src_im.w*i+j] = 1.; 
            else
                mask_ext[src_im.w*i+j] = 0.; 
        }
    }
    // compute the norm of all patches with cross-correlation
    cross_correlation_fft2_known(mask_ext, fft_sq_rgb_src_im, 
                                 patch_normsL, src_im.w, src_im.h);
    
    // free memory
    free(mask_ext);
    free(sq_rgb_src_im);
    #pragma omp critical (make_plan)
    fftwl_destroy_plan(plan_forward);
    #pragma omp critical (make_plan)
    fftwl_free(fft_sq_rgb_src_im);
}



/** @brief This function computes the squared distance between a given patch
 * and all patches of another image. The squared distance is (implicitly) 
 * weighted by a binary mask on the patch domain.It does not appear because 
 * pixels that are not yet defined (outside the L-shape corner generally) are 
 * set to zero, so that multiplying by the binary mask does nàt change the
 * patch.
 * @param src_im     The input texture sample
 * @param out_im     The currently synthesized texture image 
 *                     (unknown pixels must be set to zero)
 * @param temp       The patch from the output image to compare with those in 
 *                     the texture sample (patch to synthesize)
 * @param patch_sz   The pathc size
 * @param src_patch_normsH, src_patch_normsV, src_patch_normsL long double 
 *                     arrays corresponding to the norm of input patches for 
 *                     the 3 possible overlap regions
 * @param fft_srcR, fft_srcG,fft_srcB  long double Fourier transforms of 
 *                     the 3 input channels
 * @param distances  Image of distances, WARNING: suppposed to be set to 0
 *                     at initialization (result)
 * @author Lara Raad, Bruno Galerne
 * @date 2014/06/10
 */
static void compute_distance_to_input_patches(Image src_im, Image out_im, 
            Corner temp, int patch_sz, long double *src_patch_normsH, 
            long double *src_patch_normsV, long double *src_patch_normsL, 
            fftwl_complex *fft_srcR,fftwl_complex *fft_srcG, 
            fftwl_complex *fft_srcB, Image_LD *distances)
{
    long double *cross_corr;
    float *img1;
    int i,j,k;
    float patchij;
    long double normpatch;
    long double *src_patch_norms;
    fftwl_complex *fft2;
    
    // allocate memory
    img1 = (float *) calloc(src_im.h*src_im.w, sizeof(float));
    cross_corr = (long double *) malloc(src_im.h*src_im.w*sizeof(long double));
    
    for(k=0; k<(int)src_im.c; k++)
    { 
        normpatch = 0.;
        // extend patch channel times mask into a large image 
        // and compute squared norm of patch
        for(i=0;i<patch_sz;i++)
        {
            for(j=0;j<patch_sz;j++)
            {
                patchij = out_im.img[out_im.w*(temp.x+i)+temp.y+j+
                                     k*out_im.w*out_im.h];
                img1[src_im.w*i+j] =  patchij;
                normpatch += (long double) patchij * (long double) patchij;
            }
        }
        
        // compute cross-correlation
        // change input fft into the concatenation of the 3 arrays ?
        if (k==0) fft2 = fft_srcR; 
        else if (k==1) fft2 = fft_srcG; 
        else fft2 = fft_srcB;
        cross_correlation_fft2_known(img1, fft2, cross_corr,
                                     src_im.w, src_im.h);
        
        // add to distances
        for(i=0; i< (int) (src_im.h*src_im.w); i++)
        {
            // WARNING: distances is supposed to be allocated and be set to 0 
            // before calling this function
            distances->img[i] += normpatch - 2.0*cross_corr[i]; 
        }
    }
    // Add norm of input RGB patches
    // choose overlap case
    // L-shape overlap : default case
    src_patch_norms = src_patch_normsL;
    if (temp.y == 0) 
        src_patch_norms = src_patch_normsH;// horizontal overlap (first column)
    if (temp.x == 0) 
        src_patch_norms = src_patch_normsV;// vertical overlap (first row)
        
    for(i=0; i<(int) (src_im.h*src_im.w); i++)
    {
        distances->img[i] += src_patch_norms[i];
    }
    
    /* free memory */
    free(img1);
    free(cross_corr);

}


/**
 *   This function picks a random patch within all the patches of the source 
 *   image that verify dist(p, current) < dist_min*(1+tolerance) with dist_min
 *   the minimum distance between the current patch and the all the patches
 *   from the input image
 *   @author Lara Raad, Bruno Galerne
 *   @param distances    Image containing the distance between the current 
 *                       patch and the patch of source image (distance has 
 *                       the same size as the input image, but the bottom-right
 *                       border must be discarded)
 *   @param h, w         size of the domain of admissible patches
 *   @param patch        The position of the picked patch (Corner type)
 *   @param tol          The threshold used to select patches from the input
 *   @date 2014/06/10
 */
static void pick_patch(Image_LD *distances, int h, int w, 
                       Corner *patch, float tol)
{
    Corner *cands;//patches veryfing (dist <= min_dist*(1+tolerance))
    int rand_num, pos_dist;
    int count_cands;//amount of patches veryfing dist <= min_dist*(1+tolerance)
    long double min_dist; // minimal distance between the current patch and 
                          // those from the input image
    
    // compute the minimal distance
    min_dist = MAX_DIST;//distances->img[0];
    for(int i=0; i<h; i++)
    {
        for (int j=0; j<w; j++)
        {
            pos_dist = i*distances->w + j;
            assert(distances->img[pos_dist]==distances->img[pos_dist]);
            if(distances->img[pos_dist] < 
                            min_dist && distances->img[pos_dist] > 1.)
                {
                    min_dist = distances->img[pos_dist];
                }
        }
    }
    
    // create the list of patches verifying dist <= min_dist*(1+tolerance)    
    cands = (Corner *) calloc(h*w, sizeof(Corner));
    count_cands = 0;
    for(int i=0; i<h; i++)
    {
        for (int j=0; j<w; j++)
        {
            pos_dist = i*distances->w + j;
            assert(pos_dist < (int) (distances->w*distances->h));
            if(distances->img[pos_dist] <= min_dist*(1+tol))
            {
                cands[count_cands].x = i;
                cands[count_cands].y = j;
                count_cands++;
            }
        }
    }
    assert(count_cands>0);
    assert(count_cands<=h*w);
    // choose a random patch within the list of candidates
    rand_num = random_number(NULL) % count_cands;
    assert(rand_num<count_cands);
    patch->x = cands[rand_num].x;
    patch->y = cands[rand_num].y;
    
    //free memory
    free(cands);
}

/**
 * @brief This function executes the two steps of the patch search procedure:
 * 1. Compute the distance of the current patch 
 *    to all the admissible input patches
 * 2. Randomly pick one of the patch that is close to the minimal distance
 * @param src_im     The input texture sample
 * @param out_im     The currently synthesized texture image
 *                     (unknown pixels must be set to zero)
 * @param temp       The patch from the output image to compare with those
 *                     in the texture sample (patch to synthesize)
 * @param patch_sz   The pathc size
 * @param src_patch_normsH, src_patch_normsV, src_patch_normsL long double 
 *                     arrays corresponding to the norm of input patches for
 *                     the 3 possible overlap regions
 * @param fft_srcR, fft_srcG,fft_srcB  long double Fourier transforms 
 *                     of the 3 input channels
 * @param h, w       size of the domain of admissible patches
 * @param patch      The position of the picked patch (Corner type)
 * @param tol        The threshold used to select patches from the input image
 * @date 2014/06/10
 */
void patch_search(Image src_im, Image out_im, Corner temp, int patch_sz, 
        long double *src_patch_normsH, long double *src_patch_normsV, 
        long double *src_patch_normsL, fftwl_complex *fft_srcR,
        fftwl_complex *fft_srcG, fftwl_complex *fft_srcB, int h, int w, 
        Corner *patch, float tol)
{
    Image_LD distances;
    
    /* memory allocation*/
    initialize_image_LD(&distances, src_im.h, src_im.w, 1);
    
    /* compute distance to all input patches */
    compute_distance_to_input_patches(src_im, out_im, temp, patch_sz, 
        src_patch_normsH, src_patch_normsV, src_patch_normsL, 
        fft_srcR, fft_srcG, fft_srcB, &distances);
    
    /* pick a random patch */
    pick_patch(&distances, h, w, patch, tol);

    /* free memory */
    free(distances.img);
}
