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
#include <string.h>
#include <fftw3.h>

#include "common.h"
#include "quilting.h"
#include "random_number.h"

/**
*   This method initializes the output image taking a random patch from the 
*   input texture and placing it on the left-top corner of the output image.
*   @author Lara Raad
*   @param src_img     The input texture sample
*   @param wp          The size of the patches
*   @param wo          The size of the overlap area betwenn patches
*   @param p_row       Number of patches per row in the output image
*   @param p_col       Number of patches per column in the output image
*   @param patch_src   Top left corner position of the patch from the source
*   @date 24/05/2012
*/
Image initialize(Image src_img, int wp, int wo, int p_row, 
                 int p_col, Corner *patch_src){
    Image img;
    int i, j, k, rows, cols, channels, pos, pos_s, inc, inc_s;
    rows = p_row*wp - (p_row - 1)*wo;
    cols = p_col*wp - (p_col - 1)*wo;
    channels = src_img.c;

    initialize_image(&img, rows, cols, channels);
      
    // Paste random block from input texture in left top corner
    // of the output image.
     
    // Random positions in [0,h-wp)x[O,w-wp)
    int x_rand = random_number(NULL) % (src_img.h-wp+1);
    assert(x_rand<(int)(src_img.h-wp+1));
    int y_rand = random_number(NULL) % (src_img.w-wp+1);
    assert(y_rand<(int)(src_img.w-wp+1));
    
    patch_src->x = x_rand;
    patch_src->y = y_rand;

    inc = img.h*img.w;
    inc_s = src_img.h*src_img.w;
    
    //paste the first block into the output image
    for(i=0; i < wp; i++)
    {
        for(j=0; j < wp; j++)
        {
            pos = i*img.w + j; 
            pos_s = (x_rand + i)*src_img.w +(y_rand +j);
            assert(pos_s < (int) (src_img.w*src_img.h));
            for(k=0; k<channels; k++)
            {
                img.img[pos + k*inc] = src_img.img[pos_s + k*inc_s];
            }
        }
    }
            
    return img;
}

/**
*   This method initializes images of type float. 
*   It sets the numeber of rows, columns and channels 
*   and sets the pixels values to zero.
*   @author Lara Raad
*   @param im       The image to initialize
*   @param rows     Number of rows of the image
*   @param cols     Number of columns of the image
*   @param channels Number of channels of the image
*   @date 24/05/2012
*/
void initialize_image(Image *im, int rows, int cols, int channels)
{
    int size;
    im->h = rows;
    im->w = cols;
    im->c = channels;
    size = rows*cols*channels;
    im->img = (float *) calloc(size, sizeof(float));
}

/**
*   This method initializes images of type double. It sets the number
*   of rows, columns and channels and sets the pixels values to zero.
*   @author Lara Raad
*   @param im       The image to initialize
*   @param rows     Number of rows of the image
*   @param cols     Number of columns of the image
*   @param channels Number of channels of the image
*   @date 24/05/2012
*/
void initialize_image_D(Image_D *im, int rows, int cols, int channels)
{
    int size;
    im->h = rows;
    im->w = cols;
    im->c = channels;
    size = rows*cols*channels;
    im->img = (double *) calloc(size, sizeof(double));
}

/**
*   This method initializes images of type long double. It sets the number
*   of rows, columns and channels and sets the pixels values to zero.
*   @author Lara Raad
*   @param im       The image to initialize
*   @param rows     Number of rows of the image
*   @param cols     Number of columns of the image
*   @param channels Number of channels of the image
*   @date 24/05/2012
*/
void initialize_image_LD(Image_LD *im, int rows, int cols, int channels)
{
    int size;
    im->h = rows;
    im->w = cols;
    im->c = channels;
    size = rows*cols*channels;
    im->img = (long double *) calloc(size, sizeof(long double));
}

/**
*   This method initializes images of type unsigned long. It sets the number
*    of rows, columns and channels and sets the pixels values to zero.
*   @author Lara Raad
*   @param im       The image to initialize
*   @param rows     Number of rows of the image
*   @param cols     Number of columns of the image
*   @param channels Number of channels of the image
*   @date 24/05/2012
*/
void initialize_image_UI(Image_UI *im, int rows, int cols, int channels)
{
    int size;
    im->h = rows;
    im->w = cols;
    im->c = channels;
    size = rows*cols*channels;
    im->img = (unsigned long *) calloc(size, sizeof(unsigned long));
}

/**
*   This method creates a position map that associates a different
*    color from a continuous colormap to every pixel of the source 
*    image
*   @author Lara Raad
*   @param pos_map     Position map image
*   @date 09/05/2013
*/
void create_position_map(Image pos_map)
{
    size_t count = 0;
    size_t inc = pos_map.h*pos_map.w;
    size_t i,j;
    
    for(i=0;i<pos_map.h;i++)
    {
        for(j=0;j<pos_map.w;j++)
        {
            pos_map.img[count] = (i*255)/pos_map.h;//Red
            pos_map.img[count+inc] = (j*255)/pos_map.w;//Green
            pos_map.img[count+2*inc] = (i*j*255)/(pos_map.w*pos_map.h);//Blue
            count++;
        }
    }
}

/**
*   This method crops an image to be of size MxNxC
*   @author Lara Raad
*   @param im       Image to crop
*   @param M,N,C    Size of the cropped image
*   @param im_crop  The cropped image
*   @date 09/08/2013
*/
int crop_image(Image *im_crop, Image *im, int M, int N, int C)
{
    if( (M > (int) im->h) || (N > (int) im->w) )
    {
        return 1;
    }
    
    int i,j,k;
    initialize_image(im_crop, M, N, C);
    for(k=0;k<C;k++)
    {
        for(i=0;i<M;i++)
        {
            for(j=0;j<N;j++)
            {
                im_crop->img[i*N+j + k*M*N] = 
                                        im->img[i*im->w+j + k*im->h*im->w];
            }
        }
    }        
    return 0;
}

/**
*   This method updates the synthesis map for each patch added
*    in the synthesized image
*   @author Lara Raad
*   @param synth_map        Synthesis map
*   @param pos_map          Color map of the texture sample
*   @param leftTop_src      Position of the patch to add taken in the src image
*   @param leftTop_output   Position in the ouput image where we add a block
*   @param pathc_sz         Patch size
*    @param binIm           Binary image used to indicate which pixels are form
*                           the old patch and wich pixels are from the new patch
*   @date 26/02/2014
*/
int update_synth_map(Image * synth_map, Image * pos_map, Corner * leftTop_src, 
                     Corner * leftTop_output, int wp, Image binIm)
{
    int kk, ll, pos_posMap, pos_synthMap, x_s, y_s, x_o, y_o, inc_o, inc_s;
    x_s = leftTop_src->x;
    y_s = leftTop_src->y;
    assert( (x_s < (int)((synth_map->w-wp+1)*(synth_map->h-wp+1))) && 
            (y_s < (int)((synth_map->w-wp+1)*(synth_map->h-wp+1))) );
    x_o = leftTop_output->x;
    y_o = leftTop_output->y;
    inc_s = pos_map->w*pos_map->h;
    inc_o = synth_map->w*synth_map->h;
    
    for (kk=0; kk<wp; kk++)
    {
        for (ll=0; ll<wp; ll++)
        {
            pos_posMap = (x_s + kk)*pos_map->w + (y_s + ll);
            pos_synthMap = (x_o + kk)*synth_map->w + (y_o + ll);
            assert(pos_posMap < (int)(pos_map->h*pos_map->w));
            assert(pos_synthMap < (int)(synth_map->h*synth_map->w));
            synth_map->img[pos_synthMap] = 
            binIm.img[kk*binIm.w+ll]*synth_map->img[pos_synthMap]
            + (1-binIm.img[kk*binIm.w+ll])*pos_map->img[pos_posMap];
            synth_map->img[pos_synthMap + inc_o] = 
            binIm.img[kk*binIm.w+ll]*synth_map->img[pos_synthMap + inc_o] 
            + (1-binIm.img[kk*binIm.w+ll])*pos_map->img[pos_posMap + inc_s];
            synth_map->img[pos_synthMap + 2*inc_o] = 
            binIm.img[kk*binIm.w+ll]*synth_map->img[pos_synthMap + 2*inc_o]
            + (1-binIm.img[kk*binIm.w+ll])*pos_map->img[pos_posMap + 2*inc_s];
        }
    }

    return 0;
}

/**
*   This method removes the extension of a filename
*   @author Lara Raad
*   @param filename    File name to which the extension is removed
*   @date 25/08/2014
*/
char* remove_ext_filename(char *filename)
{
    // get the string after the last '/' (when exists) and
    // remove extension .png
    char *file;
    if ((file = strrchr( filename,  '/'))==NULL)
        file = filename;
    else
        file = file +1;
    int i = strlen(file);
    while( (file[i-1] != '.') && (i>0))
    {
        i--;    
    }
    int posDot = i-1;
    char *fileNoExt = (char *) calloc(posDot, sizeof(char));
    strncpy(fileNoExt,file,posDot);
    return fileNoExt;    
}

/**
*   This method computes the 2-D FFT of an input image
*   @author Lara Raad
*   @param img            Input image to apply the 2-D FFT
*   @param fft_im         2-D FFT of the input image img
*   @param h,w            Dimensions of img (height and width)
*   @date 26/02/2014
*/
int im_fft(float *img, fftwl_complex *fft_im, int h, int w)
{
    fftwl_plan p_im;
    Image_LD src_ld;// copy of src_im in a long double version
    int pos_im;
        
    initialize_image_LD(&src_ld, h, w, 1);
    
    for(int i=0; i<(int)h; i++)
    {
        for(int j=0; j<(int)w; j++)
        {
            pos_im = i*w + j;
            src_ld.img[pos_im] = (long double)(img[pos_im]);
        }
    }
    
    // create plan
    #pragma omp critical (make_plan)
    p_im = fftwl_plan_dft_r2c_2d(h,w,src_ld.img,fft_im,FFTW_ESTIMATE);
    // execute fft
    fftwl_execute(p_im);
    
    free(src_ld.img);
    #pragma omp critical (make_plan)
    fftwl_destroy_plan(p_im);
    
    return 0;
}
