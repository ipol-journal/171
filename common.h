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

/**
*   @file common.h
*   @brief this header file contains all definitions.
*
*   @author Lara Raad
*
*   @date   24/05/2012
*/

#ifndef  COMMON_H_
#define  COMMON_H_

#include <fftw3.h>

typedef struct image
{
    float *img;
    size_t h,w,c;
} Image;

typedef struct image_d
{
    double *img;
    size_t h,w,c;
} Image_D;

typedef struct image_ld
{
    long double *img;
    size_t h,w,c;
} Image_LD;

typedef struct image_ui
{
    unsigned long *img;
    size_t h,w,c;
} Image_UI;

typedef struct corner
{
    int x,y;
} Corner;

typedef int bool;

Image initialize(Image src_img, int block_sz, int overlap_sz, int blocks_row,
                 int blocks_col, Corner *patch_src);
void initialize_image(Image *im, int rows, int cols, int channels);
void initialize_image_D(Image_D *im, int rows, int cols, int channels);
void initialize_image_LD(Image_LD *im, int rows, int cols, int channels);
void initialize_image_UI(Image_UI *im, int rows, int cols, int channels);
void create_position_map(Image position_map);
int crop_image(Image *im_crop, Image *im, int M, int N, int C);
int update_synth_map(Image * synth_map, Image * position_map, 
                     Corner * leftTopCorner_src, Corner * leftTopCorner_output,
                     int patch_sz, Image boundary_mask);
char* remove_ext_filename(char *filename);
int im_fft(float *img, fftwl_complex *fft_im, int h, int w);

#endif
