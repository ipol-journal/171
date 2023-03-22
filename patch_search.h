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
 * @file patch_search.h
 * This file contains all the functions for the patch search procedure that
 * computes the distance of a patch to all input patches and then randomly
 * pick on of the patches that are close to the minimal distance.
 *
 * @author Lara Raad, Bruno Galerne
 * @date   2014/06/11
 */

#ifndef  PATCH_SEARCH_H_
#define  PATCH_SEARCH_H_

#include "common.h"

void compute_input_patch_norms(Image src_im, 
        long double *patch_normsH, long double *patch_normsV, 
        long double *patch_normsL, int patch_sz, int overlap_sz);

void patch_search(Image src_im, Image out_im, Corner temp, int patch_sz, 
        long double *src_patch_normsH, long double *src_patch_normsV, 
        long double *src_patch_normsL, fftwl_complex *fft_srcR,
        fftwl_complex *fft_srcG, fftwl_complex *fft_srcB, 
        int h_dict, int w_dict, Corner *patch, float tol);

#endif
