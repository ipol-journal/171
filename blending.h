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
*   @file blending2.h
*   @brief this header file contains all necessary
*   functions to perform the boundary mask to thenn 
*   blend the patch under construction into the
*   synthesized image.
*
*   @author Lara Raad
*
*   @date   24/06/2014
*/
#ifndef  BLENDING2_H_
#define  BLENDING2_H_

#include "common.h"

int blend_patch(Image boundary_mask, Image *out_im, Image src_im, 
                Corner patch_src, Corner temp, int ps, int os);

int smooth_mask(Image *boundary_mask, int ps, int os, Corner temp);

int filter_mask(int X1[],int Y1[], int X2[], int Y2[], int ps, 
                Image *boundary_mask, float kernel[], int center);

#endif
