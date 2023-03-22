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
*   @file boundary_cut.h
*   @brief this header file contains all necessary
*   functions to perform the boundary mask to thenn 
*   blend the patch under construction into the
*   synthesized image.
*
*   @author Lara Raad
*
*   @date   24/06/2014
*/
#ifndef  BOUNDARY_CUT_H_
#define  BOUNDARY_CUT_H_

#include "common.h"

int cumulative_vertical_error(Image *out_im, Image src_im, Image *E, 
            Image *paths, int ps, int os, int x, int y, Corner patch_src);

int cumulative_horizontal_error(Image *out_im, Image src_im, Image *E, 
            Image *paths, int ps, int os, int x, int y, Corner patch_src);

int compute_boundary_mask(Image src_im, Image out_im, Image *boundary_mask, 
            int os, int ps, Corner * temp, Corner patch_src);

int create_vertical_mask(Image *E_V, int ps, int os, Image *boundary_mask, 
            Image paths_V);

int create_horizontal_mask(Image *E_H, int ps, int os, Image *boundary_mask, 
            Image paths_H);

int create_Lshape_mask(Image *E_V, Image *E_H, int ps, int os, 
            Image *boundary_mask, Image paths_V, Image paths_H);

#endif
