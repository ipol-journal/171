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
*   @file quilting.h
*   @brief this header file contains the main method of
*   image quilting method.
*
*   @author Lara Raad
*
*   @date   24/05/2012
*/
#ifndef  QUILTING_H_
#define  QUILTING_H_

#include "common.h"

Image quilting_synthesis(Image src_im, int patch_sz, int overlap_sz, 
                         int blocks_row, int blocks_col, float tolerance, 
                         Image * position_map, Image * synth_map,
                         long * seed);

#endif
