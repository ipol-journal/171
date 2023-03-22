README file for image quilting for texture synthesis

 - Compilation:

   make           : normal compilation
   make openmp    : compilation with OpenMP parallelism support

 - Required libraries:

libpng, fftw3

 - Test:

./quilting test/input.png test/output.png test/input_map.png test/output_map.png

 - Usage:

./quilting [-p patch_sz] [-r ratio] [-t tolerance] [-o overlap] [-s seed] src_im.png out_im.png position_map.png synthesis_map.png

Required parameters:
   src_im.png        :   name of the input PNG image
   out_im.png        :   name of the output PNG image
   position_map.png  :   name of the position map PNG image
   synthesis_map.png :   name of the synthesis map PNG image

Optional parameters:
   patch_sz        :   int > size(src_im.png) that specifies the size of 
                       the patch used to synthesize out_im.png
   ratio           :   float >= 1 that specifies the increase of size of 
                       out_im.png relatively to the size of src_im.png
   tolerance       :   float > 0 to specify the error tolerance when comparing
                       patches in src_im.png
   overlap         :   float 0 < o < 1 to specify the patch overlap ratio
   seed            :   unsigned int to specify the seed for the random number 
                       generator (seed = time(NULL) by default)

- Copyright and License

Copyright (c) 2016 Lara Raad <lara.raad@cmla.ens-cachan.fr>,
                   Bruno Galerne <bruno.galerne@parisdescartes.fr>

Quilting is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

Smooth Contours is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
