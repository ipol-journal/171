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
#include <time.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <unistd.h>

#include "quilting.h"
#include "random_number.h"
#include "io_png.h"

//* print a message and abort the execution *//
#define FATAL(MSG)              \
do {                            \
     fprintf(stderr, MSG "\n"); \
     abort();                   \
   } while (0);

static int display_usage()
{
    printf("\nquilting_main -p wp -r ratio -t tol -o overlap -s seed "
           "src_im.png out_im.png pos_map.png synthesis_map.png\n\n");
    printf("Required parameters:\n");
    printf("   src_im.png          : name of the input PNG image\n");
    printf("   out_im.png          : name of the output PNG image\n");
    printf("   pos_map.png         : name of the position map PNG image\n");
    printf("   synthesis_map.png   : name of the synthesis map PNG image\n\n");
    printf("Optionnal parameters:\n");
    printf("   wp      :   int > size(src_im.png) that "
           "specifies the size of the patch used to synthesize "
           "out_im.png\n");
    printf("   ratio   :   float >= 1 that specifies the increase "
           "of size of out_im.png relatively to the size "
           "of src_im.png\n");
    printf("   tol     :   float > 0 to specify the distance "
           "tolerance when comparing patches in src_im.png\n");
    printf("   overlap :   float 0 < o < 1 to specify the patch "
           "overlap ratio\n");
    printf("   seed    :   unsigned int to specify the seed "
           "for the random number generator "
           "(seed = time(NULL) by default)\n");

    return (1);
}

int main(int argc, char *argv[])
{       
    //* 
    //* Getting parameters: patch size, overlap size, number 
    //* of blocks per row, number of blocks per column, tolerance 
    //* for distance
    //*

    char *fname_in;             //* input file name *//
    char *fname_out;            //* output file name *//
    char *fname_positionMap;    //* position map file name *//
    char *fname_synthesisMap;   //* synthesis map file name *//

    Image src_im;               //* input data *// 
    Image out_im;               //* output data *//
    int wp;                     //* patch size *//
    float ratio;                //* ratio for the
                                //* output size (must be >1) *//
    float tol;                  //* distance tolerance *//
    float overlap;              //* patch overlap ratio *//
    int wo;                     //* overlap size between blocks *//
    int p_row;                  //* number of patches per row *//
    int p_col;                  //* number of patches per column *//
    
    long seed = 0;              /* seed for the random number
                                 * generator */
    int p_flag = 0;
    int s_flag = 0;
    int c;                      /* getopt flag */
    
    //* "default" value initialization *//
    wp = 40;
    ratio = 1.;
    tol = .1;
    overlap = .25;
    
    //* process the options and parameters *//
    
    while ((c = getopt(argc, argv, "p:r:t:o:s:h")) != -1) {
        switch (c) {
        case 'p':
            //* patch size specified *//
            p_flag = 1;
            wp = atof(optarg);
            break;
        case 'r':
            //* ratio specified *//
            ratio = atof(optarg);
            if (ratio < 1)
                FATAL("Option usage: The ratio "
                      "specified by -r must be greater than 1.");
            break;
        case 't':
            //* distance tolerance specified *//
            tol = atof(optarg);
            if (tol <= 0.)
                FATAL("Option usage: The tolerance error  "
                      "specified by -t must be greater than 0.");
            break;
        case 'o':
            //* overlap ratio specified *//
            overlap = atof(optarg);
            if (overlap <= 0. || overlap >= 1.)
                FATAL("Option usage: The overlap ratio  "
                      "specified by -o must be between 0. and 1.");
            break;
        case 's':
            //* seed specified *//
            s_flag = 1;
            seed = atol(optarg);
            break;
        case 'h':
            //* display usage *//
            display_usage();
            return (-1);
        default:
            abort();
        }
    }
    
    //* process the non-option parameters *//
    if (4 > (argc - optind)) {
        printf("The input and output image file names are missing\n\n");
        display_usage();
        return (-1);
    }
    fname_in = argv[optind++];
    fname_out = argv[optind++];
    fname_positionMap = argv[optind++];
    fname_synthesisMap = argv[optind++];
    
    //* Reading input image *//
    src_im.img = io_png_read_f32(fname_in, &src_im.w, &src_im.h, &src_im.c);

    //* Initilization of the random number generator *//
    if (s_flag == 0) {
        long s = time(NULL);
        random_number(&s);
    }
    else {
        random_number(&seed);
    }

    //* deal with patch size,*//
    //*if higher than the min lentgh of the input abort*//
    if((p_flag==1) && ((wp > (int) src_im.h) || (wp > (int) src_im.w)))
    {
        printf("Conflict: the patch size should be smaller "
               "than the minimal length of the input image\n\n");
        display_usage();
        return (-1);
    }
    
    //* deal with patch size, if equal to zero abort*//
    if((p_flag==1) && (wp==0))
    {
        printf("Conflict: the patch size should be greater than zero\n\n");
        display_usage();
        return (-1);
    }

    //* compute overlap size *//
    if ((wo = (int) round(wp*overlap))<=0 || wo>=wp)
        FATAL("Conflict: The patch overlap size should be at least one pixel "
              "and less than the patch size");
    
    //* compute number of patches per row and per column needed*//
    //* to synthesize a (factor*M) x (factor*N) modulo (wp) image *//
    int M = roundf(ratio*src_im.h);
    int N = roundf(ratio*src_im.w);
    int C = src_im.c;
    float a;    
    if (src_im.w < src_im.h)
        a = M - wo;
    else
        a = N - wo;
    float b;
    if (src_im.h < src_im.w)
        b = N - wo;
    else
        b = M - wo;
    float d = wp - wo;
    p_row = (int) ceil(b/d);
    p_col = (int) ceil(a/d);

    Image pos_map;// each pixel in the source image is associated to
                   // a different color from a continuous colormap
    Image synth_map;// shows for each synthesized patch its position
                     // in the source image

    //* Create position map*//
    initialize_image(&pos_map, src_im.h, src_im.w, 3);
    create_position_map(pos_map);

    //* Image synthesis *//
    if (s_flag == 0)
        out_im = quilting_synthesis(src_im, wp, wo, p_row, p_col, 
                                    tol, &pos_map, &synth_map, NULL);
    else
        out_im = quilting_synthesis(src_im, wp, wo, p_row, p_col, 
                                    tol, &pos_map, &synth_map, &seed);
    
    //* Crop synthesis map *//
    Image synth_map_crop;
    if (crop_image(&synth_map_crop, &synth_map, M, N, 3))
    {
        FATAL("The size of the cropped image is bigger "
                      "than the image to crop");
    }

    //* Crop output image *//
    Image out_im_crop;
    if (crop_image(&out_im_crop, &out_im, M, N, C))
    {
        FATAL("The size of the cropped image is bigger "
                      "than the image to crop");
    }

    //* Writing position map in png file *//
    io_png_write_f32(fname_positionMap, pos_map.img, 
                        pos_map.w, pos_map.h, pos_map.c);
    
    //* Writing synthesis map in png file *//
    io_png_write_f32(fname_synthesisMap, synth_map_crop.img, 
                        synth_map_crop.w, synth_map_crop.h, synth_map_crop.c);

    //* Writing output image in png file *//
    io_png_write_f32(fname_out, out_im_crop.img, 
                        out_im_crop.w, out_im_crop.h, out_im_crop.c);
                                              
    //* Release memory *//
    free(pos_map.img);
    free(synth_map.img);
    free(synth_map_crop.img);
    free(src_im.img);
    free(out_im.img);
    free(out_im_crop.img);
    
    return 0;
}
