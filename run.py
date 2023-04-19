#!/usr/bin/env python3

import subprocess
import argparse
import sys
import cv2

# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("ps", type=int)
ap.add_argument("t", type=float)
ap.add_argument("r", type=float)
ap.add_argument("o", type=float)
args = ap.parse_args()

image = cv2.imread('input_0.png')
cv2.imwrite('input_0.png', image)

p = subprocess.run(['quilting', '-p', str(args.ps), '-r', str(args.r), '-t', str(args.t), '-o',
                str(args.o), 'input_0.png', 'output.png', 'position_map.png', 'synthesis_map.png'])
if p.returncode != 0:
    with open("demo_failure.txt", "w") as file:
        file.write("Input image is too small. The size should be greater than 50x50 pixels")
        sys.exit(0)