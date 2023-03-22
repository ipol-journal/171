#!/usr/bin/env python3

import subprocess
import argparse

# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("ps", type=float)
ap.add_argument("t", type=float)
ap.add_argument("r", type=float)
ap.add_argument("o", type=float)

args = ap.parse_args()

p = subprocess.run(['quilting', '-p', str(args.ps), '-r', str(args.r), '-t', str(args.t), '-o',
                    str(args.o), 'input_0.png', 'output.png', 'position_map.png', 'synthesis_map.png']) 
