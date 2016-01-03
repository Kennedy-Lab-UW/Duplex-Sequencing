#!/usr/bin/env python
import os
import pickle
from argparse import ArgumentParser
import shutil
import sys

parser = ArgumentParser()
parser.add_argument('--scripts_folder', dest='scripts_folder', help='Duplex-Sequencing scripts folder', required=True)
parser.add_argument('--output_folder', dest='output_folder', help='Results folder', required=True)
o = parser.parse_args()

scripts=os.path.realpath(o.scripts_folder)

if not os.path.exists(o.output_folder):
    os.makedirs(o.output_folder)
output=os.path.realpath(o.output_folder)

os.chdir(scripts)

new_folder=o.output_folder.split("/")
new_folder=new_folder[len(new_folder)-1]

new_files=os.listdir(scripts)
new_files=[ s for s in new_files if new_folder in s]

for s in new_files:
    if s != new_folder:
        shutil.move(s, output+"/"+s)

sys.exit(0)
