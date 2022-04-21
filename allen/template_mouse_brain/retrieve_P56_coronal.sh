#!/bin.sh

# Get Allen Mouse Brain atlas 602630314 (structural MRI data) with downsample 3.
python retrieve_template.py 602630314 3

# Get the annotations - the coloured map of locations.
python retrieve_template.py 602630314 3 1
