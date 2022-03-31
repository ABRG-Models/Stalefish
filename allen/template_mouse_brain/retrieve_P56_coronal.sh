#!/bin.sh

# Get Allen Mouse Brain atlas 602630314 with downsample 3.
python retrieve_template.py 602630314 3

# Get the annotations.
python retrieve_template.py 602630314 3 1
