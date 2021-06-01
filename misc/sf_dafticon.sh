#!/bin/bash

# Take original icon and make many copies at different sizes. Assume
# original is 1024x1024

for i in 16 32 128 256 512; do
    convert sf_dafticon.png -resize ${i}x${i} sf_dafticon${i}.png
done
