#!/bin/bash

# Take original icon and make many copies at different sizes. Assume
# original is 1024x1024

for i in 16 32 128 256 512; do
    convert sf_dafticon.png -resize ${i}x${i} icon_${i}x${i}.png
done

convert sf_dafticon.png -resize 32x32 icon_16x16@2x.png
convert sf_dafticon.png -resize 64x64 icon_32x32@2x.png
convert sf_dafticon.png -resize 256x256 icon_128x128@2x.png
convert sf_dafticon.png -resize 512x512 icon_256x256@2x.png
cp sf_dafticon.png icon_512x512@2x.png
