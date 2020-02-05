# Stalefish

## Dependencies

This program links to
[morphologica](https://github.com/ABRG-Models/morphologica) (and thus
to the associated libraries armadillo, opencv and hdf5). So head over
to that page and follow the instructions to build and install
morphologica.

## Compile Stalefish using cmake

```bash
mkdir build && cd build
cmake ..
make
```
## Running the program

run using e.g.,:
```
./path/to/stalefish ../data/testimg.json
```

The JSON file contains a list of the images to fit curves to, and some
other parameters.

left click points to which a curve should be fit. Bezier curves fit
best to 4 or 5 points; pressing space will create a new curve to fit;
the program will ensure they join nicely.

Press "x" to fit the curve(s)
Press "w" to save the values to an HDF5 file
Press "n" to advance to next frame
Press "1" to toggle view of Bezier handles
Press "2" to toggle view of user-provided points
Press "3" to toggle view of the fit
Press "4" to toggle view of the bins

Sliders give control over the size and shape of the bins.

Press "p" to increase the order of the polynomial fit
