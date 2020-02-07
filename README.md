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

Left-click points to which a curve should be fit. Bezier curves fit
best to 4 or 5 points; pressing space will create a new curve to fit;
the program will ensure they join nicely. Sliders give control over
the size and shape of the bins.

Press "h" to see help text.
