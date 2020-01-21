# Stalefish

## Compile using cmake

```bash
mkdir build && cd build
cmake ..
make
```

Alternatively, compile on Mac/Linux using a single command (for now):

```bash
cd src
g++ -o stalefish stalefish.cpp -i `pkg-config opencv --libs --cflags`
```

## Running the program

run using e.g.,:
```
./path/to/stalefish ../data/testimg1.tiff ../data/testimg2.tiff ../data/testimg3.tiff
```
left click points to draw segmented line, then

Press "x" to fit the curve
Press "w" to print the values
Press "n" to advance to next frame
Press "p" to increase the order of the polynomial fit
Press "b" to increase the number of bins for sampling
