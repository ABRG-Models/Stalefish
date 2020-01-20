# Stalefish


compile on mac using:

g++ -o stalefish stalefish.cpp -i `pkg-config opencv --libs --cflags`

run using e.g.,:

./stalefish ../data/testimg1.tiff ../data/testimg2.tiff ../data/testimg3.tiff

left click points to draw segmented line, then

Press "x" to fit the curve
Press "w" to print the values
Press "n" to advance to next frame
Press "p" to increase the order of the polynomial fit
Press "b" to increase the number of bins for sampling
