// The Data Manager/drawing class
#include "DM.h"
// The image-and-associated-fit-parameters class
#include "FrameData.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Please supply path to the stalefish hdf5 file" << std::endl;
        return 1;
    }
    std::string paramsfile (argv[1]);
    DM::i()->setup (paramsfile);
    DM::i()->writeFrames();
    return 0;
}
