/*
 * This is a utility to do the equivalent of opening an existing stalefish project and
 * writing it out to h5.
 */

#include "DM.h"
#include "FrameData.h"
#include <iostream>
#include <string>

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Opens, then saves a stalefish project. Two function calls:\n"
                  << "  DM::setup(argv[1]);\n"
                  << "  DM::writeFrames();\n";
        std::cout << "Please supply path to the stalefish hdf5 file" << std::endl;
        return 1;
    }
    std::string paramsfile (argv[1]);
    std::cout << "Opens, then saves a stalefish project. Two function calls:\n"
              << "  DM::setup("<<paramsfile<<");\n"
              << "  DM::writeFrames();\n";
    DM::i()->setup (paramsfile);
    DM::i()->writeFrames();
    return 0;
}