/*
 * Stalefish - sfresave utility
 * Copyright (C) 2022 Sebastian S. James, Stuart P. Wilson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
