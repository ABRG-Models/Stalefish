/*
 * Stalefish - sfgetjson utility
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
 * This is a utility to open an existing .h5 project and save its stored json into the
 * json file.
 */

#include <iostream>
#include <string>
#include <morph/tools.h>
#include <morph/HdfData.h>

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Opens a stalefish project and extracts the json config,\n"
                  << "writing it into a companion file. The filename matches the original\n"
                  << "file, but with the .h5 suffix changed to .json.\n"
                  << "Please supply a path to the stalefish hdf5 file" << std::endl;
        return 1;
    }
    std::string paramsfile (argv[1]);

    // First check that the user-supplied paramsfile looks like it's an h5 file
    std::string::size_type jsonpos = paramsfile.find(".json");
    if (jsonpos != std::string::npos) {
        std::cerr << paramsfile << " looks like it IS a json file; exiting." << std::endl;
        return -1;
    }

    std::string::size_type h5pos = paramsfile.find(".h5");
    if (h5pos == std::string::npos) {
        std::cerr << paramsfile << " seems not to be an h5 file; exiting." << std::endl;
        return -1;
    }

    std::string jsoncontent("");
    try {
        // Open the .h5 file with a morph::HdfData object
        morph::HdfData d(paramsfile, true); // true for read
        // Read out the json config
        d.read_string ("/config", jsoncontent);
        // Prepare a file for output
        std::string jsonfile(paramsfile);
        // Strip off .h5 suffix and append .json instead
        morph::Tools::stripFileSuffix (jsonfile);
        jsonfile.append(".json");
        // Check if the json file already exists & if so, bail out
        if (morph::Tools::dirExists (jsonfile)) {
            std::stringstream ee;
            ee << jsonfile << " exists as a directory. Won't replace it.";
            throw std::runtime_error (ee.str());
        }
        if (morph::Tools::fileExists (jsonfile)) {
            std::stringstream ee;
            ee << jsonfile << " already exists. Won't replace it.";
            throw std::runtime_error (ee.str());
        }
        // jsonfile path looks good, so write out the config
        std::ofstream jf;
        jf.open (jsonfile.c_str(), std::ios::out|std::ios::trunc);
        if (!jf.is_open()) {
            std::stringstream ee;
            ee << "Failed to open json file '" << jsonfile << "' for writing.";
            throw std::runtime_error (ee.str());
        }
        jf << jsoncontent;
        jf.close();
        // success :)
        std::cout << "JSON content written to file '" << jsonfile <<"'\n";

    } catch (const std::exception& e) {
        std::cerr << "Error reading h5/writing to json: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}
