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

    // Though DM::setup will test filename, lets do so explicitly here as want to exit if file is .json.
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

    // Set up from .h5
    std::string jsoncontent("");
    try {
        morph::HdfData d(paramsfile, true); // true for read
        d.read_string ("/config", jsoncontent);
        std::string jsonfile(paramsfile);
        // Strip off .h5
        morph::Tools::stripFileSuffix (jsonfile);
        jsonfile.append(".json");
        // Check if jsonfile already exists & bail
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
        std::ofstream jf;
        jf.open (jsonfile.c_str(), std::ios::out|std::ios::trunc);
        if (!jf.is_open()) {
            std::stringstream ee;
            ee << "Failed to open json file '" << jsonfile << "' for writing.";
            throw std::runtime_error (ee.str());
        }
        jf << jsoncontent;
        jf.close();

        std::cout << "JSON content written to file '" << jsonfile <<"'\n";

    } catch (const std::exception& e) {
        std::cerr << "Error reading h5/writing to json: " << e.what() << std::endl;
        exit (1);
    }

    return 0;
}
