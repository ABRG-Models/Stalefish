/*
 * A 3D viewer for the data generated by Stalefish
 */
# include <GLFW/glfw3.h>
#include <morph/Visual.h>
using morph::Visual;
#include <morph/HdfData.h>
using morph::HdfData;
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>

using namespace std;

int main (int argc, char** argv)
{
    int rtn = -1;

    Visual v(1024, 768, "Visualization");
    v.zNear = 0.001;

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " /path/to/data.h5" << endl;
        return rtn;
    }

    string datafile (argv[1]);

    try {
        array<float, 3> offset = { 0.0, 0.0, 0.0 };
        array<float, 4> scale = { 0.1, 0.0, 1.0, 0.0};

        cout << "Opening H5 file " << datafile << endl;
        HdfData d(datafile, true); // true for read
        int nf = 0;
        //d.read_val ("/nframes", nf);

        vector<array<float, 12>> quads; // Get from HDF5
        vector<float> means;

        string frameName("");
        for (int i = 0; i<17/*nf*/; ++i) {

            stringstream ss;
            ss << "/Frame";
            ss.width(3);
            ss.fill('0');
            ss << i;
            frameName = ss.str();

            vector<array<float, 12>> frameQuads;
            vector<double> frameMeans;

            // Read quads and data for each frame and add to an overall pair of vectors...
            string str = frameName+"/sboxes";
            d.read_contained_vals (str.c_str(), frameQuads);
            str = frameName+"/means";
            d.read_contained_vals (str.c_str(), frameMeans);

            // Gah,convert frameMeans to float
            vector<float> frameMeansF;
            for (unsigned int j = 0; j < frameMeans.size(); ++j) {
                frameMeansF.push_back (static_cast<float>(frameMeans[i]));
            }
            quads.insert (quads.end(), frameQuads.begin(), frameQuads.end());
            means.insert (means.end(), frameMeansF.begin(), frameMeansF.end());
        }


        unsigned int visId = v.addQuadsVisual (&quads, offset, means, scale);
        cout << "Added Visual with visId " << visId << endl;
        v.render();

        while (v.readyToFinish == false) {
            glfwWaitEventsTimeout (0.018);
            v.render();
        }

    } catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << endl;
        rtn = -1;
    }


    return rtn;
}
