/*
 * A 3D viewer for the data generated by Stalefish
 */
# include <GLFW/glfw3.h>
#include <morph/Visual.h>
using morph::Visual;
#include <morph/HdfData.h>
using morph::HdfData;
#include <morph/Vector.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <morph/PointRowsVisual.h>

using namespace std;

int main (int argc, char** argv)
{
    int rtn = -1;

    Visual v(1024, 768, "Visualization");
    v.zNear = 0.001;
    v.zFar = 40.0;

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " /path/to/data.h5" << endl;
        return rtn;
    }

    string datafile (argv[1]);

    try {
        morph::Vector<float> offset = { 0.0, 0.0, 0.0 };

        morph::Scale<float> scale;
        scale.setParams (1.0, 0.0);

        cout << "Opening H5 file " << datafile << endl;
        HdfData d(datafile, true); // true for read
        int nf = 0;
        d.read_val ("/nframes", nf);

        // FIXME: thickness is hard-coded
        float thickness = 0.1f;
        float xx = thickness;

        vector<array<float, 12>> quads; // Get from HDF5
        vector<array<float, 12>> fquads; // Flat quads, for the flat visualization
        vector<morph::Vector<float>> points; // Centres of boxes; for smooth surface (points rows)
        vector<float> means;
        vector<float> fmeans;

        string frameName("");
        for (int i = 1; i<=nf; ++i) {

            stringstream ss;
            ss << "/Frame";
            ss.width(3);
            ss.fill('0');
            ss << i;
            frameName = ss.str();

            vector<array<float, 12>> frameQuads;
            vector<double> frameMeans;
            vector<morph::Vector<float>> framePoints;

            // Read quads and data for each frame and add to an overall pair of vectors...
            string str = frameName+"/sboxes";
            d.read_contained_vals (str.c_str(), frameQuads);

            for (auto fq : frameQuads) {
                // FIXME: Use centre of box, or even each end of box, or something
                morph::Vector<float> pt = {fq[0],fq[1],fq[2]};
                framePoints.push_back (pt);
            }

            bool autoscale_per_slice = true;
            if (autoscale_per_slice) {
                // Use the auto-scaled version of the means
                str = frameName+"/means_autoscaled";
                d.read_contained_vals (str.c_str(), frameMeans);
            } else {
                // Use the raw means and autoscale them as an entire group
                str = frameName+"/means";
                d.read_contained_vals (str.c_str(), frameMeans);
                scale.do_autoscale = true;
            }

            // Gah, convert frameMeans to float (there's a better way to do this)
            vector<float> frameMeansF;
            for (unsigned int j = 0; j < frameMeans.size(); ++j) {
                frameMeansF.push_back (static_cast<float>(frameMeans[j]));
            }

            quads.insert (quads.end(), frameQuads.begin(), frameQuads.end());
            means.insert (means.end(), frameMeansF.begin(), frameMeansF.end());
            points.insert (points.end(), framePoints.begin(), framePoints.end());

            // Load in linear stuff as well, to make up flat boxes? Or easier to do at source?
            vector<float> linbins;
            str = frameName+"/sbox_linear_distance";
            d.read_contained_vals (str.c_str(), linbins);

            vector<array<float,12>> flatsurf_boxes;
            array<float, 12> sbox;
            for (unsigned int j = 1; j < linbins.size(); ++j) {
                // c1 x,y,z
                sbox[0] = xx-thickness;
                sbox[1] = linbins[j-1];  // y
                sbox[2] = 0.0;
                // c2 x,y,z
                sbox[3] = xx-thickness;               // x
                sbox[4] = linbins[j];  // y
                sbox[5] = 0.0;
                // c3 x,y,z
                sbox[6] = xx;
                sbox[7] = linbins[j];     // y
                sbox[8] = 0.0;
                // c4 x,y,z
                sbox[9] = xx;
                sbox[10] = linbins[j-1];  // y
                sbox[11] = 0.0;

                flatsurf_boxes.push_back (sbox);
            }

            fquads.insert (fquads.end(), flatsurf_boxes.begin(), flatsurf_boxes.end());
            fmeans.insert (fmeans.end(), frameMeansF.begin(), --frameMeansF.end());
            xx += thickness;
        }

        unsigned int visId = 0;
        //offset[0] -= 3.0;
        //visId = v.addQuadsVisual (&quads, offset, means, scale);
        //cout << "Added Visual with visId " << visId << endl;

        //offset[0]+=3.0;
        //cout << "fquads size: " << fquads.size() << "fmeans isze: " << fmeans.size() << endl;
        //visId = v.addQuadsVisual (&fquads, offset, fmeans, scale);
        //cout << "Added Visual with visId " << visId << endl;

        //offset[0]+=5.0;
        //visId = v.addPointRowsVisual (&points, offset, means, scale, morph::ColourMapType::Plasma);
        visId = v.addVisualModel (new morph::PointRowsVisual<float> (v.shaderprog,
                                                                     &points, offset,
                                                                     &means, scale,
                                                                     morph::ColourMapType::Plasma));

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
