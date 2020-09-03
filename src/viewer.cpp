/*
 * A 3D viewer for the data generated by Stalefish
 */
#include <GLFW/glfw3.h>
#include <morph/Visual.h>
#include <morph/HdfData.h>
#include <morph/Vector.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <morph/PointRowsVisual.h>
#include <morph/ScatterVisual.h>

using namespace std;

int main (int argc, char** argv)
{
    int rtn = -1;

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " /path/to/data.h5 [perslice|overall*] [landmark*|auto]\n"
             << "  where 'perslice' autoscales the signal for each slice before\n"
             << "  assembling them together, and 'overall' autoscales the signal AFTER the\n"
             << "  slices have been assembled together. providing a third arg of 'auto' will\n"
             << "  show the curve-autoaligned' alignment, even if landmark alignment has been\n"
             << "  carried out.\n";
        return rtn;
    }

    string datafile (argv[1]);

    morph::Visual v(1024, 768, datafile);
    v.zNear = 0.001;
    v.zFar = 40.0;
    v.setZDefault (-15.4);

    bool autoscale_per_slice = false;
    if (argc > 2) {
        string autoscale_str (argv[2]);
        if (autoscale_str == "perslice") {
            autoscale_per_slice = true;
        } else {
            autoscale_per_slice = false;
        }
    }

    bool align_lm = true;
    if (argc > 3) {
        string align_str (argv[3]);
        if (align_str == "auto") {
            align_lm = false;
        } else {
            align_lm = true;
        }
    }

    try {
        morph::Vector<float> offset = { 0.0, 0.0, 0.0 };

        morph::Scale<float> scale;
        scale.setParams (1.0, 0.0);

        // FIXME: thickness is hard-coded
        float thickness = 0.1f;
        float xx = thickness;

        vector<array<float, 12>> quads_autoaligned; // Get from HDF5
        vector<array<float, 12>> quads_lmaligned;
        vector<array<float, 12>> quads_scaled;
        vector<array<float, 12>> fquads; // Flat quads, for the flat visualization
        vector<morph::Vector<float>> points_autoaligned; // Centres of boxes; for smooth surface (points rows)
        vector<morph::Vector<float>> points_lmaligned; // Centres of boxes; for smooth surface (points rows)
        vector<morph::Vector<float>> points_scaled; // Centres of boxes; for smooth surface (points rows)

        vector<morph::Vector<float>> landmarks_autoaligned;
        vector<morph::Vector<float>> landmarks_lmaligned;
        vector<float> landmarks_id;

        vector<float> means;
        vector<float> fmeans;

        {
            cout << "Opening H5 file " << datafile << endl;
            morph::HdfData d(datafile, true); // true for read
            int nf = 0;
            d.read_val ("/nframes", nf);

            // Check first frame for alignments
            bool lmalignComputed = false;
            d.read_val ("Frame001/lmalign/computed", lmalignComputed);
            bool autoalignComputed = false;
            d.read_val ("Frame001/autoalign/computed", autoalignComputed);

            string frameName("");
            for (int i = 1; i<=nf; ++i) {

                stringstream ss;
                ss << "/Frame";
                ss.width(3);
                ss.fill('0');
                ss << i;
                frameName = ss.str();

                vector<array<float, 12>> frameQuads_scaled;
                vector<array<float, 12>> frameQuads_lmaligned;
                vector<array<float, 12>> frameQuads_autoaligned;

                vector<morph::Vector<float>> framePoints_autoaligned;
                vector<morph::Vector<float>> framePoints_lmaligned;
                vector<morph::Vector<float>> framePoints_scaled;

                // Read quads and data for each frame and add to an overall pair of vectors...
                string str = frameName+"/autoalign/sboxes";
                d.read_contained_vals (str.c_str(), frameQuads_autoaligned);
                str = frameName+"/lmalign/sboxes";
                d.read_contained_vals (str.c_str(), frameQuads_lmaligned);
                // Un-transformed:
                str = frameName+"/scaled/sboxes";
                d.read_contained_vals (str.c_str(), frameQuads_scaled);

                for (auto fq : frameQuads_autoaligned) {
                    // FIXME: Use centre of box, or even each end of box, or something
                    morph::Vector<float> pt = {fq[0],fq[1],fq[2]};
                    framePoints_autoaligned.push_back (pt);
                }
                for (auto fq : frameQuads_lmaligned) {
                    morph::Vector<float> pt = {fq[0],fq[1],fq[2]};
                    framePoints_lmaligned.push_back (pt);
                }
                for (auto fq : frameQuads_scaled) {
                    morph::Vector<float> pt = {fq[0],fq[1],fq[2]};
                    framePoints_scaled.push_back (pt);
                }

                // Landmarks
                vector<array<float, 3>> LM_autoaligned;
                vector<array<float, 3>> LM_lmaligned;
                str = frameName+"/autoalign/landmarks";
                d.read_contained_vals (str.c_str(), LM_autoaligned);
                str = frameName+"/lmalign/landmarks";
                d.read_contained_vals (str.c_str(), LM_lmaligned);
                size_t lmcount = 0;
                float lmid = 0.0f;
                float lmidmax = (float)LM_autoaligned.size();
                for (auto lm : LM_autoaligned) {
                    morph::Vector<float> _lm = {lm[0],lm[1],lm[2]};
                    landmarks_autoaligned.push_back (_lm);
                    lmid = (float)lmcount++ / lmidmax;
                    landmarks_id.push_back (lmid);
                }
                for (auto lm : LM_lmaligned) {
                    morph::Vector<float> _lm = {lm[0],lm[1],lm[2]};
                    landmarks_lmaligned.push_back (_lm);
                }

                vector<double> frameMeans;
                if (autoscale_per_slice) {
                    // Use the auto-scaled version of the means, with each slice autoscaled to [0,1]
                    str = frameName+"/signal/postproc/boxes/means_autoscaled";
                    d.read_contained_vals (str.c_str(), frameMeans);
                } else {
                    // Use the raw means and autoscale them as an entire group
                    str = frameName+"/signal/postproc/boxes/means";
                    d.read_contained_vals (str.c_str(), frameMeans);
                    // The morph::Scale object scale with autoscale the who thing.
                    scale.do_autoscale = true;
                }

                // Gah, convert frameMeans to float (there's a better way to do this)
                vector<float> frameMeansF;
                for (unsigned int j = 0; j < frameMeans.size(); ++j) {
                    frameMeansF.push_back (static_cast<float>(frameMeans[j]));
                }

                quads_autoaligned.insert (quads_autoaligned.end(), frameQuads_autoaligned.begin(), frameQuads_autoaligned.end());
                quads_lmaligned.insert (quads_lmaligned.end(), frameQuads_lmaligned.begin(), frameQuads_lmaligned.end());
                quads_scaled.insert (quads_scaled.end(), frameQuads_scaled.begin(), frameQuads_scaled.end());
                means.insert (means.end(), frameMeansF.begin(), frameMeansF.end());

                points_lmaligned.insert (points_lmaligned.end(), framePoints_lmaligned.begin(), framePoints_lmaligned.end());
                points_autoaligned.insert (points_autoaligned.end(), framePoints_autoaligned.begin(), framePoints_autoaligned.end());
                points_scaled.insert (points_scaled.end(), framePoints_scaled.begin(), framePoints_scaled.end());

                // Load in linear stuff as well, to make up flat boxes? Or easier to do at source?
                vector<float> linbins;
                //str = frameName+"/scaled/flattened/sbox_linear_distance";
                // Better:
                //str = frameName+"/lmalign/flattened/sbox_angles";
                // Even better:
                if (lmalignComputed == true && align_lm == true) {
                    str = frameName+"/lmalign/flattened/sbox_linear_distance";
                } else {
                    str = frameName+"/autoalign/flattened/sbox_linear_distance";
                }
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
            cout << "fquads size: " << fquads.size() << ", fmeans size: " << fmeans.size() << endl;
            cout << "fmeans min/max: " << fmeans.size() << endl;

            cout << "landmarks_autoaligned.size(): " << landmarks_autoaligned.size() << endl;

            unsigned int visId = 0;

#if 0
            // The 'ribbons' maps
            offset[0] = 0.0;
            visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                     &quads_autoaligned, offset,
                                                                     &means, scale,
                                                                     morph::ColourMapType::MonochromeBlue));
            offset[0] = 5.0;
            visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                     &quads_lmaligned, offset,
                                                                     &means, scale,
                                                                     morph::ColourMapType::MonochromeRed));
            offset[0] = 10.0;
            visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                     &quads_scaled, offset,
                                                                     &means, scale,
                                                                     morph::ColourMapType::MonochromeGreen));
#endif

            // This is the flattened map; showing it alongside the 3D map for now
            offset[0]=-5.0;
            visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                     &fquads, offset,
                                                                     &fmeans, scale,
                                                                     morph::ColourMapType::Greyscale));

            offset[0]=0.0;

            // Show landmark aligned for preference:
            if (lmalignComputed == true && align_lm == true) {
                visId = v.addVisualModel (new morph::PointRowsVisual<float> (v.shaderprog,
                                                                             &points_lmaligned, offset,
                                                                             &means, scale,
                                                                             morph::ColourMapType::MonochromeRed));
                // Show the landmarks with a ScatterVisual
                visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                           &landmarks_lmaligned, offset,
                                                                           &landmarks_id, 0.07f, scale,
                                                                           morph::ColourMapType::Plasma));
            } else {
                visId = v.addVisualModel (new morph::PointRowsVisual<float> (v.shaderprog,
                                                                             &points_autoaligned, offset,
                                                                             &means, scale,
                                                                             morph::ColourMapType::MonochromeBlue));
                visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                           &landmarks_autoaligned, offset,
                                                                           &landmarks_id, 0.07f, scale,
                                                                           morph::ColourMapType::Plasma));
            }

            cout << "Added Visual with visId " << visId << endl;
        }
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
