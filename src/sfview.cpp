/*
 * A 3D viewer for the data generated by Stalefish
 */
#include <morph/Visual.h>
#include <morph/HdfData.h>
#include <morph/Vector.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <morph/PointRowsVisual.h>
#include <morph/ScatterVisual.h>
#include <morph/QuadsVisual.h>
#include <morph/RodVisual.h>
#include <popt.h>

using namespace std;

// Derive Visual to add the extra sfview-specific keyhandling callback
class SFVisual : public morph::Visual
{
public:
    SFVisual (int width, int height, const std::string& title,
              const morph::Vector<float> caOffset, const morph::Vector<float> caLength, const float caThickness)
        : morph::Visual (width, height, title, caOffset, caLength, caThickness) {}

    //! Vector of VisualModel IDs for the landmarks. To hide landmarks, hide these.
    std::vector<unsigned int> landmarks;
    //! The 'angle=0' lines - in Jet colours
    std::vector<unsigned int> angle_centres;
    //! The axis (or axes) from 'axis marks'
    std::vector<unsigned int> axes;
    //! The 3D surfaces
    std::vector<unsigned int> surfaces_3d;
    //! 2D surfaces
    std::vector<unsigned int> surfaces_2d;

protected:
    //! Act on keys and toggle 'hidden' for the relevant VisualModels
    virtual void key_callback_extra (GLFWwindow* window, int key, int scancode, int action, int mods)
    {
        // Landmarks
        if (key == GLFW_KEY_F && action == GLFW_PRESS) {
            for (auto id : this->landmarks) { this->vm[id]->toggleHide(); }
        }
        // Angle markers
        if (key == GLFW_KEY_G && action == GLFW_PRESS) {
            for (auto id : this->angle_centres) { this->vm[id]->toggleHide(); }
        }
        // axes
        if (key == GLFW_KEY_D && action == GLFW_PRESS) {
            for (auto id : this->axes) { this->vm[id]->toggleHide(); }
        }
        // 2D map
        if (key == GLFW_KEY_J && action == GLFW_PRESS) {
            for (auto id : this->surfaces_2d) { this->vm[id]->toggleHide(); }
        }
        // 3D map
        if (key == GLFW_KEY_K && action == GLFW_PRESS) {
            for (auto id : this->surfaces_3d) { this->vm[id]->toggleHide(); }
        }

        // Additional help
        if (key == GLFW_KEY_H && action == GLFW_PRESS) {
            std::cout << "sfview help:\n";
            std::cout << "f: toggle (show/hide) landmarks\n";
            std::cout << "g: toggle zero angle marks\n";
            std::cout << "d: toggle user-defined brain axis\n";
            std::cout << "j: toggle 2D brain map\n";
            std::cout << "k: toggle 3D brain surface\n";
        }
    }
};

//! libpopt features - the features that are available to change on the command line.
struct CmdOptions
{
    //! If true, then autoscale the signal for each slice
    int scale_perslice;
    //! If true, then use auto align, even if landmark alignment data is present
    int use_autoalign;
    //! If true, plot ribbons for 3D instead of the smooth surface map. Useful for debugging
    int show_ribbons;
    //! Show landmarks from first model
    int show_landmarks;
    //! Show landmarks from all models
    int show_landmarks_all;
    //! Takes the index of the flattened map to show
    int show_flattened;
    //! flattened_type could be 0. aligned linear distance, starting from angle 0 1. linear distance, centered, 2. angle of thing about 0.
    int flattened_type;
    //! Temporary datafile, for adding to datafiles.
    char* datafile;
    //! The h5 files to visualize
    std::vector<string> datafiles;
};

//! Initialise a CmdOptions object
void zeroCmdOptions (CmdOptions* copts)
{
    copts->scale_perslice = 0;
    copts->use_autoalign = 0;
    copts->show_ribbons = 0;
    copts->show_landmarks = 1;
    copts->show_landmarks_all = 0;
    copts->show_flattened = 0;
    copts->flattened_type = 0;
    copts->datafile = (char*)0;
    copts->datafiles.clear();
}

//! cmdOptions is global, allowing callbacks to access this easily.
struct CmdOptions cmdOptions;

/*!
 * This callback is used when there's a -f option, to allow me to
 * collect multiple files for visualization (e.g. to vis. multiple layers)
 *
 * e.g. viewer -f /path/to/ctx_superficial.h5 -f /path/to/ctx_mid.h5
 */
void popt_option_callback (poptContext con,
                           enum poptCallbackReason reason,
                           const struct poptOption * opt,
                           const char * arg,
                           void * data)
{
    switch(reason) {
    case POPT_CALLBACK_REASON_PRE: { break; } // Doesn't occur
    case POPT_CALLBACK_REASON_POST: { break; } // Ignore
    case POPT_CALLBACK_REASON_OPTION:
    {
        // Test shortName. This means we could respond to other "multiple options"
        if (opt->shortName == 'f') {
            cmdOptions.datafiles.push_back (cmdOptions.datafile);
        }
        break;
    }
    }
}

//! Add just the landmarks in the datafile
int addLandmarks (SFVisual& v, const string& datafile, const CmdOptions& co)
{
    int rtn = 0;

    bool align_lm = co.use_autoalign > 0 ? false : true;

    try {
        morph::Vector<float> offset = { 0.0, 0.0, 0.0 };

        morph::Scale<float> scale;
        scale.setParams (1.0, 0.0);
        float xx = 0.0f;

        vector<morph::Vector<float>> landmarks_autoaligned;
        vector<morph::Vector<float>> landmarks_lmaligned;

        vector<float> landmarks_id;

        {
            morph::HdfData d(datafile, true); // true for read
            int nf = 0;
            d.read_val ("/nframes", nf);

            // Check first frame for alignments
            bool lmalignComputed = false;
            d.read_val ("/Frame001/lmalign/computed", lmalignComputed);
            bool autoalignComputed = false;
            d.read_val ("/Frame001/autoalign/computed", autoalignComputed);

            string frameName("");
            for (int i = 1; i<=nf; ++i) {

                stringstream ss;
                ss << "/Frame";
                ss.width(3);
                ss.fill('0');
                ss << i;
                frameName = ss.str();

                // x position comes from FrameNNN/class/layer_x
                string str = frameName+"/class/layer_x";
                d.read_val (str.c_str(), xx);

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
            }
            unsigned int visId = 0;

            offset[0]=0.0;

            // Show landmark aligned for preference:
            if (lmalignComputed == true && align_lm == true) {
                // Show the landmarks with a ScatterVisual
                visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                           &landmarks_lmaligned, offset,
                                                                           &landmarks_id, 0.07f, scale,
                                                                           morph::ColourMapType::Plasma));
            } else {
                visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                           &landmarks_autoaligned, offset,
                                                                           &landmarks_id, 0.07f, scale,
                                                                           morph::ColourMapType::Plasma));
            }
            v.landmarks.push_back (visId);
        }
    } catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << endl;
        rtn = -1;
    }

    return rtn;
}

//! Add a visual model for the expression surface, created from the file datafile, to the scene v
int addVisMod (SFVisual& v, const string& datafile, const CmdOptions& co, const float hue)
{
    int rtn = 0;

    bool autoscale_per_slice = co.scale_perslice > 0 ? true : false;
    bool align_lm = co.use_autoalign > 0 ? false : true;
    bool showribbons = co.show_ribbons > 0 ? true : false;

    try {
        morph::Vector<float> offset = { 0.0, 0.0, 0.0 };

        morph::Scale<float> scale;
        scale.setParams (1.0, 0.0);

        float thickness = 0.0f;
        float xx = 0.0f;

        vector<array<float, 12>> quads_autoaligned; // Get from HDF5
        vector<array<float, 12>> quads_lmaligned;
        vector<array<float, 12>> quads_scaled;
        vector<morph::Vector<float>> points_autoaligned; // Centres of boxes; for smooth surface (points rows)
        vector<morph::Vector<float>> points_lmaligned; // Centres of boxes; for smooth surface (points rows)
        vector<morph::Vector<float>> points_scaled; // Centres of boxes; for smooth surface (points rows)
        vector<float> means;

        vector<morph::Vector<float>> centres_lmaligned;
        vector<morph::Vector<float>> centres_autoaligned;
        vector<morph::Vector<float>> AM_origins_lmaligned;
        vector<morph::Vector<float>> AM_origins_autoaligned;
        vector<float> centres_id;

        {
            cout << "Opening H5 file " << datafile << endl;
            morph::HdfData d(datafile, true); // true for read
            d.read_error_action = morph::ReadErrorAction::Exception;
            int nf = 0;
            d.read_val ("/nframes", nf);

            // Check first frame for alignments
            bool lmalignComputed = false;
            d.read_val ("/Frame001/lmalign/computed", lmalignComputed);
            bool autoalignComputed = false;
            d.read_val ("/Frame001/autoalign/computed", autoalignComputed);

            string frameName("");
            for (int i = 1; i<=nf; ++i) {

                stringstream ss;
                ss << "/Frame";
                ss.width(3);
                ss.fill('0');
                ss << i;
                frameName = ss.str();

                // x position comes from FrameNNN/class/layer_x
                string str = frameName+"/class/layer_x";
                d.read_val (str.c_str(), xx);
                str = frameName+"/class/thickness";
                d.read_val (str.c_str(), thickness);

                // Centres. Get the index into the fitted_lmaligned, to get y/z coordinates of centre locations
                str = frameName+"/lmalign/centre_box_index";
                int clm_idx = 0;
                d.read_val (str.c_str(), clm_idx);

                str = frameName+"/autoalign/centre_box_index";
                int caa_idx = 0;
                d.read_val (str.c_str(), caa_idx);

                std::vector<cv::Point2d> fitted_lmaligned;
                str = frameName+"/lmalign/fitted";
                d.read_contained_vals (str.c_str(), fitted_lmaligned);

                std::vector<cv::Point2d> fitted_autoaligned;
                str = frameName+"/autoalign/fitted";
                d.read_contained_vals (str.c_str(), fitted_autoaligned);

                morph::Vector<float> cp;
                cp[0] = xx;
                cp[1] = fitted_lmaligned[clm_idx].x;
                cp[2] = fitted_lmaligned[clm_idx].y;
                centres_lmaligned.push_back (cp);
                cp[1] = fitted_autoaligned[caa_idx].x;
                cp[2] = fitted_autoaligned[caa_idx].y;
                centres_autoaligned.push_back (cp);
                centres_id.push_back (0.1f*(float)i);

                // axis alignment marks are optional, and the data may not exist in the h5 file.
                try {
                    std::vector<cv::Point2d> AM_lmaligned;
                    str = frameName+"/lmalign/alignmark_origins";
                    d.read_contained_vals (str.c_str(), AM_lmaligned);
                    cp[1] = AM_lmaligned[0].x;
                    cp[2] = AM_lmaligned[0].y;
                    AM_origins_lmaligned.push_back (cp);
                } catch (const exception& ee) {
                    // Ignore missing AM_origins
                    std::cout << "Missing alignmark_origins: " << ee.what() << std::endl;
                }

                // axis alignment marks are optional, and the data may not exist in the h5 file.
                try {
                    std::vector<cv::Point2d> AM_autoaligned;
                    str = frameName+"/autoalign/alignmark_origins";
                    d.read_contained_vals (str.c_str(), AM_autoaligned);
                    cp[1] = AM_autoaligned[0].x;
                    cp[2] = AM_autoaligned[0].y;
                    AM_origins_autoaligned.push_back (cp);
                } catch (const exception& ee) {
                    // Ignore missing AM_origins
                    std::cout << "Missing alignmark_origins: " << ee.what() << std::endl;
                }

                vector<array<float, 12>> frameQuads_scaled;
                vector<array<float, 12>> frameQuads_lmaligned;
                vector<array<float, 12>> frameQuads_autoaligned;

                vector<morph::Vector<float>> framePoints_autoaligned;
                vector<morph::Vector<float>> framePoints_lmaligned;
                vector<morph::Vector<float>> framePoints_scaled;

                // Read quads and data for each frame and add to an overall pair of vectors...
                str = frameName+"/autoalign/sboxes";
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
            }
            unsigned int visId = 0;

            offset[0]=0.0;

            // Show landmark aligned for preference:
            if (lmalignComputed == true && align_lm == true) {
                if (showribbons) {
                    visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                             &quads_lmaligned, offset,
                                                                             &means, scale,
                                                                             morph::ColourMapType::Monochrome, hue));
                } else {
                    visId = v.addVisualModel (new morph::PointRowsVisual<float> (v.shaderprog,
                                                                                 &points_lmaligned, offset,
                                                                                 &means, scale,
                                                                                 morph::ColourMapType::Monochrome, hue));
                }
                v.surfaces_3d.push_back (visId);

                visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                           &centres_lmaligned, offset,
                                                                           &centres_id, 0.03f, scale,
                                                                           morph::ColourMapType::Jet));
                v.angle_centres.push_back (visId);

                if (!AM_origins_lmaligned.empty()) {
                    size_t amo_last = AM_origins_lmaligned.size()-1;
                    std::array<float, 3> rcol = {1.0f,1.0f,1.0f};
                    visId = v.addVisualModel (new morph::RodVisual (v.shaderprog, offset,
                                                                    AM_origins_lmaligned[0], AM_origins_lmaligned[amo_last],
                                                                    0.05f, rcol));
                    v.axes.push_back (visId);
                }

            } else {

                if (showribbons) {
                    visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                             &quads_autoaligned, offset,
                                                                             &means, scale,
                                                                             morph::ColourMapType::Monochrome, hue));
                } else {
                    visId = v.addVisualModel (new morph::PointRowsVisual<float> (v.shaderprog,
                                                                                 &points_autoaligned, offset,
                                                                                 &means, scale,
                                                                                 morph::ColourMapType::Monochrome, hue));
                }
                v.surfaces_3d.push_back (visId);

                visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                           &centres_autoaligned, offset,
                                                                           &centres_id, 0.03f, scale,
                                                                           morph::ColourMapType::Jet));
                v.angle_centres.push_back (visId);

                if (!AM_origins_autoaligned.empty()) {
                    size_t amo_last = AM_origins_autoaligned.size()-1;
                    std::array<float, 3> rcol = {1.0f,1.0f,1.0f};
                    visId = v.addVisualModel (new morph::RodVisual (v.shaderprog, offset,
                                                                    AM_origins_autoaligned[0], AM_origins_autoaligned[amo_last],
                                                                    0.05f, rcol));
                    v.axes.push_back (visId);
                }
            }
        }
    } catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << endl;
        rtn = -1;
    }

    return rtn;
}

//! Add flattened map
int addFlattened (SFVisual& v, const string& datafile, const CmdOptions& co)
{
    int rtn = 0;

    bool autoscale_per_slice = co.scale_perslice > 0 ? true : false;
    bool align_lm = co.use_autoalign > 0 ? false : true;

    try {
        morph::Vector<float> offset = { 0.0, 0.0, 0.0 };

        morph::Scale<float> scale;
        scale.setParams (1.0, 0.0);

        float thickness = 0.0f;
        float xx = 0.0f;
        vector<array<float, 12>> fquads; // Flat quads, for the flat visualization
        vector<float> means;
        vector<float> fmeans;

        {
            cout << "Opening H5 file " << datafile << endl;
            morph::HdfData d(datafile, true); // true for read
            int nf = 0;
            d.read_val ("/nframes", nf);

            // Check first frame for alignments
            bool lmalignComputed = false;
            d.read_val ("/Frame001/lmalign/computed", lmalignComputed);
            bool autoalignComputed = false;
            d.read_val ("/Frame001/autoalign/computed", autoalignComputed);

            if (co.flattened_type == 1) {
                // linear distance boxes:
                std::cout << "Showing centered, linear distance as y axis for 2D map.\n";
            } else if (co.flattened_type == 2) {
                // angle based boxes:
                std::cout << "Showing angle as y axis for 2D map.\n";
            } else {
                // linear distance based on 0 angle starting point:
                if (lmalignComputed == true && align_lm == true) {
                    std::cout << "Showing linear distance unwrapped from 3D landmark aligned reconstruction as y axis for 2D map.\n";
                } else {
                    std::cout << "Showing linear distance unwrapped from 3D auto-aligned reconstruction as y axis for 2D map.\n";
                }
            }

            string frameName("");
            for (int i = 1; i<=nf; ++i) {

                stringstream ss;
                ss << "/Frame";
                ss.width(3);
                ss.fill('0');
                ss << i;
                frameName = ss.str();

                // x position comes from FrameNNN/class/layer_x
                string str = frameName+"/class/layer_x";
                d.read_val (str.c_str(), xx);
                str = frameName+"/class/thickness";
                d.read_val (str.c_str(), thickness);
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
                // Load in linear stuff as well, to make up flat boxes? Or easier to do at source?
                vector<float> linbins;
                if (co.flattened_type == 1) {
                    // linear distance boxes:
                    str = frameName+"/scaled/flattened/sbox_linear_distance";
                } else if (co.flattened_type == 2) {
                    // angle based boxes:
                    str = frameName+"/lmalign/flattened/sbox_angles";
                } else {
                    // linear distance based on 0 angle starting point:
                    if (lmalignComputed == true && align_lm == true) {
                        str = frameName+"/lmalign/flattened/sbox_linear_distance";
                    } else {
                        str = frameName+"/autoalign/flattened/sbox_linear_distance";
                    }
                }
                d.read_contained_vals (str.c_str(), linbins);

                // now - if the linbins we loaded were the sbox_angles, then we need to
                // sort linbins, while sorting, at the same time, the signals.

                vector<array<float,12>> flatsurf_boxes;
                array<float, 12> sbox;
                for (unsigned int j = 1; j < linbins.size(); ++j) {
                    // c1 x,y,z
                    sbox[0] = xx-thickness; // x
                    sbox[1] = linbins[j-1]; // y
                    sbox[2] = 0.0;          // z
                    // c2 x,y,z
                    sbox[3] = xx-thickness;
                    sbox[4] = linbins[j];
                    sbox[5] = 0.0;
                    // c3 x,y,z
                    sbox[6] = xx;
                    sbox[7] = linbins[j];
                    sbox[8] = 0.0;
                    // c4 x,y,z
                    sbox[9] = xx;
                    sbox[10] = linbins[j-1];
                    sbox[11] = 0.0;

                    // For angle based view, have to tweak the box angles by +-2pi to avoid boxes that span the whole ribbon.
                    if (co.flattened_type == 2) {
                        if (std::abs(sbox[1] - sbox[4])  > morph::PI_F && std::abs(sbox[7] - sbox[10]) > morph::PI_F) {
                            //std::cout << "PROBLEM BOX\n";
                            if (sbox[4] < 0) {
                                sbox[4] += morph::TWO_PI_F;
                                sbox[7] += morph::TWO_PI_F;
                            } else {
                                sbox[4] -= morph::TWO_PI_F;
                                sbox[7] -= morph::TWO_PI_F;
                            }
                        }
                    }
#if 0 // DEBUG
                    std::cout << "c1: " << sbox[0] << "," << sbox[1] << "," << sbox[2];
                    std::cout << "   c2: " << sbox[3] << "," << sbox[4] << "," << sbox[5];
                    std::cout << "   c3: " << sbox[6] << "," << sbox[7] << "," << sbox[8];
                    std::cout << "   c4: " << sbox[9] << "," << sbox[10] << "," << sbox[11] << std::endl;
#endif
                    flatsurf_boxes.push_back (sbox);
                }
                fquads.insert (fquads.end(), flatsurf_boxes.begin(), flatsurf_boxes.end());
                fmeans.insert (fmeans.end(), frameMeansF.begin(), --frameMeansF.end());
            }
            unsigned int visId = 0;

            // This is the flattened map; showing it alongside the 3D map for now
            offset[0]=-5.5;
            visId = v.addVisualModel (new morph::QuadsVisual<float> (v.shaderprog,
                                                                     &fquads, offset,
                                                                     &fmeans, scale,
                                                                     morph::ColourMapType::Greyscale));
            v.surfaces_2d.push_back (visId);

            // Add a row of points for the centre marker, for debugging
            vector<morph::Vector<float>> centres_;
            vector<float> centres_id;
            for (int i = 1; i < nf; ++i) {
                stringstream ss;
                ss << "/Frame";
                ss.width(3);
                ss.fill('0');
                ss << i;
                frameName = ss.str();
                string str = frameName+"/class/layer_x";
                d.read_val (str.c_str(), xx);
                morph::Vector<float> cp;
                cp[0] = xx;
                cp[1] = 0;
                cp[2] = 0;
                centres_.push_back (cp);
                centres_id.push_back (0.1f*(float)i);
            }
            visId = v.addVisualModel (new morph::ScatterVisual<float> (v.shaderprog,
                                                                       &centres_, offset,
                                                                       &centres_id, 0.03f, scale,
                                                                       morph::ColourMapType::Jet));
            v.angle_centres.push_back (visId);
        }
    } catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << endl;
        rtn = -1;
    }

    return rtn;
}

int main (int argc, char** argv)
{
    int rtn = -1;

    // popt command line argument processing setup
    zeroCmdOptions (&cmdOptions);

    struct poptOption opt[] = {
        POPT_AUTOHELP

        {"scale_perslice", 'p',
         POPT_ARG_NONE, &(cmdOptions.scale_perslice), 0, // 0 is 'val' which is available in callbacks
         "If set, auto-scale the signal for each slice."},

        {"use_autoalign", 'a',
         POPT_ARG_NONE, &(cmdOptions.use_autoalign), 0,
         "If set, prefer the auto-alignment (based on the curve points only) rather than landscape alignment."},

        {"show_ribbons", 'r',
         POPT_ARG_NONE, &(cmdOptions.show_ribbons), 0,
         "If set, display the ribbon-like surface boxes, rather than the smoothed surface."},

        {"show_landmarks", 'l',
         POPT_ARG_INT, &(cmdOptions.show_landmarks), 0,
         "Display landmarks from the indexed (counting from 1) data file."},

        {"show_landmarks_all", 'L',
         POPT_ARG_NONE, &(cmdOptions.show_landmarks), 0,
         "If set, display landmarks from ALL first data files."},

        {"show_flattened", 'm',
         POPT_ARG_INT, &(cmdOptions.show_flattened), 0,
         "Display flattened image from the data file with this index (counting from 1), to compare with 3D."},

        {"flattened_type", 't',
         POPT_ARG_INT, &(cmdOptions.flattened_type), 0,
         "Selects the type of flattened map to show: 0 (default): use aligned 3D image and compute linear "
         "distance along the curve starting from the surface box that is on the zero degree line wrt the "
         "origin. 1: Use centered surface box linear distance. 2: Plot vs. the angle of the surface box."},

        // options following this will cause the popt_option_callback to be executed.
        { "callback", '\0',
          POPT_ARG_CALLBACK|POPT_ARGFLAG_DOC_HIDDEN, (void*)&popt_option_callback, 0,
          NULL, NULL },

        {"datafile", 'f',
         POPT_ARG_STRING, &(cmdOptions.datafile), 0,
         "Add a data file to visualise in 3D. Provide an argument like /path/to/file.h5. "
         "This option can be used multiple times, and you can even leave the -f out; any "
         "'non-option' strings on your command line will be interpreted as data files."},

        POPT_AUTOALIAS
        POPT_TABLEEND
    };
    poptContext con;
    con = poptGetContext (argv[0], argc, (const char**)argv, opt, 0);
    while (poptGetNextOpt(con) != -1) {}
    const char* argg = (char*)0;
    while ((argg = poptGetArg(con)) != (char*)0) {
        // Treat any extra args as files.
        //string df(argg);
        cmdOptions.datafiles.push_back (string(argg));
    }

    try {
        if (cmdOptions.datafiles.empty()) {
            throw std::runtime_error ("Please supply at least one HDF5 file with the -f option.");
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        rtn = -1;
        poptFreeContext(con);
        return rtn;
    }
    // End processing options. Can now access options via cmdOptions

    // Visual scene for all models
    morph::Vector<float> coordArrowLocn = {0,0,0};
    morph::Vector<float> coordArrowLengths = {6,2,2};
    SFVisual v(1024, 768, cmdOptions.datafiles[0], coordArrowLocn, coordArrowLengths, 0.4);
    v.zNear = 0.001;
    v.zFar = 40.0;
    v.setZDefault (-15.4);

    // For each file in cmdOptions.datafiles:
    // 0 red .2 yellow .3 green .4 cyan green .5 cyan .6 blue .7 blue .8 purple .9 red
    vector<float> hues = {0.0f, 0.7f, 0.8f, 0.1f, 0.5f, 0.6f, 0.1f, 0.8f};
    float hue = 0.0f;
    for (unsigned int ii = 0; ii < cmdOptions.datafiles.size(); ++ii) {
        if (ii < hues.size()) {
            hue = hues[ii];
        } else {
            hue = 0.8f;
        }
        rtn += addVisMod (v, cmdOptions.datafiles[ii], cmdOptions, hue);
    }

    // And landmarks
    std::cout << "show_landmarks: " << cmdOptions.show_landmarks << "\n";
    if (cmdOptions.show_landmarks_all) {
        for (unsigned int ii = 0; ii < cmdOptions.datafiles.size(); ++ii) {
            rtn += addLandmarks (v, cmdOptions.datafiles[ii], cmdOptions);
        }
    } else {
        if (cmdOptions.show_landmarks > 0 && static_cast<int>(cmdOptions.datafiles.size()) >= cmdOptions.show_landmarks) {
            rtn += addLandmarks (v, cmdOptions.datafiles[cmdOptions.show_landmarks-1], cmdOptions);
        }
    }

    // Add requested flattened map
    std::cout << "show_flattened: " << cmdOptions.show_flattened << "\n";
    if (cmdOptions.show_flattened > 0 && static_cast<int>(cmdOptions.datafiles.size()) >= cmdOptions.show_flattened) {
        rtn += addFlattened (v, cmdOptions.datafiles[cmdOptions.show_flattened-1], cmdOptions);
    }

    try {
        v.render();
        while (v.readyToFinish == false) {
            glfwWaitEventsTimeout (0.018);
            v.render();
        }
    } catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << endl;
        rtn = -1;
    }

    poptFreeContext(con);
    return rtn;
}
