#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <stdexcept>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>
#include <bitset>
#include <limits>
#include <morph/BezCurvePath.h>
#include <morph/BezCurve.h>
#include <morph/BezCoord.h>
#define BUILD_HDFDATA_WITH_OPENCV 1
#include <morph/HdfData.h>
#include <morph/MathAlgo.h>
#include <morph/NM_Simplex.h>
#include <morph/MathConst.h>
#include <morph/Winder.h>

//! This is the input mode (or you could think of them as tools). So in Bezier mode, you
//! add points for the curve fitting; in freehand mode you draw a loop, and in Landmark
//! mode, you give landmarks for slice alignment.
enum class InputMode
{
    Bezier,    // Cubic Bezier: "curve drawing mode"
    Freehand,  // A freehand drawn loop enclosing a region
    Landmark,  // User provides alignment landmark locations on each slice
    GlobalLandmark, // User provides a number of landmarks, numbered across the entire slice
                    // set, for linear transformations between individual datasets.
    ReverseBezier, // Curve drawing mode but adding/deleting points at the start of the curve
    Circlemark, // Circular landmarks that are large and require 3 points to estimate
                // their centre. Developed to handle needle alignment holes.
    Axismark   // Allows user to define two locations that mark a linear axis through
               // the brain. To help construct nice 2D maps from 3D reconstructions.
};

//! What sort of colour model is in use?
enum class ColourModel
{
    Greyscale,     // Regular greyscale image
    AllenDevMouse, // Coloured images as found on the Allen Developing Mouse Brain Atlas
    Sfview         // Signal data loaded from a sfview-generated 'unwrapped 2D map'
};

enum FrameFlag
{
    ShowBoxes, // Show the yellow boxes?
    ShowUsers, // Show the user points?
    ShowCtrls, // Show the ctrl points of the fits?
    ShowFits,  // Show the fits?
    Mirrored,  // Is the image mirrored L-R?
    Flipped    // Is the image flipped U-D?
};

struct AllenColourParams
{
    std::array<float, 9> colour_rot; // Colour space rotation
    std::array<float, 3> colour_trans; // Colour space pre-translation
    std::array<float, 2> ellip_axes; // red-green ellipse for "elliptical tube of expressing colours
    float luminosity_factor; // The slope of the linear luminosity vs signal fit.
    float luminosity_cutoff; // at what luminosity does the signal cut off to zero?
};

/*!
 * A class to hold a cortical section image, user-supplied cortex edge points and the
 * resulting fit (either polynomial or Bezier curved). Also stores information about
 * freehand drawn loops (whose content mean luminance can be saved out) and
 * user-supplied landmarks, allowing for the alignment of multiple brain slices.
 */
class FrameData
{
private:
    //! The 'previous' frame in the stack of frames (index into)
    int previous = -1;
    std::vector<FrameData>* parentStack;

    //! Number of bins to create for the fit (one less than nFit)
    int nBins;
    //! Number of points to create in the fit
    int nFit;

public:
    //! What input mode is default?
    InputMode ct = InputMode::Bezier;

    // Bezier curve attributes

    //! A Bezier curve path to fit the cortex.
    morph::BezCurvePath<double> bcp;
    //! The user-supplied points from which to make a curve fit. Added to end of PP.
    std::deque<cv::Point> P;
    //! user-supplied points for potential addition to PP at its *start*.
    std::deque<cv::Point> sP;
    //! A deque or deques of points for multi-section Bezier curves
    std::deque<std::deque<cv::Point>> PP;
    //! Index into PP
    int pp_idx = 0;

    // Landmark attributes

    //! The landmark points for this frame.
    std::vector<cv::Point> LM;
    //! The landmark points scaled by pixels_per_mm
    std::vector<cv::Point2d> LM_scaled;
    //! We'll also save out the autoaligned landmarks
    std::vector<cv::Point2d> LM_autoalign_translated;
    std::vector<cv::Point2d> LM_autoaligned;
    //! As part of alignment, have to hold a copy of the aligned landmarks
    std::vector<cv::Point2d> LM_lmaligned;

    //! A point, on the fitted curve, that is defined as being the 'centre' of the
    //! curve. This is used as a starting point from which to unwrap the 3D map into a
    //! 2D map. This is an index into fitted_lmaligned
    int centre_lmaligned = 0;
    //! Corresponding centre for auto-aligned curves - index into fitted_autoaligned
    int centre_autoaligned = 0;

    // Circlemarks. For every 3 circlemarks, we effectively add a landmark. So
    // circlemark mode simply adds landmarks, but by placing 3 points on a circle. The
    // difference comes when deleting, I suppose. Usually don't expect to see more than
    // 3 of these? Do I save the circle marks for each landmark or have the system
    // re-set them each time? If I don't save, then the off-screen landmarks will be
    // missing.

    //! Current circle mark points, user-supplied
    std::vector<cv::Point> CM_points;
    //! Holds a collection of finished circlemark triplets, ordered by the key, which is
    //! the index into LM for the associated landmark, which is the centre of the circle
    //! defined by the triplets.
    std::map<size_t, std::vector<cv::Point>> CM;

    //! Axismark points. For linear axis, just have 2 of these in the whole set of
    //! slices, so may be unset for most FrameData instances. This is a vector a) to
    //! match the other 'marks' and b) in case in the future it might be desirable to
    //! define more than one axis through the brain.
    std::vector<cv::Point> AM;
    //! Axismark points scaled by pixels_per_mm
    std::vector<cv::Point2d> AM_scaled;
    //! We'll also save out the autoaligned axismarks, just like we save out autoaligned landmarks
    std::vector<cv::Point2d> AM_autoalign_translated;
    std::vector<cv::Point2d> AM_autoaligned;
    //! As part of alignment, have to hold a copy of the aligned axismarks
    std::vector<cv::Point2d> AM_lmaligned;

    //! This holds the offset, on this frame, computed from axis marks on other
    //! frames. Computed by DM before frames are written. This offset is used when
    //! generating 2D maps in FrameData::write()
    std::vector<cv::Point2d> AM_origins_lmaligned;
    std::vector<cv::Point2d> AM_origins_autoaligned;

    //! Global landmarks are a bit like axismarks, in that there are some of them, across
    //! the whole slice set. It's possible that there are >1 on a given slice, too.
    std::vector<cv::Point> GLM;
    //! Global landmark points scaled by pixels_per_mm
    std::vector<cv::Point2d> GLM_scaled;
    //! We'll also save out the autoaligned global landmarks, just like we save out autoaligned landmarks
    std::vector<cv::Point2d> GLM_autoalign_translated;
    std::vector<cv::Point2d> GLM_autoaligned;
    //! As part of alignment, have to hold a copy of the aligned global landmarks
    std::vector<cv::Point2d> GLM_lmaligned;

    // May need:
    //std::vector<cv::Point2d> GLM_origins_lmaligned;
    //std::vector<cv::Point2d> GLM_origins_autoaligned;

    //! The means computed for the boxes. This is "mean_signal".
    std::vector<float> box_signal_means;
    //! And corresponding standard deviations
    std::vector<float> box_signal_sds;
    //! Mean pixel value in a box.
    std::vector<unsigned int> box_pixel_means;
    //! And corresponding standard deviations
    std::vector<unsigned int> box_pixel_sds;
    //! The raw pixel values for each box as a vector of unsigned ints for each box.
    std::vector<std::vector<unsigned int>> boxes_pixels;
    //! The coordinates, in pixels, of the pixels in each box.
    std::vector<std::vector<cv::Point>> box_coords_pixels;
    //! The coordinates, in mm, of the pixels in each box in the autoalign coordinate system
    std::vector<std::vector<cv::Point2d>> box_coords_autoalign;
    //! The coordinates, in mm, of the pixels in each box in the lmalign coordinate system
    std::vector<std::vector<cv::Point2d>> box_coords_lmalign;
    //! The signal values for each box as a vector of floats for each box.
    std::vector<std::vector<float>> boxes_signal;
    //! Raw colours of boxes in RGB
    std::vector<std::vector<std::array<float, 3>>> boxes_pixels_bgr;

    //! Bin A sets the number of pixels along the normal to the fit to start the sample box
    int binA = 0;
    //! Bin A sets the number of pixels along the normal to the fit to end the sample box
    int binB = 100;

    //! A set of points created from the fit. Units: pixels, but stored in double precision.
    std::vector<cv::Point2d> fitted;
    //! Take FrameData::fitted and scale by pixels_per_mm. Units: mm.
    std::vector<cv::Point2d> fitted_scaled;

    //! The centroid of the curve, after lm_translation or autoalign_translation have
    //! been applied to move the user points This will be 0,0 in the latter case (in
    //! fact, it's only initially going to be used for LM alignment).
    cv::Point2d curve_centroid;

    // Attributes to do with the auto-transformed slices, where the centroid of each
    // curve is used to arrange the slices on an axis, then each slice's curve is
    // rotated to best-fit the previous curve.

    //! (-) the centroid of fitted. Units: mm.
    cv::Point2d autoalign_translation;
    //! This holds scaled then translated fitted points: fitted_scaled - autoalign_translation
    std::vector<cv::Point2d> fitted_autoalign_translated;
    //! The rotation of this slice, with respect to the previous slice. Used by slice alignment algorithms.
    double autoalign_theta = 0.0;
    //! This holds the fitted_autoalign_translated points after they have been rotated to be in line
    //! with fitted_autoaligned points in the 'previous' frame. 'in line' means the
    //! smallest sum-of-square distances between the two fitted_autoaligned sets. Depends
    //! on nBins being the same in each.
    std::vector<cv::Point2d> fitted_autoaligned;
    //! Set true if centroid-and-rotate alignment was completed before write()
    bool autoalignComputed = false;

    //! The translation applied to fitted_scaled to get to fitted_lmalign_translated
    cv::Point2d lm_translation;
    //! Holds the scaled, fitted points after they have been translated by lm_translation.
    std::vector<cv::Point2d> fitted_lmalign_translated;
    //! The rotation applied to fitted_lmalign_translated to get to fitted_lmaligned
    double lm_theta = 0.0;
    //! Holds the final, translated and rotated points aligned with landmarks.
    std::vector<cv::Point2d> fitted_lmaligned;
    //! If true, and there are >1 landmark per slice, apply the "rotate slices about
    //! landmark 1" alignment procedure anyway. This rotational alignment is applied by
    //! default if there is ONLY 1 landmark per slice.
    bool rotateLandmarkOne = false;
    //! If true, in "rotate about landmark 1 mode" align the other landmarks, instead of the curves.
    bool rotateButAlignLandmarkTwoPlus = false;
    //! Set true if landmark alignment was completed before write()
    bool lmalignComputed = false;

    //! For point in fitted, the tangent at that location. Units: pixels.
    std::vector<cv::Point2d> tangents;
    //! For point in fitted, the normal at that location. Units: pixels.
    std::vector<cv::Point2d> normals;

    //! origins for the lines making the box sides (some distance from the curve)
    std::vector<cv::Point> pointsInner;
    //! endpoints for the lines making the box sides (a greater distance from the curve)
    std::vector<cv::Point> pointsOuter;
    //! The boxes that are drawn and from which to sample the gene expression
    std::vector<std::vector<cv::Point> > boxes;

    //! A vector of user-supplied points for the Freehand drawn loop
    std::vector<cv::Point> FL;
    //! All loops - boundary pixels
    std::vector<std::vector<cv::Point>> FLB;
    //! Extents of the loop FL
    std::array<cv::Point, 2> extents_FL;
    //! vector of vectors containing the points enclosed by the path FL
    std::vector<std::vector<cv::Point>> FLE; // FL_coords_pixels
    //! To hold the  points enclosed by the path FL transformed according to autoalign
    std::vector<std::vector<cv::Point2d>> FLE_autoalign;
    //! To hold the  points enclosed by the path FL transformed according to lmalign
    std::vector<std::vector<cv::Point2d>> FLE_lmalign;

    //! The mean luminance of each freehand loop enclosed region in FLE.
    std::vector<float> FL_signal_means;
    //! The mean pixel value (0-255) in a freehand loop
    std::vector<unsigned int> FL_pixel_means;
    //! The raw pixel values for each region as a vector of unsigned ints
    std::vector<std::vector<unsigned int>> FL_pixels;
    //! The signal values for each region as a vector of floats
    std::vector<std::vector<float>> FL_signal;
    //! Raw values for each region in colour
    std::vector<std::vector<std::array<float, 3>>> FL_pixels_bgr;
    //! Set true if a loop has been drawn, then completed. Reset to false once mouse button is released
    bool loopFinished = false;

    //! A bit set containing flags
    std::bitset<8> flags;
    //! The image data, required when sampling the image in one of the boxes. CV_8UC3
    cv::Mat frame;
    //! If set true, then the read() function will load FrameNNN/class/frame from the H5
    //! file into this->frame and then call this->setupFrames()
    bool loadFrameFromH5 = false;
    //! If true, then save FrameData::frame contents to H5 file on write()
    bool saveFrameToH5 = true;
    //! Maximum and minimum pixel values of frame
    std::pair<int, int> frame_maxmin;
    //! A blurred copy of the image data. CV_32FC3.
    cv::Mat blurred;
    //! Sets the width of the blurring Gaussian's sigma
    double bgBlurScreenProportion = 0.1667;
    //! An offset used when subtracting blurred BG from image.
    float bgBlurSubtractionOffset = 255.0f;
    //! The signal frame, CV_32FC3 format. 1.0f - frame_bgoff. In previous code I've
    //! converted this to CV_8UC3 format because that means that text and lines drawn
    //! onto the photo can be antialiased. To save memory, I've foregone this.
    cv::Mat frame_signal;
    //! Holds the max and min of the signal in frame_signal.
    std::pair<float, float> frame_signal_maxmin;

    //! The frame image filename from which frame was loaded. Stored so it can be
    //! recorded when writing out.
    std::string filename;
    //! The 'x' position of the brain slice in this frame (coordinates in the plane of
    //! the slice are y/z (units: mm).
    float layer_x = 0.0f;
    //! The thickness of this slice (mm)
    float thickness = 0.05;
    //! The scaling factor
    double pixels_per_mm = 100.0;
    //! The index of the frame
    int idx;

    // Colour space parameters.

    //! What kind of colour model is in use?
    ColourModel cmodel = ColourModel::Greyscale;
    //! Colour space rotation to apply to [b g r] colour vectors
    std::array<float, 9> colour_rot;
    //! Colour space pre-translation. [x y z] is [b g r]
    std::array<float, 3> colour_trans;
    //! red-green ellipse for "elliptical tube of expressing colours
    std::array<float, 2> ellip_axes;
    //! The slope of the linear luminosity vs signal fit. Only used for ColourModel::AllenDevMouse
    float luminosity_factor = -0.00392f; // -1/255 aka -1.0f
    //! at what luminosity does the signal cut off to zero? Only used for ColourModel::AllenDevMouse
    float luminosity_cutoff = 255.0f; // aka 1.0f

public:
    FrameData() { throw std::runtime_error ("Default constructor is not allowed"); }

    //! Constructor for creating a FrameData with NO image data. The image data is to be
    //! read from the .h5 file.
    FrameData (const double _bgBlurScreenProportion,
               const float _bgBlurSubtractionOffset,
               const ColourModel _cmodel,
               const AllenColourParams& acparams)
    {
        this->previous = -1;
        this->parentStack = nullptr;
        this->cmodel = _cmodel;
        this->colour_rot = acparams.colour_rot;
        this->colour_trans = acparams.colour_trans;
        this->ellip_axes = acparams.ellip_axes;
        this->luminosity_factor = acparams.luminosity_factor;
        this->luminosity_cutoff = acparams.luminosity_cutoff;
        // An offset so that we don't lose very small amounts of signal above the
        // background when we subtract the blurred version of the image. To become a
        // parameter for the user to modify at application level. Really? Want this here?
        this->bgBlurScreenProportion = _bgBlurScreenProportion;
        if (_bgBlurSubtractionOffset < 0.0f || _bgBlurSubtractionOffset > 255.0f) {
            throw std::runtime_error ("The bg blur subtraction offset should be in range [0,255]");
        }
        this->bgBlurSubtractionOffset = _bgBlurSubtractionOffset;
        // setupFrames will have to wait until later.
        this->loadFrameFromH5 = true;
    }

    //! Constructor initializes default values and calls setupFrames()
    FrameData (const cv::Mat& fr,
               const double _bgBlurScreenProportion,
               const float _bgBlurSubtractionOffset,
               const ColourModel _cmodel,
               const AllenColourParams& acparams)
    {
        std::cout << __FUNCTION__ << " called\n";
        // init previous to null.
        this->previous = -1;
        this->parentStack = nullptr;
        this->cmodel = _cmodel;
        this->colour_rot = acparams.colour_rot;
        this->colour_trans = acparams.colour_trans;
        this->ellip_axes = acparams.ellip_axes;
        this->luminosity_factor = acparams.luminosity_factor;
        this->luminosity_cutoff = acparams.luminosity_cutoff;

        if (this->cmodel == ColourModel::Sfview) {
            // In this case, we should have been passed in a CV_32FC3 image to write
            // into BOTH frame_signal and to convert into frame.
            this->frame_signal = fr.clone(); // fr should be 3 channels to start with.
            cv::Mat temp;
            // Copy frame into temp, multiplying by 255:
            fr.convertTo (temp, CV_8UC3, 255.0);
            // Check number of channels and format of this->frame
            if (temp.channels() == 3) {
                this->frame = temp.clone();
            } else {
                throw std::runtime_error ("Expecting 3 channels");
            }
        } else {
            this->frame = fr.clone();
        }
        this->bgBlurScreenProportion = _bgBlurScreenProportion;
        if (_bgBlurSubtractionOffset < 0.0f || _bgBlurSubtractionOffset > 255.0f) {
            throw std::runtime_error ("The bg blur subtraction offset should be in range [0,255]");
        }
        this->bgBlurSubtractionOffset = _bgBlurSubtractionOffset;
        this->setupFrames();
    }

    //! From the initial frame, as loaded from an image file, or as retreived from the
    //! .h5 file, generate the signal frame
    void setupFrames()
    {
        std::cout << __FUNCTION__ << " called\n";
        this->frame_maxmin = this->showMaxMinU (this->frame, "frame (original)");
        // Scale and convert frame to float format. frameF is a copy of the image data
        // in float format (in range 0 to 1, rather than 0 to 255). CV_32FC3
        cv::Mat frameF;
        this->frame.convertTo (frameF, CV_32FC3, 1/255.0);
        this->showMaxMin (frameF, "frameF (float)");
        // NB: Init these before the next three resize() calls
        this->setBins (100);
        // Init flags
        this->flags.set (ShowFits);
        this->flags.set (ShowUsers);
        this->flags.set (ShowBoxes);

        // Make a blurred copy of the floating point format frame, for estimating lighting background
        this->blurred = cv::Mat::zeros (frameF.rows, frameF.cols, CV_32FC3);
        cv::Size ksz;
        ksz.width = frameF.cols * 2.0 * this->bgBlurScreenProportion;
        ksz.width += (ksz.width%2 == 1) ? 0 : 1; // ensure ksz.width is odd
        ksz.height = frameF.rows/3;
        ksz.height += (ksz.height%2 == 1) ? 0 : 1;
        double sigma = (double)frameF.cols * this->bgBlurScreenProportion;
        cv::GaussianBlur (frameF, this->blurred, ksz, sigma);
        std::cout << "Blur mean = " << cv::mean (this->blurred) << std::endl;

        // Now subtract the blur from the original
        cv::Mat suboffset_minus_blurred;
        cv::subtract (this->bgBlurSubtractionOffset/255.0f, this->blurred,
                      suboffset_minus_blurred, cv::noArray(), CV_32FC3);
        // show max mins (For debugging)
        this->showMaxMin (this->blurred, "this->blurred");
        this->showMaxMin (suboffset_minus_blurred, "suboffset_minus_blurred (const-blurred)");
        // Add suboffset_minus_blurred to frameF to get frame_bgoff.
        cv::Mat frame_bgoff;
        cv::add (frameF, suboffset_minus_blurred, frame_bgoff, cv::noArray(), CV_32FC3);
        this->showMaxMin (frame_bgoff, "frame_bgoff");

        // This is where we distinguish between possible different ColourModels. Apply
        // some conversion to go from the input (frame/frame_bgoff) to the signal:
        // frame_signal
        if (this->cmodel == ColourModel::AllenDevMouse) {
            // Pass in this->frame, which is CV_8UC3 format; the Allen colour conversion
            // depends on RBG values being in range 0-255. Using this->frame, means that
            // we don't subtract an overall background offset, as we do with our own
            // data, but we assume that the Allen process did not have any significant
            // spatial heterogeneity in the illumination.
            this->allenColourConversion (this->frame, this->frame_signal);
        } else if (this->cmodel == ColourModel::Sfview) {
            // In this case, we have data read from a sfview .h5 file (.TF.*.h5) format
            // frame_signal and frame will have been set up in the constructor.
        } else {
            // Assume it's grayscale: Invert frame_bgoff to create the signal frame.
            std::cout << "Treating image as greyscale\n";
            cv::subtract (1.0f, frame_bgoff, this->frame_signal, cv::noArray(), CV_32FC3);
        }

        this->frame_signal_maxmin = this->showMaxMin (this->frame_signal, "frame_signal");
    }

    //! Apply the Allen colour conversion to the frame data to obtain the signal. We
    //! bring in parameters from JSON for transforming the colours and determining if
    //! they're on the "expressing" axis. This will include a translate matrix, a
    //! rotation matrix and ellipse parameters, obtained from the octave script
    //! plotcolour.m. See DM.h for the settings of colour_trans and colour_rot.
    void allenColourConversion (const cv::Mat& inframe, cv::Mat& outframe)
    {
        // inframe should be:
        //          this->frame, CV_8UC3 format
        // outframe is going to be:
        //          this->frame_signal, CV_32FC3 format

        if (inframe.rows != outframe.rows || inframe.cols != outframe.cols) {
            // Make outframe same size as inframe
            cv::Mat z(cv::Size(inframe.cols, inframe.rows), CV_32FC3);
            cv::resize (z, outframe, cv::Size(inframe.cols, inframe.rows));
        }

        float ellip_maj_sq = ellip_axes[0] * ellip_axes[0];
        float ellip_min_sq = ellip_axes[1] * ellip_axes[1];

        // This is "for each pixel in inframe"
        for (int rw = 0; rw < inframe.rows; ++rw) {
            for (int cl = 0; cl < inframe.cols; ++cl) {
                // Perform colour transform here, so that we get a transformed blue value
                cv::Vec3b bgr = inframe.at<cv::Vec3b>(rw,cl);
                //std::cout << "bgr: " << bgr << std::endl;
                // 1. Translate rgb colour. NB: It's this->colour_trans
                float b_t = (float)bgr(0) - colour_trans[0];
                float g_t = (float)bgr(1) - colour_trans[1];
                float r_t = (float)bgr(2) - colour_trans[2];
                //std::cout << "bgr_translated: " << b_t << "," << g_t << "," << r_t << std::endl;

                // 2. Rotate colour. NB: It's this->colour_rot
                float b_r = colour_rot[0]*b_t + colour_rot[1]*g_t + colour_rot[2]*r_t;
                float g_r = colour_rot[3]*b_t + colour_rot[4]*g_t + colour_rot[5]*r_t;
                float r_r = colour_rot[6]*b_t + colour_rot[7]*g_t + colour_rot[8]*r_t;
                //std::cout << "bgr transformed: " << b_r << "," << g_r << "," << r_r << std::endl;

                // if (r_r,g_r) lies inside ellipse (given by ellip_axes) then blue_transf equals b_r
                float blue_transf = this->luminosity_cutoff;
                float erad = ((g_r*g_r)/ellip_maj_sq) + ((r_r*r_r)/ellip_min_sq);
                // If erad <= 1, then we're inside ellipse:
                if (erad <= 1.0f) { blue_transf = b_r; }

                // Now apply signal conversion
                float signal = blue_transf - this->luminosity_cutoff;
                // m * (x - x_0):
                signal *= this->luminosity_factor;
                // Any signal <0 is 0.
                if (signal > 0.0f) {
#ifdef __DEBUG
                    std::cout << "blue_transf = " << blue_transf << std::endl;
                    std::cout << "Set outframe[" << rw << "][" << cl << "] to signal  "
                              << (signal > 0.0f ? signal : 0.0f) << std::endl;
                    std::cout << "cf: bgr: " << bgr << std::endl;
#endif
                    outframe.at<cv::Vec<float, 3> >(rw, cl) = {(signal), (signal), (signal)};
                }
                //outframe.at<float>(rw, cl, 1) = signal > 0.0f ? signal : 0.0f;
                //outframe.at<float>(rw, cl, 2) = signal > 0.0f ? signal : 0.0f;
            }
        }
    }

    //! Show the max and the min of a
    std::pair<float, float> showMaxMin (const cv::Mat& m, const std::string& matlabel = "(unknown)")
    {
        float minm = 100.0f;
        float maxm = -100.0f;
        for (int r = 0; r < m.rows; ++r) {
            for (int c = 0; c < m.cols; ++c) {
                //std::cout << "Frame("<<r<<","<<c<<") = " << m.at< cv::Vec<float, 3> >(r,c) << "\n";
                float val = (float)m.at< cv::Vec<float, 3> >(r,c)[0];
                minm = val < minm ? val : minm;
                maxm = val > maxm ? val : maxm;
            }
        }
        std::cout << "The matrix " << matlabel << "  has min/max: " << minm << "/" << maxm << std::endl;
        return std::make_pair(maxm, minm);
    }
    //! Show the max and the min of a, which should be in U8 format
    std::pair<int, int> showMaxMinU (const cv::Mat& m, const std::string& matlabel = "(unknown)")
    {
        int minm = 256;
        int maxm = -256;
        for (int r = 0; r < m.rows; ++r) {
            for (int c = 0; c < m.cols; ++c) {
                //std::cout << "Frame("<<r<<","<<c<<") = " << m.at< cv::Vec<float, 3> >(r,c) << "\n";
                int val = (int)m.at< cv::Vec<unsigned char, 3> >(r,c)[0];
                minm = val < minm ? val : minm;
                maxm = val > maxm ? val : maxm;
            }
        }
        std::cout << "The matrix " << matlabel << "  has min/max: " << minm << "/" << maxm << std::endl;
        return std::make_pair(maxm, minm);
    }

    //! Getter for the blurred image
    cv::Mat* getBlur() { return &this->blurred; }

    //! This getter returns the _target_ number of bins; the corresponding setter sets
    //! both target and nBins itself.
    int getBins() const { return this->nBins; }

    //! Set the number of bins and update the size of the various containers
    void setBins (int num)
    {
        if (num > 5000) { throw std::runtime_error ("Too many bins..."); }
        this->nBins = num;
        this->nFit = num + 1;
        this->fitted_scaled.resize (this->nFit);
        this->fitted_autoalign_translated.resize (this->nFit);
        this->fitted_autoaligned.resize (this->nFit);
        this->fitted_lmalign_translated.resize (this->nFit);
        this->fitted_lmaligned.resize (this->nFit);
    }

    //! Setter for FrameData::previous, which indexes the previous frame in the stack of
    //! frames (FrameData::parentStack)
    void setPrevious (int prev) { this->previous = prev; }
    //! Setter for the stack of frames (FrameData::parentStack)
    void setParentStack (std::vector<FrameData>* parentSt) { this->parentStack = parentSt; }

    //! Get information about the Bezier fit and the input mode for this frame.
    std::string getFitInfo() const
    {
        std::stringstream ss;
        std::stringstream bb;
        bool first = true;
        if (this->bcp.curves.empty()) {
            bb << "none";
        } else {
            for (auto cv : this->bcp.curves) {
                if (first) {
                    bb << cv.getOrder();
                    first = false;
                } else {
                    bb << "/" << cv.getOrder();
                }
            }
        }
        ss << "Bezier order: " << bb.str() << ", Bins: " << this->nBins << "/A"<< this->binA << ":B" << this->binB;

        if (this->ct == InputMode::Bezier) {
            ss << ". Curve mode (end)";
        } else if (this->ct == InputMode::ReverseBezier) {
            ss << ". Curve mode (start)";
        } else if (this->ct == InputMode::Freehand) {
            // Get any fit info for a freehand loop (e.g. is it contiguous; how many pixels)
            ss << ". Freehand mode";
        } else if (this->ct == InputMode::Landmark) {
            ss << ". Landmark mode";
        } else if (this->ct == InputMode::GlobalLandmark) {
            ss << ". Globalmark mode";
        } else if (this->ct == InputMode::Circlemark) {
            ss << ". Circlemark mode";
        } else if (this->ct == InputMode::Axismark) {
            ss << ". Axismark mode";
        } else {
            ss << ". unknown mode";
        }
        return ss.str();
    }

    //! This is a candidate for MathAlgo; filling squares in between two randomly
    //! chosen squares on a grid. Each square is 1x1 on a grid, with its centre
    //! specified by its cv::Point.  Draw a line between firstSquare and
    //! endSquare. Fill in all pixels which are crossed by the line. Do this with a
    //! recursive algorithm, as it's easy and we're very unlikely to exceed the
    //! recursion limit.
    void fillFL (cv::Point& firstSquare, const cv::Point& endSquare)
    {
        // Finished when the last element of FL is pt.
        if (!this->FL.empty() && this->FL.back() == endSquare) { return; }
        // Start at firstSquare
        if (firstSquare == endSquare) { return; }

        // firstSquare is not the same as endSquare, so fill in between them.  Which of the 8
        // adjoining squares contains the line specified by firstSquare.xy and m?  What's
        // the intersection of the line and the square perimeter of firstSquare? Note, we
        // also need to track currentSquare.

        // Don't need to test if .x==.x AND .y==.y

        // First, push back the firstSquare itself
        this->FL.push_back (firstSquare);

        int xdiff = std::abs(firstSquare.x - endSquare.x);
        int ydiff = std::abs(firstSquare.y - endSquare.y);

        if (firstSquare.y == endSquare.y) {
            // Dirn is east or west
            firstSquare.x += (endSquare.x > firstSquare.x ? 1 : -1);
            return this->fillFL (firstSquare, endSquare);

        } else if (firstSquare.x == endSquare.x) {
            // Dirn is n or s
            firstSquare.y += (endSquare.y > firstSquare.y ? 1 : -1);
            return this->fillFL (firstSquare, endSquare);

        } else if (firstSquare.y < endSquare.y) {
            // SouthWest or SouthEast
            if (firstSquare.x < endSquare.x) {
                // SouthEast
                if (xdiff > ydiff) {
                    // Mark E and move SE
                    firstSquare.x += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y += 1;

                } else if (xdiff < ydiff) {
                    // Mark S and move SE
                    firstSquare.y += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x += 1;

                } else { // xdiff == ydiff
                    // Move SE only
                    firstSquare.x += 1;
                    firstSquare.y += 1;
                }

            } else {
                // SouthWest
                if (xdiff > ydiff) {
                    // Mark W and move SW
                    firstSquare.x -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y += 1;

                } else if (xdiff < ydiff) {
                    // Mark S and move SW
                    firstSquare.y += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x -= 1;

                } else { // xdiff == ydiff
                    // Move SE only
                    firstSquare.x -= 1;
                    firstSquare.y += 1;
                }
            }
        } else {
            // std::cout << "Dirn is NW/NE, computing...\n";
            // North West/East
            if (firstSquare.x < endSquare.x) {
                // NorthEast
                if (xdiff > ydiff) {
                    // Mark E and move NE
                    firstSquare.x += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y -= 1;

                } else if (xdiff < ydiff) {
                    // Mark N and move NE
                    firstSquare.y -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x += 1;

                } else { // xdiff == ydiff
                    // Move NE only
                    firstSquare.x += 1;
                    firstSquare.y -= 1;
                }

            } else {
                // NorthWest
                if (xdiff > ydiff) {
                    // Mark W and move NW
                    firstSquare.x -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y -= 1;

                } else if (xdiff < ydiff) {
                    // Mark N and move NW
                    firstSquare.y -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x -= 1;

                } else { // xdiff == ydiff
                    // Move NW only
                    firstSquare.x -= 1;
                    firstSquare.y -= 1;
                }
            }
        }

        // Recurse
        return this->fillFL (firstSquare, endSquare);
    }

    //! Find rectangular region enclosing loop, return as two coordinates of top left
    //! (minimum x and y, with y going down) and bottom right (max x and y) corners
    std::array<cv::Point, 2> getExtents (const std::vector<cv::Point>& loop)
    {
        std::array<cv::Point, 2> extents;
        extents[0] =  {10000000, 10000000};   // MIN values for x and y
        extents[1] = {-10000000,-10000000}; // MAX values for x and y
        for (auto p : loop) {
            if (p.x < extents[0].x) { extents[0].x = p.x; }
            if (p.y < extents[0].y) { extents[0].y = p.y; }
            if (p.x > extents[1].x) { extents[1].x = p.x; }
            if (p.y > extents[1].y) { extents[1].y = p.y; }
        }
        return extents;
    }

    //! Find all pixels enclosed by the pixels in this->FL which define a loop
    std::vector<cv::Point> getEnclosedByFL()
    {
        // FIXME: Prefer not to have to uniquify here:
        auto last = std::unique(this->FL.begin(), this->FL.end());
        this->FL.erase (last, this->FL.end());

        // First, find extents of the loop
        this->extents_FL = this->getExtents (this->FL);

        // Create a winder object to compute winding numbers
        morph::Winder w (this->FL);

        // It's perhaps inefficient to compute the winding number of EVERY pixel here,
        // but I'll leave it for now (computers are fast).
        std::vector<cv::Point> rtn;
        std::vector<cv::Point> bdry;
        for (int x = this->extents_FL[0].x; x <= this->extents_FL[1].x; ++x) {
            for (int y = this->extents_FL[0].y; y <= this->extents_FL[1].y; ++y) {
                cv::Point px (x, y);
                auto inloop = std::find (this->FL.begin(), this->FL.end(), px);
                if (inloop == this->FL.end()) {
                    // Compute winding number
                    int winding_number = w.wind (px);
                    //std::cout << "Winding number of pixel " << px << " = " << winding_number << std::endl;
                    if (winding_number != 0) {
                        rtn.push_back (px);
                    }
                } else {
                    // else current pixel is a member of the loop itself.
                    int winding_number = w.wind (px);
                    //std::cout << "Winding number of BOUNDARY pixel " << px << " = " << winding_number << std::endl;
                    if (winding_number != 0) {
                        bdry.push_back (px);
                    }
                }
            }
        }

        // push the boundary onto FLB and return the enclosed pixels
        this->FLB.push_back (bdry);
        return rtn;
    }

    //! Add a pixel that was under the mouse pointer to the freehand points FL. Also
    //! add the pixels between the last pixel in FL and endSquare.
    void addToFL (cv::Point& pt)
    {
        // If FL empty then task is simple:
        if (this->FL.empty()) {
             this->FL.push_back (pt);
             return;
        }

        auto existing = std::find (this->FL.begin(), this->FL.end(), pt);
        if (existing != this->FL.end() && existing != this->FL.begin()) {
            // If pt is in FL already, then we closed the loop
            // std::cout << "An existing point has been found, so close the loop.\n";
            // We want to cut off any extraneous pixels in FL.
            // Delete from this->FL from start to one before existing.
            this->FL.erase (this->FL.begin(), existing);
            // Fill the gap between last point on curve and the existing curve point.
            cv::Point fp = this->FL.back();
            this->fillFL (fp, pt);
            // Get enclosed and add to FLE:
            std::vector<cv::Point> inside = this->getEnclosedByFL();
            this->FLE.push_back (inside);
            this->FL.clear();
            this->loopFinished = true;
            this->computeFreehandMeans();
            std::cout << "Loop complete. Mean: " << FL_signal_means.back() << ", Number of pixels: " << FLE.back().size() << "\n";

        } else { // The new point pt is NOT in FL already

            std::pair<float, float> p1, p2;
            float snap_threshold = 3.0;

            // Otherwise, check if we can close the loop...
            p1.first = (float)this->FL.begin()->x;
            p1.second = (float)this->FL.begin()->y;
            p2.first = (float)pt.x;
            p2.second = (float)pt.y;

            if (this->FL.size() > 20 // Avoid joining a just-started loop
                && morph::MathAlgo::distance<float>(p1, p2) < snap_threshold) {
                // The point pt close to the start of the freehand-drawn loop, so "autosnap" to it.
                pt = *(this->FL.begin());
                //std::cout << "Joining the loop!\n";
                cv::Point fp = this->FL.back();
                this->fillFL (fp, pt);
                // Get the enclosed points and add:
                std::vector<cv::Point> inside = this->getEnclosedByFL();
                this->FLE.push_back (inside);
                this->FL.clear();
                this->loopFinished = true;
                // Now loop is finished, compute the means, so they can be displayed in UI
                this->computeFreehandMeans();
                std::cout << "Loop complete. Mean: " << FL_signal_means.back() << ", Number of pixels: " << FLE.back().size() << "\n";

            } else {
                // Can't close the loop; just add to it
                cv::Point fp = this->FL.back();
                this->fillFL (fp, pt);
            }
        }
    }

    //! Remove the last freehand drawn region or the last point, depending on mode
    void removeLastThing()
    {
        if (this->ct == InputMode::Freehand) {
            this->removeLastRegion();
        } else if (this->ct == InputMode::Landmark) {
            this->removeLastLandmark();
        } else if (this->ct == InputMode::GlobalLandmark) {
            this->removeLastGlobalLandmark();
        } else if (this->ct == InputMode::ReverseBezier) {
            this->removeFirstPoint();
        } else if (this->ct == InputMode::Circlemark) {
            // Also removes the associated landmark
            this->removeLastCirclepoint();
        } else if (this->ct == InputMode::Axismark) {
            this->removeLastAxismark();
        } else {
            this->removeLastPoint();
            this->updateFit();
            this->refreshBoxes (-this->binA, this->binB);
        }
    }

    //! Remove the last freehand drawn region
    void removeLastRegion()
    {
        // If there's a half-finished boundary, get rid of that first:
        if (!this->FL.empty()) {
            this->FL.clear();
        } else {
            // Otherwise, remove the last completed freehand drawn region (enclosed pixels and boundary pixels)
            if (!this->FLE.empty()) { this->FLE.pop_back(); }
            if (!this->FLB.empty()) { this->FLB.pop_back(); }
        }
    }

    //! Clear the information for the Bezier fit
    void clearFitBezier()
    {
        this->bcp.reset();
        this->fitted.clear();
        this->tangents.clear();
        this->normals.clear();
    }

    //! Remove all user points - go back to a blank slate
    void removeAllPoints()
    {
        this->clearFitBezier();
        this->clearBoxes();
        if (!this->P.empty()) { this->P.clear(); }
        if (!this->PP.empty()) { this->PP.clear(); }
        this->pp_idx = 0;
    }

    //! Remove the first user point
    void removeFirstPoint()
    {
        if (this->PP.empty() && this->sP.size() == 1) {
            // Normal behaviour, just remove point from sP
            this->sP.pop_back(); // don't need to pop_front.

        } else if (!this->PP.empty() && this->sP.size() == 1) {
            // Remove point from sP and...
            this->sP.pop_back();
            // Because it is the same locn, the first point from PP.front(), too
            this->removeFirstPoint();

        } else if (!this->sP.empty()) {
            this->sP.pop_front();

        } else {
            // sP is empty, go to first curve in PP and remove a point from that
            if (this->ct == InputMode::ReverseBezier && this->pp_idx>0) {
                if (!this->PP.empty()) {
                    this->sP = this->PP[0];
                    this->PP.pop_front();
                    this->pp_idx--;
                    this->sP.pop_front();
                } else {
                    // Catch pathological case where PP is empty, but pp_idx != 0
                    this->pp_idx = 0;
                }
            }
        }
    }

    //! Remove the last user point
    void removeLastPoint()
    {
        if (this->PP.empty() && this->P.size() == 1) {
            // Normal behaviour, just remove point from P
            this->P.pop_back();

        } else if (!this->PP.empty() && this->P.size() == 1) {
            // Remove point from P and...
            this->P.pop_back();
            // Because it is the same locn, the last point from PP.back(), too
            this->removeLastPoint();

        } else if (!this->P.empty()) {
            this->P.pop_back();

        } else {
            // P is empty, go to previous curve and remove a point from that
            if (this->ct == InputMode::Bezier && this->pp_idx>0) {
                if (!this->PP.empty()) {
                    this->P = this->PP[--this->pp_idx];
                    this->PP.pop_back();
                    this->P.pop_back();
                } else {
                    // Catch pathological case where PP is empty, but pp_idx != 0
                    this->pp_idx = 0;
                }
            }
        }
    }

    //! addAxismark. Note this enforces 1 axismark in AM for now.
    void addAxismark (const cv::Point& pt)
    {
        if (this->AM.empty()) { this->AM.push_back (pt); }
    }

    //! Remove the last axis mark coordinate
    void removeLastAxismark()
    {
        if (!this->AM.empty()) {
            this->AM.pop_back();
            size_t AM_sz = this->AM.size();
            while (this->AM_scaled.size() > AM_sz) { this->AM_scaled.pop_back(); }
            while (this->AM_lmaligned.size() > AM_sz) { this->AM_lmaligned.pop_back(); }
            while (this->AM_autoalign_translated.size() > AM_sz) { this->AM_autoalign_translated.pop_back(); }
            while (this->AM_autoaligned.size() > AM_sz) { this->AM_autoaligned.pop_back(); }
            while (this->AM_origins_lmaligned.size() > AM_sz) { this->AM_origins_lmaligned.pop_back(); }
            while (this->AM_origins_autoaligned.size() > AM_sz) { this->AM_origins_autoaligned.pop_back(); }
        }
    }

    //! Add a global landmark
    void addGlobalLandmark (const cv::Point& pt) { this->GLM.push_back (pt); }

    //! Remove the last global landmark
    void removeLastGlobalLandmark()
    {
        if (!this->GLM.empty()) {
            this->GLM.pop_back();
            size_t GLM_sz = this->GLM.size();
            while (this->GLM_scaled.size() > GLM_sz) { this->GLM_scaled.pop_back(); }
            while (this->GLM_autoalign_translated.size() > GLM_sz) { this->GLM_autoalign_translated.pop_back(); }
            while (this->GLM_autoaligned.size() > GLM_sz) { this->GLM_autoaligned.pop_back(); }
            while (this->GLM_lmaligned.size() > GLM_sz) { this->GLM_lmaligned.pop_back(); }
        }
    }

    //! Remove the last landmark coordinate
    void removeLastLandmark()
    {
        if (!this->LM.empty()) {
            // If LM.back() is in this->CM, then delete those, too
            size_t lm_last = this->LM.size()-1;
            std::map<size_t, std::vector<cv::Point>>::iterator cmi = this->CM.find(lm_last);
            if (cmi != this->CM.end()) {
                this->CM.erase (cmi);
            }
            this->LM.pop_back();
        }
    }

    //! Remove the last circlemark point. If there are no current circlepoints in
    //! CM_points, then a) move the largest key entry in this->CM over to CM_points, b)
    //! Remove the landmark associated with it and c) remove the 3rd circle point from
    //! CM_points.
    void removeLastCirclepoint()
    {
        if (this->CM_points.empty()) {
            if (!this->CM.empty()) {
                std::map<size_t, std::vector<cv::Point>>::reverse_iterator ei = this->CM.rbegin();
                std::cout << "End entry in CM has key " << ei->first << std::endl;
                std::vector<cv::Point>::iterator li = this->LM.begin();
                li += ei->first;
                if (li != this->LM.end()) {
                    this->LM.erase (li);
                }
                this->CM_points = ei->second;
                this->CM.erase (this->CM.find(ei->first));
                this->CM_points.pop_back();
            }
        } else {
            this->CM_points.pop_back();
        }
    }

private:
    //! Compute determinant for 3x3 matrix @cm arranged in col major format
    double determinant3x3 (std::array<double, 9> cm) const
    {
        double det = (cm[0]*cm[4]*cm[8])
        + (cm[3]*cm[7]*cm[2])
        + (cm[6]*cm[1]*cm[5])
        - (cm[6]*cm[4]*cm[2])
        - (cm[0]*cm[7]*cm[5])
        - (cm[3]*cm[1]*cm[8]);
        return det;
    }

    //! From 3 points, compute the centre of the circle which passes through all 3.
    //! Matrix based recipe from https://mathworld.wolfram.com/Circumcircle.html
    cv::Point circumcentre (const std::vector<cv::Point>& threepoints)
    {
        // zero is the centre
        cv::Point2d zero;

        cv::Point2d one = cv::Point2d(threepoints[0]);
        cv::Point2d two = cv::Point2d(threepoints[1]);
        cv::Point2d three = cv::Point2d(threepoints[2]);

        std::array<double, 9> Dx, Dy, _a;

        Dx[0] = one.x*one.x + one.y*one.y;         Dx[3] = one.y;   Dx[6] = 1.0;
        Dx[1] = two.x*two.x + two.y*two.y;         Dx[4] = two.y;   Dx[7] = 1.0;
        Dx[2] = three.x*three.x + three.y*three.y; Dx[5] = three.y; Dx[8] = 1.0;

        Dy[0] = one.x*one.x + one.y*one.y;         Dy[3] = one.x;   Dy[6] = 1.0;
        Dy[1] = two.x*two.x + two.y*two.y;         Dy[4] = two.x;   Dy[7] = 1.0;
        Dy[2] = three.x*three.x + three.y*three.y; Dy[5] = three.x; Dy[8] = 1.0;

        _a[0] = one.x; _a[3] = one.y; _a[6] = 1.0;
        _a[1] = two.x; _a[4] = two.y; _a[7] = 1.0;
        _a[2] = three.x; _a[5] = three.y; _a[8] = 1.0;

        double minusbx = this->determinant3x3 (Dx);
        double by = this->determinant3x3 (Dy);
        double twoa = 2.0 * this->determinant3x3 (_a);

        zero.x = minusbx / twoa;
        zero.y = -by / twoa;

        return cv::Point (zero);
    }

public:
    void addCirclepoint (cv::Point& pt)
    {
        this->CM_points.push_back (pt);

        if (this->CM_points.size() == 3) {
            // We just added the last point to CM_points
            cv::Point cent = this->circumcentre (this->CM_points);
            this->LM.push_back (cent);
            this->CM[this->LM.size()-1] = this->CM_points;
            this->CM_points.clear();
        }
    }

    //! In Bezier mode, store the current set of user points (P) into PP and clear P.
    void nextCurve()
    {
        // Or if "cursor was closest to start of curve?"
        if (this->ct == InputMode::ReverseBezier) {
            if (this->sP.size() < 3) { return; }
            this->PP.push_front (this->sP);
            this->sP.clear();
            this->sP.push_front (this->PP.front().front());
            this->pp_idx++;

        } else {
            // Don't add unless at least 3 points to fit:
            if (this->P.size() < 3) { return; }
            this->PP.push_back (this->P);
            this->P.clear();
            this->P.push_back (this->PP.back().back());
            this->pp_idx++;
        }
    }

    //! Read landmark and circlemark points from file (or rather, morph::HdfData object)
    void importLandmarks (morph::HdfData& df)
    {
        std::string frameName = this->getFrameName();

        // Landmark points
        std::string dname = frameName + "/class/LM";
        this->LM.clear();
        try {
            df.read_contained_vals (dname.c_str(), this->LM);
        } catch (...) {
            // Do nothing on exception. Move on to next.
            std::cout << "No landmarks to read for this frame" << std::endl;
        }

        // Circlemarks
        for (size_t i = 0; i < this->LM.size(); ++i) {
            dname = frameName + "/class/CM/lm" + std::to_string(i);
            try {
                std::vector<cv::Point> vpts;
                df.read_contained_vals (dname.c_str(), vpts);
                this->CM[i] = vpts;
            } catch (...) {
                // Do nothing on exception. Move on to next.
                std::cout << "Could not read circlemarks" << std::endl;
            }
        }
        try {
            dname = frameName + "/class/CM_points";
            this->CM_points.clear();
            df.read_val (dname.c_str(), this->CM_points);
        } catch (...) {
            // Do nothing on exception. Move on to next.
            std::cout << "No CM_points to read" << std::endl;
        }

        // Axismark points
        dname = frameName + "/class/AM";
        this->AM.clear();
        try {
            df.read_contained_vals (dname.c_str(), this->AM);
        } catch (...) {
            // Do nothing on exception. Move on to next.
            std::cout << "No axismarks to read for this frame" << std::endl;
        }

        // Global landmark points
        dname = frameName + "/class/GLM";
        this->GLM.clear();
        try {
            df.read_contained_vals (dname.c_str(), this->GLM);
        } catch (...) {
            // Do nothing on exception. Move on to next.
            std::cout << "No global landmarks to read for this frame" << std::endl;
        }
    }

    //! Import the user-supplied coordinates used to create fitted Bezier curves as well
    //! as the sample box info (size and number)
    void importCurves (morph::HdfData& df)
    {
        std::cout << __FUNCTION__ << " called\n";
        std::string frameName = this->getFrameName();
        std::string dname = frameName + "/class/P";
        this->P.clear();
        df.read_contained_vals (dname.c_str(), this->P);
        dname = frameName + "/class/sP";
        this->sP.clear();
        df.read_contained_vals (dname.c_str(), this->sP);

        dname = frameName + "/class/PP_n";
        unsigned int pp_size = 0;
        df.read_val (dname.c_str(), pp_size);
        std::cout << "pp_size = " << pp_size << std::endl;

        this->PP.resize(pp_size);
        for (size_t i = 0; i<pp_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/PP";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.read_contained_vals (ss.str().c_str(), this->PP[i]);
        }

        dname = frameName + "/class/pp_idx";
        df.read_val (dname.c_str(), this->pp_idx);

        dname = frameName + "/class/nBins";
        int _nBins = 0;
        try {
            df.read_val (dname.c_str(), _nBins);
        } catch (...) {
            // In case nBins stored as 'nBinsTarg'
            dname = frameName + "/class/nBinsTarg";
            df.read_val (dname.c_str(), _nBins);
        }
        this->setBins (_nBins);
        dname = frameName + "/class/binA";
        df.read_val (dname.c_str(), this->binA);
        dname = frameName + "/class/binB";
        df.read_val (dname.c_str(), this->binB);
    }

    void importFreehand (morph::HdfData& df)
    {
        std::string frameName = this->getFrameName();
        // Freehand-drawn regions
        std::string dname = frameName + "/class/FL";
        this->FL.clear();
        df.read_contained_vals (dname.c_str(), this->FL);
        dname = frameName + "/class/FLE_n";
        unsigned int fle_size = 0;
        df.read_val (dname.c_str(), fle_size);
        this->FLE.resize(fle_size);
        this->FLB.resize(fle_size);
        for (size_t i = 0; i<fle_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/FLE";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.read_contained_vals (ss.str().c_str(), this->FLE[i]);

            std::stringstream ss1;
            ss1 << frameName + "/class/FLB";
            ss1.width(3);
            ss1.fill('0');
            ss1 << i;
            df.read_contained_vals (ss1.str().c_str(), this->FLB[i]);
        }
    }

    //! Read important data from file
    void read (morph::HdfData& df, bool oldformat=false)
    {
        // Note this file assumes idx has been set for the frame.
        std::string frameName = this->getFrameName();

        if (oldformat == true) {
            throw std::runtime_error ("Note: there is currently no old format conversion code.");
            // NB: This code is left as a place holder in case we need to read in an old
            // format and save in a new format.
        }

        this->importCurves (df);
        this->importFreehand (df);
        this->importLandmarks (df);

        std::string dname = frameName + "/class/flags";
        df.read_val (dname.c_str(), this->flags);
        dname = frameName + "/class/filename";
        df.read_string (dname.c_str(), this->filename);

        if (this->loadFrameFromH5 == true) {
            std::cout << "Attempting to read frame from H5 file...\n";
            try {
                dname = frameName + "/class/frame";
                df.read_contained_vals (dname.c_str(), this->frame);
            } catch (...) {
                std::stringstream ee;
                ee << "Failed to load frame image data from HDF5 file at path " << dname;
                throw std::runtime_error (ee.str());
            }
            this->setupFrames();
        }
        // NB: Don't read back from H5 file those variables which are specified in json,
        // thus changes in json will override the h5 file.
    }

    //! If false, then don't save the per-pixel data
    bool savePerPixelData = false;
    //! If true, save data relating to the autoaligned slices
    bool saveAutoAlignData = true;
    //! If true, save data relating to the landmark-aligned slices
    bool saveLMAlignData = true;

    //! Export the user-supplied points for curve drawing
    void exportCurves (morph::HdfData& df) const
    {
        std::string frameName = this->getFrameName();
        std::string dname = frameName + "/class/P";
        df.add_contained_vals (dname.c_str(), this->P);
        dname = frameName + "/class/sP";
        df.add_contained_vals (dname.c_str(), this->sP);

        dname = frameName + "/class/PP_n";
        unsigned int pp_size = this->PP.size();
        df.add_val (dname.c_str(), pp_size);
        for (size_t i = 0; i<pp_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/PP";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.add_contained_vals (ss.str().c_str(), this->PP[i]);
        }
        dname = frameName + "/class/pp_idx";
        df.add_val (dname.c_str(), this->pp_idx);
        dname = frameName + "/class/nBins";
        df.add_val (dname.c_str(), this->nBins);
        dname = frameName + "/class/binA";
        df.add_val (dname.c_str(), this->binA);
        dname = frameName + "/class/binB";
        df.add_val (dname.c_str(), this->binB);
    }

    //! Export freehand drawn regions
    void exportFreehand (morph::HdfData& df) const
    {
        std::string frameName = this->getFrameName();
        std::string dname = frameName + "/class/FL";
        df.add_contained_vals (dname.c_str(), this->FL);
        dname = frameName + "/class/FLE_n";
        unsigned int fle_size = this->FLE.size();
        df.add_val (dname.c_str(), fle_size);
        for (size_t i = 0; i<fle_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/FLE";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.add_contained_vals (ss.str().c_str(), this->FLE[i]);

            std::stringstream ss1;
            ss1 << frameName + "/class/FLB";
            ss1.width(3);
            ss1.fill('0');
            ss1 << i;
            df.add_contained_vals (ss1.str().c_str(), this->FLB[i]);
        }
    }

    //! Export landmark, axismark and circlemark points for this frame
    void exportLandmarks (morph::HdfData& df) const
    {
        std::string frameName = this->getFrameName();

        // The landmark points
        std::string dname = frameName + "/class/LM";
        df.add_contained_vals (dname.c_str(), this->LM);
        dname = frameName + "/class/LM_scaled";
        df.add_contained_vals (dname.c_str(), this->LM_scaled);

        // Circlemark points
        if (!this->CM_points.empty()) {
            dname = frameName + "/class/CM_points";
            df.add_contained_vals (dname.c_str(), this->CM_points);
        }
        // CM is map<size_t, vector<cv::Point>>. Unpack here and save.
        std::map<size_t, std::vector<cv::Point>>::const_iterator cmi = this->CM.cbegin();
        while (cmi != this->CM.cend()) {
            dname = frameName + "/class/CM/lm" + std::to_string(cmi->first);
            df.add_contained_vals (dname.c_str(), cmi->second);
            ++cmi;
        }

        // The axismark points
        dname = frameName + "/class/AM";
        df.add_contained_vals (dname.c_str(), this->AM);
        dname = frameName + "/class/AM_scaled";
        df.add_contained_vals (dname.c_str(), this->AM_scaled);

        dname = frameName + "/lmalign/alignmark_origins";
        df.add_contained_vals (dname.c_str(), this->AM_origins_lmaligned);
        dname = frameName + "/autoalign/alignmark_origins";
        df.add_contained_vals (dname.c_str(), this->AM_origins_autoaligned);

        // Global landmark points
        dname = frameName + "/class/GLM";
        df.add_contained_vals (dname.c_str(), this->GLM);
        dname = frameName + "/class/GLM_scaled";
        df.add_contained_vals (dname.c_str(), this->GLM_scaled);
    }

    //! A subroutine of FrameData::write.
    void saveAutoAlignAngles (morph::HdfData& df, const double theta_middle,
                              const std::vector<std::array<float,3>>& surface_box_centroids_autoaligned)
    {
        std::string frameName = this->getFrameName();

        std::vector<double> autoalign_angles (this->fitted_autoaligned.size()-1, 0.0);
        for (size_t i = 1; i < this->fitted_autoaligned.size(); ++i) {
            if (this->AM_origins_autoaligned.empty()) {
                autoalign_angles[i-1] = std::atan2 (0.5 * (this->fitted_autoaligned[i].x + this->fitted_autoaligned[i-1].x),
                                                    0.5 * (this->fitted_autoaligned[i].y + this->fitted_autoaligned[i-1].y)) - theta_middle;
            } else {
                // We have a user-supplied brain axis, so compute angles about that axis.
                cv::Point2d angle_origin = this->AM_origins_autoaligned[0];
                cv::Point2d a = this->fitted_autoaligned[i] - angle_origin;
                cv::Point2d b = this->fitted_autoaligned[i-1] - angle_origin;
                autoalign_angles[i-1] = std::atan2 (0.5 * (a.x + b.x), 0.5 * (a.y + b.y)) - theta_middle;

                // Correct angles to lie in range (-pi, pi]
                if (autoalign_angles[i-1] > morph::PI_D) {
                    autoalign_angles[i-1] -= morph::TWO_PI_D;
                } else if (autoalign_angles[i-1] <= -morph::PI_D) {
                    autoalign_angles[i-1] += morph::TWO_PI_D;
                }
            }
        }
        std::string dname = frameName + "/autoalign/flattened/sbox_angles";
        df.add_contained_vals (dname.c_str(), autoalign_angles);

        // Linear distance, but centered using the auto aligned slices
        float total_linear_distance = 0.0f;
        float middle_distance = 0.0f;
        std::vector<float> linear_distances (this->nBins, 0.0f);
        for (int i=1; i<this->nBins; ++i) {
            // Compute distance from Previous to current
            float d = morph::MathAlgo::distance<float> (surface_box_centroids_autoaligned[i-1],
                                                        surface_box_centroids_autoaligned[i]);
            total_linear_distance += d;
            linear_distances[i] = total_linear_distance;
            if (i == this->centre_autoaligned) { middle_distance = total_linear_distance; }
        }
        // Now offset the linear distances so that the middle is the 0 degree location
        for (int i=0; i<this->nBins; ++i) { linear_distances[i] -= middle_distance; }

        dname = frameName + "/autoalign/flattened/sbox_linear_distance";
        df.add_contained_vals (dname.c_str(), linear_distances);
    }

    //! A subroutine of FrameData::write.
    void saveLmAlignAngles (morph::HdfData& df, const double theta_middle,
                            const std::vector<std::array<float,3>>& surface_box_centroids_lmaligned)
    {
        std::string frameName = this->getFrameName();
        // Compute the angles about the coordinate zero axis
        std::vector<double> lmalign_angles (this->fitted_lmaligned.size()-1, 0.0);
        for (size_t i = 1; i < this->fitted_lmaligned.size(); ++i) {
            if (this->AM_origins_lmaligned.empty()) {
                lmalign_angles[i-1] = std::atan2 (0.5 * (this->fitted_lmaligned[i].x + this->fitted_lmaligned[i-1].x),
                                                  0.5 * (this->fitted_lmaligned[i].y + this->fitted_lmaligned[i-1].y)) - theta_middle;
            } else {
                // We have a user-supplied brain axis, so compute angles about that axis.
                cv::Point2d angle_origin = this->AM_origins_lmaligned[0];
                cv::Point2d a = this->fitted_lmaligned[i] - angle_origin;
                cv::Point2d b = this->fitted_lmaligned[i-1] - angle_origin;
                lmalign_angles[i-1] = std::atan2 (0.5 * (a.x + b.x), 0.5 * (a.y + b.y)) - theta_middle;

                // Correct angles to lie in range (-pi, pi]
                if (lmalign_angles[i-1] > morph::PI_D) {
                    lmalign_angles[i-1] -= morph::TWO_PI_D;
                } else if (lmalign_angles[i-1] <= -morph::PI_D) {
                    lmalign_angles[i-1] += morph::TWO_PI_D;
                }
                //std::cout << "Set lmalign_angles[" << (i-1) << "] to " << lmalign_angles[i-1] << std::endl;
            }
        }
        std::string dname = frameName + "/lmalign/flattened/sbox_angles";
        df.add_contained_vals (dname.c_str(), lmalign_angles);

        // Linear distance, but centered using the landmark aligned slices
        float total_linear_distance = 0.0f;
        float middle_distance = 0.0f;
        std::vector<float> linear_distances (this->nBins, 0.0f);
        for (int i=1; i<this->nBins; ++i) {
            // Compute distance from Previous to current
            float d = morph::MathAlgo::distance<float> (surface_box_centroids_lmaligned[i-1],
                                                        surface_box_centroids_lmaligned[i]);
            total_linear_distance += d;
            linear_distances[i] = total_linear_distance;
            if (i == this->centre_lmaligned) { middle_distance = total_linear_distance; }
        }
        // Now offset the linear distances so that the middle is the 0 degree location
        for (int i=0; i<this->nBins; ++i) { linear_distances[i] -= middle_distance; }

        dname = frameName + "/lmalign/flattened/sbox_linear_distance";
        df.add_contained_vals (dname.c_str(), linear_distances);
    }

    //! Write the data out to an HdfData file \a df. theta_middle is an angular offset
    //! for the 'zero angle' about either the origin x axis, or the axis specified by
    //! axismarks.
    void write (morph::HdfData& df, const double theta_middle)
    {
        // Update box means. not const
        this->computeBoxMeans();
        // And any freehand regions. not const
        this->computeFreehandMeans();

        // Export curve points, freehand regions, landmarks and circlemarks
        this->exportCurves (df);
        this->exportFreehand (df);
        this->exportLandmarks (df);

        std::string frameName = this->getFrameName();
        std::string dname = frameName + "/class/flags";
        df.add_val (dname.c_str(), this->flags);
        dname = frameName + "/class/filename";
        df.add_string (dname.c_str(), this->filename);
        dname = frameName + "/class/layer_x";
        df.add_val (dname.c_str(), this->layer_x);
        dname = frameName + "/class/thickness";
        df.add_val (dname.c_str(), this->thickness);
        dname = frameName + "/class/pixels_per_mm";
        df.add_val (dname.c_str(), this->pixels_per_mm);
        dname = frameName + "/class/idx";
        df.add_val (dname.c_str(), this->idx);
        // Write the BG blurring parameters
        dname = frameName + "/class/bg_blur_screen_proportion";
        df.add_val (dname.c_str(), this->bgBlurScreenProportion);
        dname = frameName + "/class/bg_blur_subtraction_offset";
        df.add_val (dname.c_str(), this->bgBlurSubtractionOffset);

        if (this->cmodel == ColourModel::AllenDevMouse) {
            dname = frameName + "/class/luminosity_factor";
            df.add_val (dname.c_str(), this->luminosity_factor);
            dname = frameName + "/class/luminosity_cutoff";
            df.add_val (dname.c_str(), this->luminosity_cutoff);
        }

        // The frame data itself
        if (this->saveFrameToH5 == true) {
            dname = frameName + "/class/frame";
            df.add_contained_vals (dname.c_str(), this->frame);
        }

        /*
         * The rest of the methods write out data that WON'T be read by the
         * FrameData::read method (these would all be re-computed before being
         * re-written in a later run of the program).
         */
        if (this->savePerPixelData == true) {
            for (size_t bi = 0; bi < this->boxes_pixels.size(); ++bi) {
                dname = frameName + "/signal/bits8/boxes/box" + std::to_string(bi);
                df.add_contained_vals (dname.c_str(), this->boxes_pixels[bi]);
                dname = frameName + "/signal/postproc/boxes/box" + std::to_string(bi);
                df.add_contained_vals (dname.c_str(), this->boxes_signal[bi]);
            }
        }
        dname = frameName + "/nboxes";
        df.add_val (dname.c_str(), static_cast<unsigned int>(this->boxes_pixels.size()));

        // Save the mean and standard deviation of the signal values AND the pixel
        // values for each sampling box in the frame.
        dname = frameName + "/signal/postproc/boxes/means";
        df.add_contained_vals (dname.c_str(), this->box_signal_means);
        dname = frameName + "/signal/bits8/boxes/means";
        df.add_contained_vals (dname.c_str(), this->box_pixel_means);
        dname = frameName + "/signal/postproc/boxes/sds";
        df.add_contained_vals (dname.c_str(), this->box_signal_sds);
        dname = frameName + "/signal/bits8/boxes/sds";
        df.add_contained_vals (dname.c_str(), this->box_pixel_sds);

        // Autoscale box_signal_means and save a copy
        dname = frameName + "/signal/postproc/boxes/means_autoscaled";
        // this->means is vector<double>
        std::vector<float> means_autoscaled = morph::MathAlgo::autoscale (this->box_signal_means, 0.0, 1.0);
        df.add_contained_vals (dname.c_str(), means_autoscaled);

        dname = frameName + "/lmalign/centre_box_index";
        df.add_val (dname.c_str(), this->centre_lmaligned);
        dname = frameName + "/autoalign/centre_box_index";
        df.add_val (dname.c_str(), this->centre_autoaligned);

        // Save box_coords_pixels, box_coords_autoalign, box_coords_lmalign
        // Inflates the output files significantly. May be worth making this optional.
        if (this->savePerPixelData == true) {
            for (size_t bi = 0; bi < this->boxes_pixels.size(); ++bi) {

                dname = frameName + "/pixels/coords/boxes/box" + std::to_string(bi);
                df.add_contained_vals (dname.c_str(), this->box_coords_pixels[bi]);

                if (this->saveAutoAlignData == true) {
                    dname = frameName + "/autoalign/coords/boxes/box" + std::to_string(bi);
                    df.add_contained_vals (dname.c_str(), this->box_coords_autoalign[bi]);
                }

                if (this->saveLMAlignData == true) {
                    dname = frameName + "/lmalign/coords/boxes/box" + std::to_string(bi);
                    df.add_contained_vals (dname.c_str(), this->box_coords_lmalign[bi]);
                }

                if (this->saveAutoAlignData == true) {
                    // Plus also box_depths, where box_depths = (box_coords_* - fitted_*) . (unit normal_*)
                    std::vector<double> box_depth_autoalign (this->box_coords_pixels[bi].size());
                    for (size_t p_i = 0; p_i < box_depth_autoalign.size(); ++p_i) {
                        cv::Point2d diffvec = box_coords_autoalign[bi][p_i] - fitted_autoaligned[bi];
                        box_depth_autoalign[p_i] = diffvec.dot(this->normals[bi]); // Should have length 1, so ok
                    }
                    dname = frameName + "/autoalign/box_depth/box" + std::to_string(bi);
                    df.add_contained_vals (dname.c_str(), box_depth_autoalign);
                }
                // LM aligned
                if (this->saveLMAlignData == true) {
                    std::vector<double> box_depth_lmalign (this->box_coords_pixels[bi].size());
                    for (size_t p_i = 0; p_i < box_depth_lmalign.size(); ++p_i) {
                        cv::Point2d diffvec = box_coords_lmalign[bi][p_i] - fitted_lmaligned[bi];
                        box_depth_lmalign[p_i] = diffvec.dot(this->normals[bi]); // Should have length 1, so ok
                    }
                    dname = frameName + "/lmalign/box_depth/box" + std::to_string(bi);
                    df.add_contained_vals (dname.c_str(), box_depth_lmalign);
                }
            }

            // Transform FLE[i] according to autoalign and lmalign and save
            for (size_t i = 0; i<this->FLE.size(); ++i) {

                std::vector<cv::Point2d> FLE_scaled (this->FLE[i].size());
                if (this->saveAutoAlignData == true || this->saveLMAlignData == true) {
                    this->scalePoints (this->FLE[i], FLE_scaled);
                }

                if (this->saveAutoAlignData == true) {
                    std::string s1 = frameName + "/autoalign/coords/freehand/loop" + std::to_string(i);
                    std::vector<cv::Point2d> FLE_aa;
                    this->transform (FLE_scaled, FLE_aa, this->autoalign_translation, this->autoalign_theta);
                    df.add_contained_vals (s1.c_str(), FLE_aa);
                }

                if (this->saveLMAlignData == true) {
                    std::string s2 = frameName + "/lmalign/coords/freehand/loop" + std::to_string(i);
                    std::vector<cv::Point2d> FLE_lm;
                    this->transform (FLE_scaled, FLE_lm, this->lm_translation, this->lm_theta);
                    df.add_contained_vals (s2.c_str(), FLE_lm);
                }
            }
        }

        // Freehand loop areas
        for (size_t i = 0; i<this->FLE.size(); ++i) {
            std::string s3 = frameName + "/scaled/freehand/area" + std::to_string(i);
            // Area is area per pixel * num pixels
            double area = static_cast<double>(this->FLE[i].size()) / (this->pixels_per_mm * this->pixels_per_mm);
            df.add_val (s3.c_str(), area);
        }

        // Record the normal vectors
        dname = frameName + "/unit_normals";
        df.add_contained_vals (dname.c_str(), this->normals);

        // Freehand drawn regions - results
        for (size_t ri = 0; ri < this->FL_pixels.size(); ++ri) { // FL_pixels.size() is the outer container: num loops
            dname = frameName + "/signal/bits8/freehand/loop" + std::to_string(ri);
            df.add_contained_vals (dname.c_str(), this->FL_pixels[ri]);
            dname = frameName + "/signal/postproc/freehand/loop" + std::to_string(ri);
            df.add_contained_vals (dname.c_str(), this->FL_signal[ri]);
        }
        dname = frameName + "/nfreehand";
        df.add_val (dname.c_str(), static_cast<unsigned int>(this->FL_pixels.size()));

        dname = frameName + "/signal/postproc/freehand/means";
        df.add_contained_vals (dname.c_str(), this->FL_signal_means);
        dname = frameName + "/signal/bits8/freehand/means";
        df.add_contained_vals (dname.c_str(), this->FL_pixel_means);

        // Here, compute centroid of freehand regions, and then save this in autoalign
        // coordinate system and in the lmalign coord system.
        //
        // Add the centroid of the freehand regions (in the y-z or 'in-slice' plane)
        unsigned int fle_size = this->FLE.size();
        for (size_t i = 0; i<fle_size; ++i) {
            cv::Point cntroid = morph::MathAlgo::centroid (this->FLE[i]);

            std::cout << "centroid in screen pix: " << cntroid << std::endl;

            cv::Point2d cntroid_scaled = this->scalePoint (cntroid);
            cv::Point2d cntroid_autoaligned = this->transform (cntroid_scaled, this->autoalign_translation, this->autoalign_theta);
            cv::Point2d cntroid_lmaligned = this->transform (cntroid_scaled, this->lm_translation, this->lm_theta);

            if (this->saveAutoAlignData == true) {
                std::stringstream cntss1;
                cntss1 << frameName + "/autoalign/freehand/loop" << std::to_string(i) << "_centroid";
                df.add_contained_vals (cntss1.str().c_str(), cntroid_autoaligned);
            }

            if (this->saveLMAlignData == true) {
                std::stringstream cntss2;
                cntss2 << frameName + "/lmalign/freehand/loop" << std::to_string(i) << "_centroid";
                df.add_contained_vals (cntss2.str().c_str(), cntroid_lmaligned);
            }
        }

        // I'm storing all coordinates of the fitted points here.
        dname = frameName + "/pixels/fitted";
        df.add_contained_vals (dname.c_str(), this->fitted);

        // Also save the scaled version of the fitted points
        dname = frameName + "/scaled/fitted";
        df.add_contained_vals (dname.c_str(), this->fitted_scaled);

        if (this->saveAutoAlignData == true) {
            // These are the fitted points in the final coordinate system for the
            // auto-rotated (and transformed) slices
            dname = frameName + "/autoalign/fitted";
            df.add_contained_vals (dname.c_str(), this->fitted_autoaligned);
            dname = frameName + "/autoalign/computed";
            df.add_val (dname.c_str(), this->autoalignComputed);
        }

        if (this->saveLMAlignData == true) {
            // Fitted points where slices have been landmark aligned
            dname = frameName + "/lmalign/fitted";
            df.add_contained_vals (dname.c_str(), this->fitted_lmaligned);
            // Was the landmark alignment computed?
            dname = frameName + "/lmalign/computed";
            df.add_val (dname.c_str(), this->lmalignComputed);
        }

        if (this->saveAutoAlignData == true) {
            // The parameters of the translation which takes us from (fitted/LM)_scaled to (fitted/LM)_autoaligned
            dname = frameName + "/autoalign/translation";
            std::pair<double, double> aa_t = std::make_pair(this->autoalign_translation.x, this->autoalign_translation.y);
            df.add_contained_vals (dname.c_str(), aa_t);
            dname = frameName + "/autoalign/theta";
            df.add_val (dname.c_str(), this->autoalign_theta);
        }

        if (this->saveLMAlignData == true) {
            // The parameters of the translation which takes us from (fitted/LM)_scaled to (fitted/LM)_autoaligned
            dname = frameName + "/lmalign/translation";
            std::pair<double, double> lm_t = std::make_pair(this->lm_translation.x, this->lm_translation.y);
            df.add_contained_vals (dname.c_str(), lm_t);
            dname = frameName + "/lmalign/theta";
            df.add_val (dname.c_str(), this->lm_theta);
        }

        // Save autoalign and lmalign translated landmark coordinates.
        if (this->saveAutoAlignData == true) {
            std::vector<std::array<float,3>> LM_autoaligned_3d (this->LM_autoaligned.size(), {this->layer_x,0.0f,0.0f});
            for (size_t i = 0; i < this->LM_autoaligned.size(); ++i) {
                LM_autoaligned_3d[i][1] = static_cast<float>(this->LM_autoaligned[i].x);
                LM_autoaligned_3d[i][2] = static_cast<float>(this->LM_autoaligned[i].y);
            }
            dname = frameName + "/autoalign/landmarks";
            df.add_contained_vals (dname.c_str(), LM_autoaligned_3d);
        }

        if (this->saveLMAlignData == true) {
            std::vector<std::array<float,3>> LM_lmaligned_3d(this->LM_lmaligned.size(), {this->layer_x,0.0f,0.0f});
            for (size_t i = 0; i < this->LM_lmaligned.size(); ++i) {
                LM_lmaligned_3d[i][1] = static_cast<float>(this->LM_lmaligned[i].x);
                LM_lmaligned_3d[i][2] = static_cast<float>(this->LM_lmaligned[i].y);
            }
            dname = frameName + "/lmalign/landmarks";
            df.add_contained_vals (dname.c_str(), LM_lmaligned_3d);
        }

        // Save autoalign and lmalign translated global landmark coordinates.
        if (this->saveAutoAlignData == true) {
            std::vector<std::array<float,3>> GLM_autoaligned_3d (this->GLM_autoaligned.size(), {this->layer_x,0.0f,0.0f});
            for (size_t i = 0; i < this->GLM_autoaligned.size(); ++i) {
                GLM_autoaligned_3d[i][1] = static_cast<float>(this->GLM_autoaligned[i].x);
                GLM_autoaligned_3d[i][2] = static_cast<float>(this->GLM_autoaligned[i].y);
            }
            dname = frameName + "/autoalign/global_landmarks";
            df.add_contained_vals (dname.c_str(), GLM_autoaligned_3d);
        }

        if (this->saveLMAlignData == true) {
            std::vector<std::array<float,3>> GLM_lmaligned_3d(this->GLM_lmaligned.size(), {this->layer_x,0.0f,0.0f});
            for (size_t i = 0; i < this->GLM_lmaligned.size(); ++i) {
                GLM_lmaligned_3d[i][1] = static_cast<float>(this->GLM_lmaligned[i].x);
                GLM_lmaligned_3d[i][2] = static_cast<float>(this->GLM_lmaligned[i].y);
            }
            dname = frameName + "/lmalign/global_landmarks";
            df.add_contained_vals (dname.c_str(), GLM_lmaligned_3d);
        }

        // Need to get from fitted to y and z. Note that fitted is in (integer) pixels...
        // vector<cv::Point> fitted;
        //
        // Make up the boxes. A box (in 3d space) can be a vector of 12 floats. Thus
        // we should be able to write a vector of boxes as a vector<vector<float>>
        // These are "surface_boxes" because they're the box thats in the plane of the
        // cortical sheet (roughly xy) rather than the box in the slice plane (yz).
        std::array<float, 12> sbox;

        std::vector<std::array<float,12>> surface_boxes_autoaligned;
        std::vector<std::array<float,12>> surface_boxes_lmaligned;
        std::vector<std::array<float,12>> surface_boxes_scaled;
        std::vector<std::array<float,3>> surface_box_centroids_autoaligned;
        std::vector<std::array<float,3>> surface_box_centroids_lmaligned;
        std::vector<std::array<float,3>> surface_box_centroids_scaled;
        //std::cout << "Surface boxes extend from " << layer_x << " to " << (layer_x + thickness) << std::endl;
        for (int i = 1; i < this->nFit; ++i) {
            // Create auto-aligned surface box
            // c1 x,y,z
            sbox[0] = this->layer_x;                    // x
            sbox[1] = this->fitted_autoaligned[i-1].x;  // y
            sbox[2] = this->fitted_autoaligned[i-1].y;  // z
            // c2 x,y,z
            sbox[3] = this->layer_x;                    // x
            sbox[4] = this->fitted_autoaligned[i].x;    // y
            sbox[5] = this->fitted_autoaligned[i].y;    // z
            // c3 x,y,z
            sbox[6] = this->layer_x+this->thickness;    // x
            sbox[7] = this->fitted_autoaligned[i].x;    // y
            sbox[8] = this->fitted_autoaligned[i].y;    // z
            // c4 x,y,z
            sbox[9] = this->layer_x+this->thickness;    // x
            sbox[10] = this->fitted_autoaligned[i-1].x; // y
            sbox[11] = this->fitted_autoaligned[i-1].y; // z

            std::array<float, 3> sbox_centroid = morph::MathAlgo::centroid3D (sbox);//<float>
            surface_boxes_autoaligned.push_back (sbox);
            surface_box_centroids_autoaligned.push_back (sbox_centroid);

            // Create landmark-aligned surface box
            // c1 x,y,z
            sbox[1] = this->fitted_lmaligned[i-1].x;  // y
            sbox[2] = this->fitted_lmaligned[i-1].y;  // z
            // c2 x,y,z
            sbox[4] = this->fitted_lmaligned[i].x;    // y
            sbox[5] = this->fitted_lmaligned[i].y;    // z
            // c3 x,y,z
            sbox[7] = this->fitted_lmaligned[i].x;    // y
            sbox[8] = this->fitted_lmaligned[i].y;    // z
            // c4 x,y,z
            sbox[10] = this->fitted_lmaligned[i-1].x; // y
            sbox[11] = this->fitted_lmaligned[i-1].y; // z

            sbox_centroid = morph::MathAlgo::centroid3D (sbox);
            surface_boxes_lmaligned.push_back (sbox);
            surface_box_centroids_lmaligned.push_back (sbox_centroid);

            // Create un-aligned surface box for debugging
            // c1 x,y,z
            sbox[1] = this->fitted_scaled[i-1].x;  // y
            sbox[2] = this->fitted_scaled[i-1].y;  // z
            // c2 x,y,z
            sbox[4] = this->fitted_scaled[i].x;    // y
            sbox[5] = this->fitted_scaled[i].y;    // z
            // c3 x,y,z
            sbox[7] = this->fitted_scaled[i].x;    // y
            sbox[8] = this->fitted_scaled[i].y;    // z
            // c4 x,y,z
            sbox[10] = this->fitted_scaled[i-1].x; // y
            sbox[11] = this->fitted_scaled[i-1].y; // z

            sbox_centroid = morph::MathAlgo::centroid3D (sbox);
            surface_boxes_scaled.push_back (sbox);
            surface_box_centroids_scaled.push_back (sbox_centroid);
        }

        if (this->saveAutoAlignData == true) {
            // sboxes are 'surface boxes' - they lie in the plane of the cortical surface
            // and are not to be confused with the yellow boxes drawn in the UI in the y-z
            // plane.
            dname = frameName + "/autoalign/sboxes";
            df.add_contained_vals (dname.c_str(), surface_boxes_autoaligned);

            dname = frameName + "/autoalign/sbox_centers";
            df.add_contained_vals (dname.c_str(), surface_box_centroids_autoaligned);
        }

        if (this->saveLMAlignData == true) {
            dname = frameName + "/lmalign/sboxes";
            df.add_contained_vals (dname.c_str(), surface_boxes_lmaligned);

            dname = frameName + "/lmalign/sbox_centers";
            df.add_contained_vals (dname.c_str(), surface_box_centroids_lmaligned);
        }

        dname = frameName + "/scaled/sboxes";
        df.add_contained_vals (dname.c_str(), surface_boxes_scaled);

        dname = frameName + "/scaled/sbox_centers";
        df.add_contained_vals (dname.c_str(), surface_box_centroids_scaled);

        // From surface_box_centroids, can compute linear distance along curve. Could
        // be useful for making naive maps that unroll the cortex in one dimension.
        dname = frameName + "/scaled/flattened/sbox_linear_distance";
        {
            float total_linear_distance = 0.0f;
            std::vector<float> linear_distances (this->nBins, 0.0f);
            for (int i=1; i<this->nBins; ++i) {
                // Compute distance from Previous to current
                float d = morph::MathAlgo::distance<float> (surface_box_centroids_autoaligned[i-1],
                                                            surface_box_centroids_autoaligned[i]);
                total_linear_distance += d;
                linear_distances[i] = total_linear_distance;
            }
            // Now offset the linear distances so that the middle is 0.
            float halftotal = total_linear_distance / 2.0f;
            for (int i=0; i<this->nBins; ++i) {
                linear_distances[i] -= halftotal;
            }
            df.add_contained_vals (dname.c_str(), linear_distances);
        }

        // From 3D data, compute a 2D map. 3D data always centred around 0! Cool.
        // So, for each slice, compute each surface box's angle from the centroid and save these values.
        if (this->saveAutoAlignData == true && !this->fitted_autoaligned.empty()) {
            this->saveAutoAlignAngles (df, theta_middle, surface_box_centroids_autoaligned);
        }

        if (this->saveLMAlignData == true && !this->fitted_lmaligned.empty()) {
            this->saveLmAlignAngles (df, theta_middle, surface_box_centroids_lmaligned);
        }

        std::cout << "write() completed for " << frameName << std::endl;
    }

    //! Returns the cv::Point2d coordinate of the point on FrameData::fitted_lmaligned
    //! that has been deemed the 'middle' or 'centre' point for the curve.
    cv::Point2d getCentreLmaligned() const
    {
        cv::Point2d rtn(.0, .0);
        if (this->centre_lmaligned > -1 && static_cast<size_t>(this->centre_lmaligned) < this->fitted_lmaligned.size()) {
            rtn = this->fitted_lmaligned[this->centre_lmaligned];
        }
        return rtn;
    }

    //! Returns the cv::Point2d coordinate of the point on FrameData::fitted_autoaligned
    //! that has been deemed the 'middle' or 'centre' point for the curve.
    cv::Point2d getCentreAutoaligned() const
    {
        cv::Point2d rtn(.0, .0);
        if (this->centre_autoaligned > -1 && static_cast<size_t>(this->centre_autoaligned) < this->fitted_autoaligned.size()) {
            rtn = this->fitted_autoaligned[this->centre_autoaligned];
        }
        return rtn;
    }

    //! Given a centre point location that is obtained from \a otherframe, find the
    //! points on fitted_lmaligned and fitted_autoaligned which are closest and mark
    //! these as centre_lmaligned/centre_autoaligned.
    void setMiddle (const FrameData& otherframe)
    {
        double diff_min = 1e9;
        int i_min = 0;
        cv::Point2d opoint = otherframe.getCentreLmaligned();
        //std::cout << "setting middle from other point " << opoint << std::endl;
        for (int i=0; i<this->nBins; ++i) {
            cv::Point2d diff = fitted_lmaligned[i] - opoint;
            double diff_len = std::sqrt(diff.x*diff.x + diff.y*diff.y);
            if (diff_len < diff_min) {
                diff_min = diff_len;
                i_min = i;
            }
        }
        //std::cout << "set centre_lmaligned = " << i_min << std::endl;
        this->centre_lmaligned = i_min;

        diff_min = 1e9;
        i_min = 0;
        for (int i=0; i<this->nBins; ++i) {
            cv::Point2d diff = fitted_autoaligned[i] - otherframe.getCentreAutoaligned();
            double diff_len = std::sqrt(diff.x*diff.x + diff.y*diff.y);
            if (diff_len < diff_min) {
                diff_min = diff_len;
                i_min = i;
            }
        }
        this->centre_autoaligned = i_min;
    }

    //! Set a centre location in the yz plane using the given angle (in radians) for the
    //! given set of fitted points. Place result in \a middle_index.
    void setMiddle (double theta_middle,
                    std::vector<cv::Point2d>& fitted_points, int& middle_index,
                    const cv::Point2d& angle_origin)
    {
        //std::cout << "setMiddle(theta_middle=" << theta_middle<<", fitted_points, middle_index) called\n";

        // For landmark alignment, find the index of the box whose angle is closest to
        // theta_middle (or where one of its corners is closest). BUT make sure that if
        // there are two options, one closer than the other, we prefer the one that's
        // further away.
        int n_points = fitted_points.size();

        std::vector<std::pair<int, morph::Vector<double, 2>>> candidates;
        double max_cand_rad = -1e9;
        double min_cand_rad = 1e9;
        for (int i=0; i<n_points; ++i) {

            // Offset the point wrt the angle origin
            cv::Point2d pt = fitted_points[i] - angle_origin;
            // Find angle diff to neighbouring bin (d_to_n)
            int i_n = i>0 ? i-1 : i+1;
            cv::Point2d pt_n = fitted_points[i_n] - angle_origin;
            double angle = std::atan2 (pt.x, pt.y);
            double d_to_n = std::abs(std::atan2 (pt_n.x, pt_n.y) - angle);
            d_to_n = (d_to_n > morph::PI_x3_OVER_2_D) ? std::abs(d_to_n - morph::TWO_PI_D) : d_to_n;

            // Prevent the above from containing multiples of 2pi
            double radius = std::sqrt(pt.x * pt.x + pt.y * pt.y);
            double diff = std::abs (theta_middle - angle);
            diff = (diff > morph::PI_x3_OVER_2_D) ? std::abs(diff - morph::TWO_PI_D) : diff;

            if (diff < d_to_n) {
                // It's a candidate
                //std::cout << "Candidate point i = " << i << " diff = " << diff
                //          << " (diff to neigb i_n = " << i_n << ": " << d_to_n << "), radius = " << radius << std::endl;
                morph::Vector<double, 2> cand_diffradius = {diff, radius};
                std::pair<int, morph::Vector<double, 2>> candidate = std::make_pair(i, cand_diffradius);
                candidates.push_back (candidate);
                if (radius > max_cand_rad) { max_cand_rad = radius; }
                if (radius < min_cand_rad) { min_cand_rad = radius; }
            }
        }

#ifdef USE_MEAN_RADIUS_IDEA
        // This use of a mean radius worked ok for cortex but not so good for spiral hippocampus
        double mean_radius = (max_cand_rad + min_cand_rad) / 2.0;
        double diff_min = 1e9;
        int i_chosen = 0;
        for (std::pair<int, morph::Vector<double, 2>> candidate : candidates) {
            // Test if the candidates radius is greater than the mean:
            if (candidate.second[1] >= mean_radius) {
                // The find the best difference:
                if (candidate.second[0] < diff_min) {
                    diff_min = candidate.second[0];
                    i_chosen = candidate.first;
                }
            }
        }
#else
        // Instead, ANY candidate is good enough, because any candidate is closer to the
        // ideal than the neighbouring box. So here, can simply choose the furthest one.
        double rad_max = -1.0;
        int i_chosen = 0;
        for (std::pair<int, morph::Vector<double, 2>> candidate : candidates) {
            if (candidate.second[1] > rad_max) {
                rad_max = candidate.second[1];
                i_chosen = candidate.first;
            }
        }
#endif

        // The fit point at i_chosen is the one, mark it as such.
        //std::cout << "From theta_middle, setting middle_index to " << i_chosen << " (diff_min: " << diff_min << ")" <<  std::endl;
        middle_index = i_chosen;
    }

    //! Set a centre location in the yz plane (for the purpose of unwrapping the 3D
    //! image into 2D) using the given angle (in radians).
    void setMiddle (double theta_middle)
    {
        cv::Point2d ao(0,0);
        this->setMiddle (theta_middle, fitted_lmaligned, this->centre_lmaligned, ao);
        this->setMiddle (theta_middle, fitted_autoaligned, this->centre_autoaligned, ao);
    }

    void setMiddle (double theta_middle, const cv::Point2d& angle_origin_lm, const cv::Point2d& angle_origin_auto)
    {
        this->setMiddle (theta_middle, fitted_lmaligned, this->centre_lmaligned, angle_origin_lm);
        this->setMiddle (theta_middle, fitted_autoaligned, this->centre_autoaligned, angle_origin_auto);
    }

    //! Mirror the image and mark in the flags that it was mirrored
    void mirror()
    {
        this->mirror_image_only();
        this->flags.flip (Mirrored);
    }
    //! Carry out the actual mirroring operation on its own, leaving flags unchanged
    void mirror_image_only()
    {
        this->compute_mirror (this->frame);
        this->compute_mirror (this->frame_signal);
    }

    //! Flip the image & mark as such in flags
    void flip()
    {
        this->flip_image_only();
        this->flags.flip (Flipped);
    }
    //! Flip the image without marking as flipped in flags.
    void flip_image_only()
    {
        this->compute_flip (this->frame);
        this->compute_flip (this->frame_signal);
    }

    //! Check landmarks on each slice. If there are at least n landmarks on every slice,
    //! and n >= 1, return n, else return 0 (or perhaps -n where n is the median number
    //! of landmarks).
    int landmarkCheck()
    {
        // Find minimum number of landmarks on any one slice
        int n = 100000;
        for (auto fr : (*this->parentStack)) {
            int _n = static_cast<int>(fr.LM.size());
            n = (_n < n) ? _n : n;
        }
        // Perhaps there were no frames, in which case, had better return 0:
        if (n == 100000) { n = 0; }
        return n;
    }

private:
    void updateAutoAlignments()
    {
        std::cout << __FUNCTION__ << " called\n";
        if (this->previous < 0) {
            // get centroid of this->fitted_scaled (the 0th slice
            cv::Point2d slice0centroid = morph::MathAlgo::centroid (this->fitted_scaled);
            // The translation to record for the first slice is the inverse of the centroid
            this->autoalign_translation = -slice0centroid;
            this->autoalign_theta = 0.0;

        } else {
            // now, if we're aligning partly to the mid slice, shouldn't we centroid it
            // first? Perhaps that's why the computeSos3d function is heavily weighted
            // to the neighbouring (i.e. to the previous) slice?
            size_t alignment_slice = this->parentStack->size()/2; // align to mid slice
            std::cout << "Alignment slice: " << alignment_slice << std::endl;

            // set number of bins in the alignment_slice and in the previous slice to
            // match the number of bins in this frame.
            int alignment_slice_bins = (*this->parentStack)[alignment_slice].getBins();
            int previous_slice_bins = (*this->parentStack)[this->previous].getBins();
            if (alignment_slice_bins != this->nBins) {
                (*this->parentStack)[alignment_slice].setBins (this->nBins);
                (*this->parentStack)[alignment_slice].updateFit();
            }
            if (previous_slice_bins != this->nBins) {
                (*this->parentStack)[this->previous].setBins (this->nBins);
                (*this->parentStack)[this->previous].updateFit();
            }

            // Call the optimization function
            this->alignOptimally (this->fitted_scaled,
                                  (*this->parentStack)[alignment_slice].fitted_autoaligned,
                                  (*this->parentStack)[this->previous].fitted_autoaligned,
                                  this->autoalign_translation, this->autoalign_theta);

            // Restore number of bins.
            if (alignment_slice_bins != this->nBins) {
                (*this->parentStack)[alignment_slice].setBins (alignment_slice_bins);
                (*this->parentStack)[alignment_slice].updateFit();
            }
            if (previous_slice_bins != this->nBins) {
                (*this->parentStack)[this->previous].setBins (previous_slice_bins);
                (*this->parentStack)[this->previous].updateFit();
            }
        }
        this->transform (this->fitted_scaled, this->fitted_autoaligned, this->autoalign_translation, this->autoalign_theta);
        this->transform (this->LM_scaled, this->LM_autoaligned, this->autoalign_translation, this->autoalign_theta);
        this->transform (this->AM_scaled, this->AM_autoaligned, this->autoalign_translation, this->autoalign_theta);
        this->transform (this->GLM_scaled, this->GLM_autoaligned, this->autoalign_translation, this->autoalign_theta);

        this->autoalignComputed = true;
    }

private:
    void landmarkAlignByRotation()
    {
        // First job is to translate landmark 1 of the slice so it's on the origin of
        // the slice plane. That's just the inverse of the location of (the scaled
        // version of) landmark 1.
        this->lm_translation = -this->LM_scaled[0];

        // If there's no previous frame, there's no rotation to apply - we'll rotate to match this one
        if (this->previous < 0) {
            this->lm_theta = 0.0;

            // get centroid of this->fitted_scaled (the 0th slice)
            cv::Point2d slice0centroid = morph::MathAlgo::centroid (this->fitted_scaled);
            this->curve_centroid = LM_scaled[0] - slice0centroid;

        } else {

            // Translate the points ready for the rotational optimization
            std::vector<cv::Point2d> LM_translated_unrotated (this->LM_scaled);
            this->transform (this->LM_scaled, LM_translated_unrotated, this->lm_translation, 0.0);

            size_t alignment_slice = 0;

            if (this->rotateButAlignLandmarkTwoPlus == true) {
                // alignment coords are rest of the landmarks
                this->lm_theta = this->alignOptimallyByRotation (LM_translated_unrotated,
                                                                 (*this->parentStack)[alignment_slice].LM_lmaligned,
                                                                 (*this->parentStack)[this->previous].LM_lmaligned);

            } else {
                // Alignment against curve points
                this->lm_theta = this->alignOptimallyByRotation (LM_translated_unrotated,
                                                                 (*this->parentStack)[alignment_slice].fitted_lmaligned,
                                                                 (*this->parentStack)[this->previous].fitted_lmaligned);
            }

        }

        // Apply the rotation.
        // The transform includes the lm_translation of the very first slice, so that
        // the centroid of this slice is roughly 0? Or better to use the actual centroid
        // of the first slice.
        this->transform (this->LM_scaled, this->LM_lmaligned,
                         this->lm_translation + (*this->parentStack)[0].curve_centroid,
                         this->lm_theta);

        this->transform (this->AM_scaled, this->AM_lmaligned,
                         this->lm_translation + (*this->parentStack)[0].curve_centroid,
                         this->lm_theta);

        this->transform (this->GLM_scaled, this->GLM_lmaligned,
                         this->lm_translation + (*this->parentStack)[0].curve_centroid,
                         this->lm_theta);

        this->transform (this->fitted_scaled, this->fitted_lmaligned,
                         this->lm_translation + (*this->parentStack)[0].curve_centroid,
                         this->lm_theta);

        this->lmalignComputed = true;
    }

private:
    void updateLandmarkAlignments()
    {
        std::cout << __FUNCTION__ << " called\n";
        // Landmark scheme, if we have >=2 landmarks on each slice and same number of
        // landmarks on each slice the we can compute alignment with the landmarks.

        // If 1 landmark, then do an
        // "align-on-the-landmark-and-rotate-optimally-thereafter" process. This could
        // either rotate to align landmarks 2 and above, but what it does at the moment
        // is rotate to align the curves.

        if (this->landmarkCheck() == 1) {
            // align-by-single-landmark and rotate to fit
            std::cout << "align by rotation...\n";
            this->landmarkAlignByRotation();
            return;

        } else if (this->landmarkCheck() < 2) { // 0. If 2
            return;
        } else {
            // 2 or higher. See if configuration mandates rotate-about-landmark-1
            if (this->rotateLandmarkOne == true) {
                this->landmarkAlignByRotation();
                return;
            }
        }

        // If there's no previous frame, we just translate the frame to lie around the origin (with no rotation)
        if (this->previous < 0) {
            // get centroid of this->fitted_scaled (the 0th slice
            cv::Point2d slice0centroid = morph::MathAlgo::centroid (this->fitted_scaled);
            // The translation to record for the first slice is the inverse of the centroid
            this->lm_translation = -slice0centroid;
            this->lm_theta = 0.0;

        } else {
            // NB: This is the landmark scheme
            size_t alignment_slice = 0; // Align to the 0th slice

            // Unify the number of sample bins on the curves in the alignment_slice frame and in this frame.
            int alignment_slice_bins = (*this->parentStack)[alignment_slice].getBins();
            if (alignment_slice_bins != this->nBins) {
                (*this->parentStack)[alignment_slice].setBins (this->nBins);
                (*this->parentStack)[alignment_slice].updateFit();
            }

            this->alignOptimally (this->LM_scaled,
                                  (*this->parentStack)[alignment_slice].LM_lmaligned,
                                  //(*this->parentStack)[this->previous].LM_lmaligned, // keep landmark aligning simple
                                  this->lm_translation, this->lm_theta);

            // Restore number of bins.
            if (alignment_slice_bins != this->nBins) {
                (*this->parentStack)[alignment_slice].setBins (alignment_slice_bins);
                (*this->parentStack)[alignment_slice].updateFit();
            }
        }

        this->transform (this->LM_scaled, this->LM_lmaligned, this->lm_translation, this->lm_theta);
        std::cout << "Transform AM_scaled size " << this->AM_scaled.size() << " into this->AM_lmaligned size " << this->AM_lmaligned.size() << std::endl;
        this->transform (this->AM_scaled, this->AM_lmaligned, this->lm_translation, this->lm_theta);
        this->transform (this->GLM_scaled, this->GLM_lmaligned, this->lm_translation, this->lm_theta);
        this->transform (this->fitted_scaled, this->fitted_lmaligned, this->lm_translation, this->lm_theta);

        this->lmalignComputed = true;
    }

public:
    //! Compute the align-centroid-and-rotate slice alignments and if possible, the
    //! landmark-based alignment.
    void updateAlignments()
    {
        if (this->fitted.empty()) { return; }

        this->scalePoints (this->fitted, this->fitted_scaled);

        // Need to resize the container for the autoalign-translated copies of any landmarks:
        // Re-size landmark containers.
        if (this->landmarkCheck() > 0) {
            std::cout << "LM (landmark container) size: " << this->LM.size() << std::endl;
            this->LM_scaled.resize(this->LM.size());
            this->LM_lmaligned.resize(this->LM.size());
            this->LM_autoalign_translated.resize(this->LM.size());
            this->LM_autoaligned.resize(this->LM.size());

            // Scale the landmarks
            this->scalePoints (this->LM, this->LM_scaled);
        }

        if (!this->AM.empty()) {
            this->AM_scaled.resize(this->AM.size());
            this->AM_lmaligned.resize(this->AM.size());
            this->AM_autoalign_translated.resize(this->AM.size());
            this->AM_autoaligned.resize(this->AM.size());
            std::cout << "Scaling AM into AM_scaled...\n";
            this->scalePoints (this->AM, this->AM_scaled);
            if (!this->AM.empty()) {
                std::cout << "AM[0]: " << AM[0] << std::endl;
                std::cout << "AM_scaled[0]: " << AM_scaled[0] << std::endl;
            }
        }

        if (!this->GLM.empty()) {
            this->GLM_scaled.resize(this->GLM.size());
            this->GLM_lmaligned.resize(this->GLM.size());
            this->GLM_autoalign_translated.resize(this->GLM.size());
            this->GLM_autoaligned.resize(this->GLM.size());
            std::cout << "Scaling GLM into GLM_scaled...\n";
            this->scalePoints (this->GLM, this->GLM_scaled);
            if (!this->GLM.empty()) {
                std::cout << "GLM[0]: " << GLM[0] << std::endl;
                std::cout << "GLM_scaled[0]: " << GLM_scaled[0] << std::endl;
            }
        }

        // Set variables saying aligned_with_centroids = false; aligned_with_landmarks =
        // false; for writing into the h5 file. Then set these true as appropriate,
        // below.
        this->autoalignComputed = false;
        this->lmalignComputed = false;

        this->updateAutoAlignments();
        this->updateLandmarkAlignments();
    }

    //! Public wrapper around updateFitBezier()
    void updateFit() { this->updateFitBezier(); }

    //! Re-compute the boxes from the curve (taking ints)
    void refreshBoxes (const int lenA, const int lenB) { this->refreshBoxes ((double)lenA, (double)lenB); }

    //! Re-compute the boxes from the curve (double version)
    void refreshBoxes (const double lenA, const double lenB)
    {
        // Don't refresh boxes for Freehand mode or if this->fitted is empty
        if (this->ct == InputMode::Freehand) { return; }
        if (this->fitted.empty()) { return; }

        this->pointsInner.resize(this->nFit);
        this->pointsOuter.resize(this->nFit);
        for (int i=0; i<this->nFit; i++) {
            cv::Point2d normLenA = this->normals[i]*lenA;
            cv::Point2d normLenB = this->normals[i]*lenB;
            this->pointsInner[i] = cv::Point(this->fitted[i] + cv::Point2d(normLenA.x, normLenA.y));
            this->pointsOuter[i] = cv::Point(this->fitted[i] + cv::Point2d(normLenB.x, normLenB.y));
        }

        // Make the boxes from pointsInner and pointsOuter
        this->boxes.resize (this->nBins);
        for (int i=0; i<this->nBins; i++) {
            std::vector<cv::Point> pts(4);
            pts[0] = this->pointsInner[i];
            pts[1] = this->pointsInner[i+1];
            pts[2] = this->pointsOuter[i+1];
            pts[3] = this->pointsOuter[i];
            this->boxes[i] = pts;
        }
    }

    //! Clear the boxes, making pointsInner and pointsOuter all 0
    void clearBoxes()
    {
        this->pointsInner.clear();
        this->pointsOuter.clear();
        this->boxes.clear();
    }

    //! For each freehand drawn loop, compute the mean luminance within the loop,
    //! storing in this->FL_signal_means
    void computeFreehandMeans()
    {
        // Loop through FLE. For each set of points, output the points as a list and
        // also compute the mean.
        this->FL_pixel_means.resize (this->FLE.size());
        this->FL_signal_means.resize (this->FLE.size());
        this->FL_pixels.resize (this->FLE.size());
        this->FL_signal.resize (this->FLE.size());
        this->FL_pixels_bgr.resize (this->FLE.size());
        for (size_t i=0; i<this->FLE.size(); i++) {
            // region is FLE[i]
            this->FL_signal_means[i] = 0.0;
            this->FL_pixel_means[i] = 0;

            this->FL_pixels[i] = this->getRegionPixelVals (this->FLE[i]);
            this->FL_signal[i] = this->getRegionSignalVals (this->FLE[i]);
            morph::MathAlgo::compute_mean_sd<unsigned int> (this->FL_pixels[i], FL_pixel_means[i]);
            morph::MathAlgo::compute_mean_sd<float> (this->FL_signal[i], FL_signal_means[i]);
        }
    }

    // Toggle controls

    void toggleShowBoxes() { this->flags[ShowBoxes] = this->flags.test(ShowBoxes) ? false : true; }
    void setShowBoxes (bool t) { this->flags[ShowBoxes] = t; }
    bool getShowBoxes() { return this->flags.test(ShowBoxes); }

    void toggleShowFits() { this->flags[ShowFits] = this->flags.test(ShowFits) ? false : true; }
    void setShowFits (bool t) { this->flags[ShowFits] = t; }
    bool getShowFits() { return this->flags.test(ShowFits); }

    void toggleShowUsers() { this->flags[ShowUsers] = this->flags.test(ShowUsers) ? false : true; }
    void setShowUsers (bool t) { this->flags[ShowUsers] = t; }
    bool getShowUsers() { return this->flags.test(ShowUsers); }

    void toggleShowCtrls() { this->flags[ShowCtrls] = this->flags.test(ShowCtrls) ? false : true; }
    void setShowCtrls (bool t) { this->flags[ShowCtrls] = t; }
    bool getShowCtrls() { return this->flags.test(ShowCtrls); }

private:

    //! Mirror cv::Mat \a m about the vertical axis
    void compute_mirror (cv::Mat& m)
    {
        cv::Mat mirrored (m.rows, m.cols, m.type());
        cv::flip (m, mirrored, 1);
        m = mirrored;
    }

    //! Flip cv::Mat \a m about the horizontal axis
    void compute_flip (cv::Mat& m)
    {
        cv::Mat flipped (m.rows, m.cols, m.type());
        cv::flip (m, flipped, 1);
        m = flipped;
    }

    //! Get the signal values from the region from the mRNA signal window (frame_signal, red channel).
    std::vector<float> getRegionSignalVals (const std::vector<cv::Point>& region)
    {
        std::vector<float> regionSignalVals (region.size(), 0.0f);
        size_t i = 0;
        for (auto px : region) {
            cv::Vec<float, 3> val = this->frame_signal.at<cv::Vec<float, 3>>(px.y, px.x);
            regionSignalVals[i++] = static_cast<float>(val[0]);
        }
        return regionSignalVals;
    }

    //! Get the pixel values from the region from the original image window (frame, red channel).
    std::vector<unsigned int> getRegionPixelVals (const std::vector<cv::Point>& region)
    {
        std::vector<unsigned int> regionPixelVals (region.size(), 0);
        size_t i = 0;
        for (auto px : region) {
            int _col = px.x;
            int _row = px.y;
            cv::Vec<unsigned char, 3> pixelval = this->frame.at<cv::Vec<unsigned char, 3>>(_row, _col);
            regionPixelVals[i++] = static_cast<unsigned int>(pixelval[0]); // Result between 0 and 255
        }
        return regionPixelVals;
    }

    std::vector<cv::Point> getBoxRegion (const std::vector<cv::Point> boxVtxs)
    {
        // Four cv::Points define a rectangle. Convert from vector to an array of Points
        cv::Point pts[4] = {boxVtxs[0],boxVtxs[1],boxVtxs[2],boxVtxs[3]};
        // Create a mask of (initially) zeros.
        cv::Mat mask = cv::Mat::zeros(this->frame.rows, this->frame.cols, CV_8UC3);
        // Set the box defined by pts to ones.
        cv::fillConvexPoly (mask, pts, 4, cv::Scalar(255,255,255));
        // Make a grayscale version of mask.
        cv::Mat maskGray;
        cv::cvtColor (mask, maskGray, cv::COLOR_BGR2GRAY);
        // Find the nonzero locations in maskGray.
        std::vector<cv::Point> maskpositives;
        cv::findNonZero (maskGray, maskpositives);
        return maskpositives;
    }

    //! Get the raw pixel values in the box defined by boxVtxs. \a maskpositives: Output
    //! the coordinates (in pixels) of all the non-zero pixels
    std::vector<unsigned int> getBoxedPixelVals (const std::vector<cv::Point> boxVtxs,
                                                 std::vector<cv::Point>& maskpositives)
    {
        maskpositives = this->getBoxRegion (boxVtxs);
        // Save maskpositives into the relevant box? These are the pixel coords for box[i], whatever i is.
        return this->getRegionPixelVals (maskpositives);
    }

    //! Find a grayscale value for each pixel of image in FrameData::frame_bgoff within
    //! the box defined by \a boxVtxs and return this in a vector of floats, without
    //! conversion
    std::vector<float> getBoxedSignalVals (const std::vector<cv::Point> boxVtxs)
    {
        std::vector<cv::Point> maskpositives = this->getBoxRegion (boxVtxs);
        return this->getRegionSignalVals (maskpositives);
    }

    //! Compute the mean values for the bins. Not const. But means don't need to be a
    //! member as they're only computed to be written out to file.
    void computeBoxMeans()
    {
        this->boxes_pixels.resize (this->boxes.size());
        this->boxes_signal.resize (this->boxes.size());
        this->box_coords_pixels.resize (this->boxes.size());
        this->box_coords_autoalign.resize (this->boxes.size());
        this->box_coords_lmalign.resize (this->boxes.size());
        this->boxes_pixels_bgr.resize (this->boxes.size());
        this->box_signal_means.resize (this->boxes.size());
        this->box_pixel_means.resize (this->boxes.size());
        this->box_signal_sds.resize (this->boxes.size());
        this->box_pixel_sds.resize (this->boxes.size());

        for (size_t i=0; i<this->boxes.size(); i++) {

            // Zero the means value
            this->box_pixel_means[i] = 0;
            this->box_signal_means[i] = 0.0;
            this->box_pixel_sds[i] = 0;
            this->box_signal_sds[i] = 0.0;

            // We get the box pixel VALUES in the return value and the COORDS in the output argument
            this->boxes_pixels[i] = this->getBoxedPixelVals (this->boxes[i], this->box_coords_pixels[i]);
            this->boxes_signal[i] = this->getBoxedSignalVals (this->boxes[i]);
            this->box_pixel_sds[i] = morph::MathAlgo::compute_mean_sd<unsigned int> (this->boxes_pixels[i], this->box_pixel_means[i]);
            this->box_signal_sds[i] = morph::MathAlgo::compute_mean_sd<float> (this->boxes_signal[i], this->box_signal_means[i]);
            // transform box_coords_pixels into box_coords_autoalign and box_coords_lmalign
            std::vector<cv::Point2d> bcpix(this->box_coords_pixels[i].size());
            this->scalePoints (this->box_coords_pixels[i], bcpix);
            this->transform (bcpix, this->box_coords_autoalign[i], this->autoalign_translation, this->autoalign_theta);
            this->transform (bcpix, this->box_coords_lmalign[i], this->lm_translation, this->lm_theta);
        }
    }

    //! Update the fit, scale and rotate by \a _theta.
    void updateFit (double _theta)
    {
        this->updateFitBezier();
        this->scalePoints (this->fitted, this->fitted_scaled);
        this->autoalign_translation = -morph::MathAlgo::centroid (this->fitted_scaled);
        this->translate (this->fitted_scaled, this->fitted_autoalign_translated, this->autoalign_translation);
        this->rotate (this->fitted_autoalign_translated, this->fitted_autoaligned, _theta);
    }

    //! Recompute the Bezier fit
    void updateFitBezier()
    {
        if (this->PP.empty()) {
            //std::cout << "Too few points to fit" << std::endl;
            this->bcp.reset();
            this->fitted.clear();
            this->tangents.clear();
            this->normals.clear();
            this->clearBoxes();
            return;
        }

        this->bcp.reset();

        // Loop over PP first
        for (auto _P : this->PP) {
            std::vector<std::pair<double,double>> user_points;
            user_points.clear();
            for (auto pt : _P) { user_points.push_back (std::make_pair(pt.x, pt.y)); }

            morph::BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                this->bcp.addCurve (bc);
            } else {
                // Have previous curve, use last control of previous curve to make
                // smooth transition.
                morph::BezCurve<double> last = this->bcp.curves.back();
                bc.fit (user_points, last);
                this->bcp.removeCurve();
                this->bcp.addCurve (last);
                this->bcp.addCurve (bc);
            }
        }

        // Update this->fitted
        this->bcp.computePoints (static_cast<unsigned int>(this->nFit));
        std::vector<morph::BezCoord<double>> coords = this->bcp.getPoints();
        // tangents and normals coming out of bcp are already unit-normalized.
        std::vector<morph::BezCoord<double>> tans = this->bcp.getTangents();
        std::vector<morph::BezCoord<double>> norms = this->bcp.getNormals();
        // Point2d fitsum;
        this->fitted.resize (this->nFit);
        this->tangents.resize (this->nFit);
        this->normals.resize (this->nFit);
        for (int i = 0; i < this->nFit; ++i) {
            this->fitted[i] = cv::Point2d(coords[i].x(),coords[i].y());
            this->tangents[i] = cv::Point2d(tans[i].x(),tans[i].y());
            this->normals[i] = cv::Point2d(norms[i].x(),norms[i].y());
        }
    }

    //! Scale \a points by FrameData::pixels_per_mm to give \a scaled_points
    void scalePoints (const std::vector<cv::Point>& points, std::vector<cv::Point2d>& scaled_points)
    {
        size_t sz = points.size();
        size_t sz2 = scaled_points.size();
        if (sz != sz2) { throw std::runtime_error ("scalePoints: size mismatch"); }
        for (size_t i = 0; i < sz; ++i) {
            scaled_points[i] = cv::Point2d(points[i]) / this->pixels_per_mm;
        }
    }

    //! Scale \a points by FrameData::pixels_per_mm to give \a scaled_points
    void scalePoints (const std::vector<cv::Point2d>& points, std::vector<cv::Point2d>& scaled_points)
    {
        size_t sz = points.size();
        size_t sz2 = scaled_points.size();
        if (sz != sz2) { throw std::runtime_error ("scalePoints: size mismatch"); }
        for (size_t i = 0; i < sz; ++i) {
            scaled_points[i] = points[i] / this->pixels_per_mm;
        }
    }

    //! Scale the coordinate \a pt by FrameData::pixels_per_mm and return the result
    cv::Point2d scalePoint (const cv::Point& pt)
    {
        cv::Point2d rtn = cv::Point2d(pt)/this->pixels_per_mm;
        return rtn;
    }

    //! Rotation only
    double computeSos1d (const std::vector<cv::Point2d> target_coords,
                         const std::vector<cv::Point2d> coords,
                         const double rotn)
    {
        std::vector<double> vtx3d = {0.0, 0.0, rotn};
        return this->computeSos3d (target_coords, coords, vtx3d);
    }

    /*!
     * vertex contains x,y,theta values and should have size 3 translate coords by
     * vertex[0],vertex[1] and rotate by vertex[2]. Compute SOS compared with target
     * coordinates
     */
    double computeSos3d (const std::vector<cv::Point2d> target_coords,
                         const std::vector<cv::Point2d> coords,
                         const std::vector<double>& vertex)
    {
        if (this->previous < 0) { return std::numeric_limits<double>::max(); }

        if (coords.size() != target_coords.size()) {
            throw std::runtime_error ("computeSos3d: Number of elements in coords and target_coords must be the same");
        }

        //std::cout << "computeSos3d: x,y,theta = (" << vertex[0] << "," << vertex[1] << "," << vertex[2] << ")" << std::endl;

        // Transform coords
        std::vector<cv::Point2d> tmp_coords (coords);
        cv::Point2d tr(vertex[0], vertex[1]);
        this->translate (coords, tmp_coords, tr);
        this->rotate (tmp_coords, tmp_coords, vertex[2]);

        // Compute SOS wrt to target_coords
        double sos = 0.0;
        for (size_t i = 0; i < target_coords.size(); ++i) {
            double xdiff = tmp_coords[i].x - target_coords[i].x;
            double ydiff = tmp_coords[i].y - target_coords[i].y;
            double d_ = (xdiff*xdiff + ydiff*ydiff);
            sos += d_;
        }

        //std::cout << "sos=" << sos << std::endl;
        return sos;
    }

    //! Rotation only
    double computeSos1d (const std::vector<cv::Point2d> target_coords,
                         const std::vector<cv::Point2d> neighbour_coords,
                         const std::vector<cv::Point2d> coords,
                         const double rotn)
    {
        std::vector<double> vtx3d = {0.0, 0.0, rotn};
        return this->computeSos3d (target_coords, neighbour_coords, coords, vtx3d);
    }

    /*!
     * vertex contains x,y,theta values and should have size 3 translate coords by
     * vertex[0],vertex[1] and rotate by vertex[2]. Compute SOS compared with target
     * coordinates AND neighbour coords.
     */
    double computeSos3d (const std::vector<cv::Point2d> target_coords,
                         const std::vector<cv::Point2d> neighbour_coords,
                         const std::vector<cv::Point2d> coords,
                         const std::vector<double>& vertex)
    {
        if (this->previous < 0) { return std::numeric_limits<double>::max(); }

        if (coords.size() != target_coords.size()) {
            throw std::runtime_error ("computeSos3d: Number of elements in coords and target_coords must be the same");
        }

        // Transform coords
        std::vector<cv::Point2d> tmp_coords (coords);
        cv::Point2d tr(vertex[0], vertex[1]);
        this->translate (coords, tmp_coords, tr);
        this->rotate (tmp_coords, tmp_coords, vertex[2]);

        // Alternative: take centroid of target_coords and use centroid of this coord
        // Compute SOS wrt to target_coords
        double cost = 0.0;

        // Output debug messages about cost?
        constexpr bool debug_cost = false;

        // This computes sos for coords wrt target coords (probably the first slice, but
        // maybe the central slice).
        for (size_t i = 0; i < target_coords.size(); ++i) {
            double xdiff = tmp_coords[i].x - target_coords[i].x;
            double ydiff = tmp_coords[i].y - target_coords[i].y;
            double d_ = (xdiff*xdiff + ydiff*ydiff);
            cost += d_;
        }
        cost /= (double)target_coords.size();

        // Weighting this by a number determined from looking at relative contribution
        // to cost
        double w_targ = 0.01;
        cost *= w_targ;

        if constexpr (debug_cost==true) { std::cout << "Costs. AlignTarg: " << cost; }

        size_t sz = neighbour_coords.size();
        double sos_  = 0.0;
        for (size_t i = 0; i < sz; ++i) {
            double xdiff = tmp_coords[i].x - neighbour_coords[i].x;
            double ydiff = tmp_coords[i].y - neighbour_coords[i].y;
            double d_ = xdiff*xdiff + ydiff*ydiff;
            sos_ += d_;
        }
        if constexpr (debug_cost==true) { std::cout << " AlignNeighbour: " << (sos_/(double)sz); }

        // Weight for aligning vs. neighbour coordinates
        double w_neigh = 1;
        cost += w_neigh * sos_ / (double)(sz);

        // If theta is big (in either dirn), then penalize this.
        double w_thet = 0.1;
        // Parameters of sigmoid cause ramp in cost around 0.2 rads
        double kxminusx0 = 30.0 * (std::abs(vertex[2]) - 0.2);
        double sigthet = (w_thet / (1.0 + std::exp(-kxminusx0)));

        if constexpr (debug_cost==true) { std::cout << " Angle: " << sigthet << std::endl; }

        cost += sigthet;

        return cost;
    }

    //! Compute a sum of squared distances between the points in this fit and the
    //! points in the previous fit
    double computeSosWithPrev (double _theta)
    {
        if (this->previous < 0) { return std::numeric_limits<double>::max(); }

        if ((*this->parentStack)[this->previous].getBins() != this->nBins) {
            // Number of bins has to be same
            return std::numeric_limits<double>::max();
        }

        this->rotate (this->fitted_autoalign_translated, this->fitted_autoaligned, _theta);

        double sos = 0.0;
        for (int i = 0; i < this->nFit; ++i) {
            // Distance^2 from this->fitted_autoaligned[i] to this->previous->fitted_autoaligned[i]
            double xdiff = this->fitted_autoaligned[i].x - (*this->parentStack)[this->previous].fitted_autoaligned[i].x;
            double ydiff = this->fitted_autoaligned[i].y - (*this->parentStack)[this->previous].fitted_autoaligned[i].y;
            double d_ = (xdiff*xdiff + ydiff*ydiff);
            //std::cout << "sos += " << d_ << ", ";
            sos += d_;
        }
        //std::cout << "\nFor rotation angle " << _theta << " returning sos=" << sos << std::endl;
        return sos;
    }

    //! Transform a single coordinate by a translation and a rotation
    cv::Point2d transform (const cv::Point2d& coord, const cv::Point2d& translation, const double _theta)
    {
        cv::Point2d trans = coord + translation;
        return this->rotate (trans, _theta);
    }

    //! Transform a vector of coordinates
    void transform (const std::vector<cv::Point2d>& coords, std::vector<cv::Point2d>& transformed,
                    const cv::Point2d& translation, const double _theta)
    {
        if (coords.size() != transformed.size()) { transformed.resize (coords.size()); }
        std::vector<cv::Point2d> translated (coords.size());
        this->translate (coords, translated, translation);
        this->rotate (translated, transformed, _theta);
    }

    //! Do a translation operation
    void translate (const std::vector<cv::Point2d>& coords, std::vector<cv::Point2d>& translated,
                    const cv::Point2d& translation)
    {
        if (coords.size() != translated.size()) {
            throw std::runtime_error ("translate(): input and output vectors must have same size");
        }
        for (size_t i = 0; i < coords.size(); ++i) {
            translated[i] = coords[i] + translation;
        }
    }

    //! Rotate the points in this->fitted_autoalign_translated by theta about the origin and store
    //! the result in this->fitted_autoaligned
    void rotate (const std::vector<cv::Point2d>& coords, std::vector<cv::Point2d>& rotated, const double _theta)
    {
        if (coords.size() != rotated.size()) {
            throw std::runtime_error ("rotate(): input and output vectors must have same size");
        }

        if (_theta == 0.0) {
            for (size_t i = 0; i < coords.size(); ++i) {
                rotated[i] = coords[i];
            }
            return;
        }

        for (size_t i = 0; i < coords.size(); ++i) {
            double xi = coords[i].x;
            double yi = coords[i].y;
            double sin_theta = sin (_theta);
            double cos_theta = cos (_theta);
            rotated[i].x = xi * cos_theta - yi * sin_theta;
            rotated[i].y = xi * sin_theta + yi * cos_theta;
        }
    }

    //! Rotate the point \a pt by theta about the origin, returning rotated point
    cv::Point2d rotate (const cv::Point2d& pt, const double _theta)
    {
        cv::Point2d rtn = pt;
        if (_theta == 0.0) { return rtn; }

        double xi = pt.x;
        double yi = pt.y;
        double sin_theta = sin (_theta);
        double cos_theta = cos (_theta);
        rtn.x = xi * cos_theta - yi * sin_theta;
        rtn.y = xi * sin_theta + yi * cos_theta;

        return rtn;
    }

    /*!
     * Translation is the translation to apply to alignment_coords prior to aligning
     * with target_aligned_alignment_coords. Does that make sense?
     *
     * By applying only rotations about the origin, find the best rotation to apply to
     * alignment_coords so that they line up with some combination of
     * target_aligned_alignment_coords adn neighbour_aligned_alignment_coords
     *
     * \return The rotation that best lines the points up.
     */
    double alignOptimallyByRotation (const std::vector<cv::Point2d>& alignment_coords,
                                     const std::vector<cv::Point2d>& target_aligned_alignment_coords,
                                     const std::vector<cv::Point2d>& neighbour_aligned_alignment_coords)
    {
         // Check that target frame has same number of coords, if not warn and return.
        if (target_aligned_alignment_coords.size() != alignment_coords.size()
            || neighbour_aligned_alignment_coords.size() != alignment_coords.size()) {
            std::stringstream ee;
            ee << "alignOptimallyByRotation: Same number of alignment coords in both frames please! target_aligned..: "
               << target_aligned_alignment_coords.size()
               << ", alignment_coords: " << alignment_coords.size()
               << ", neighbour_aligned_alignment_coords: " << neighbour_aligned_alignment_coords.size();
            std::cout << ee.str() << std::endl;
            return 0.0;
        }

        double thet1 = 0.0;
        double thet2 = 0.5;
        morph::NM_Simplex<double> simp (thet1, thet2);
        // Set a termination threshold for the SD of the vertices of the simplex
        simp.termination_threshold = 2.0 * std::numeric_limits<double>::epsilon();
        // Set a 10000 operation limit, in case the above threshold can't be reached
        simp.too_many_operations = 1000;
        while (simp.state != morph::NM_Simplex_State::ReadyToStop) {

            if (simp.state == morph::NM_Simplex_State::NeedToComputeThenOrder) {
                // 1. apply objective to each vertex
                for (unsigned int i = 0; i <= simp.n; ++i) {
                    simp.values[i] = this->computeSos1d (target_aligned_alignment_coords,
                                                         neighbour_aligned_alignment_coords,
                                                         alignment_coords,
                                                         simp.vertices[i][0]);
                }
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToOrder) {
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeReflection) {
                double val = this->computeSos1d (target_aligned_alignment_coords,
                                                 neighbour_aligned_alignment_coords,
                                                 alignment_coords,
                                                 simp.xr[0]);
                simp.apply_reflection (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeExpansion) {
                double val = this->computeSos1d (target_aligned_alignment_coords,
                                                 neighbour_aligned_alignment_coords,
                                                 alignment_coords,
                                                 simp.xe[0]);
                simp.apply_expansion (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeContraction) {
                double val = this->computeSos1d (target_aligned_alignment_coords,
                                                 neighbour_aligned_alignment_coords,
                                                 alignment_coords,
                                                 simp.xc[0]);
                simp.apply_contraction (val);
            }
        }
        std::vector<double> vP = simp.best_vertex();
        double min_sos = simp.best_value();

        std::cout << "Best sos value: " << min_sos << " and best theta: " << vP[0] << std::endl;

        // Return the best rotation
        return vP[0];
    }

    /*!
     * By optimally aligning (by 2d translate and rotate only) the \a alignment_coords
     * with \a target_aligned_alignment_coords, find a \a translation and \a rotation
     * which are the outputs of this function.
     *
     * \param alignment_coords The coordinates which will be aligned by the
     * optimization. These could be user-supplied landmarks or the Bezier curve points.
     *
     * \param target_aligned_alignment_coords The target for alignment. alignment_coords
     * should be transformed until they match target_aligned_alignment_coords as closely
     * as possible.
     *
     * \param neighbour_aligned_alignment_coords The neighbouring slice alignment
     * coords. These are passed into the cost function of the optimization.
     *
     * \param translation Output. The x/y translation applied to alignment_coords and
     * coords to give (after \a rotation has also been applied) \a
     * target_aligned_alignment_coords
     *
     * \param rotation Output. The rotation (theta) applied to alignment_coords and
     * coords to give (as long as translation was applied)
     * target_aligned_alignment_coords.
     */
    void alignOptimally (const std::vector<cv::Point2d>& alignment_coords,
                         const std::vector<cv::Point2d>& target_aligned_alignment_coords,
                         const std::vector<cv::Point2d>& neighbour_aligned_alignment_coords,
                         cv::Point2d& translation,
                         double& rotation)
    {
        // Check that target frame has same number of coords, if not warn and return.
        if (target_aligned_alignment_coords.size() != alignment_coords.size()
            || neighbour_aligned_alignment_coords.size() != alignment_coords.size()) {
            std::stringstream ee;
            ee << "alignOptimally: Same number of alignment coords in both frames please! target_aligned..: "
               << target_aligned_alignment_coords.size()
               << ", alignment_coords: " << alignment_coords.size()
               << ", neighbour_aligned_alignment_coords: " << neighbour_aligned_alignment_coords.size();
            std::cout << ee.str() << std::endl;
            return;
        }

        // Now we have a fit which has same number of bins as the previous frame,
        // this means we can compute an objective function
        //
        // Store initial_vertices elements in order (x,y,theta).
        //
        // Use some proportion of the distance between the first landmark and its
        // equivalent on the previous slice to determine the values for the initial
        // vertices.
        cv::Point2d l1_offset = alignment_coords[0] - target_aligned_alignment_coords[0];
        double d = 0.3 * std::sqrt (l1_offset.x * l1_offset.x + l1_offset.y * l1_offset.y);
        //std::cout << "distance used for initial vertices in xytheta space: " << d << std::endl;
        morph::vVector<morph::vVector<double>> initial_vertices;
        initial_vertices.push_back ({0.0, 0.0, 0.0});
        initial_vertices.push_back ({d, 0.0, 0.0});
        initial_vertices.push_back ({d, d, 0.0});
        initial_vertices.push_back ({d, 0.0, 0.1});

        morph::NM_Simplex<double> simp (initial_vertices);
        // Set a termination threshold for the SD of the vertices of the simplex
        simp.termination_threshold = 2.0 * std::numeric_limits<double>::epsilon();
        // Set an operation limit, in case the above threshold can't be reached
        simp.too_many_operations = 1000;

        while (simp.state != morph::NM_Simplex_State::ReadyToStop) {

            if (simp.state == morph::NM_Simplex_State::NeedToComputeThenOrder) {

                // 1. apply objective to each vertex
                for (unsigned int i = 0; i <= simp.n; ++i) {
                    simp.values[i] = this->computeSos3d (target_aligned_alignment_coords,
                                                         neighbour_aligned_alignment_coords,
                                                         alignment_coords, simp.vertices[i]);
                }
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToOrder) {
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeReflection) {
                double val = this->computeSos3d (target_aligned_alignment_coords,
                                                 neighbour_aligned_alignment_coords,
                                                 alignment_coords, simp.xr);
                simp.apply_reflection (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeExpansion) {
                double val = this->computeSos3d (target_aligned_alignment_coords,
                                                 neighbour_aligned_alignment_coords,
                                                 alignment_coords, simp.xe);
                simp.apply_expansion (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeContraction) {
                double val = this->computeSos3d (target_aligned_alignment_coords,
                                                 neighbour_aligned_alignment_coords,
                                                 alignment_coords, simp.xc);
                simp.apply_contraction (val);
            }
        }
        std::vector<double> vP = simp.best_vertex();
        double min_sos = simp.best_value();

        std::cout << "After " << simp.operation_count << " operations, best sos value: " << min_sos
                  << " and best x,y,theta: (" << vP[0] << "," << vP[1] << "," << vP[2] << ")" << std::endl;

        // Set translation and rotation
        translation.x = vP[0];
        translation.y = vP[1];
        rotation = vP[2];
    }

    /*!
     * By optimally aligning (by 2d translate and rotate only) the \a alignment_coords
     * with \a target_aligned_alignment_coords, find a translation and rotation to
     * transform \a alignment_coords into \a aligned_alignment_coords and \a coords into
     * \a aligned_coords.
     *
     * \param alignment_coords The coordinates which will be aligned by the
     * optimization. These could be user-supplied landmarks or the Bezier curve points.
     *
     * \param target_aligned_alignment_coords The target for alignment. alignment_coords
     * should be transformed until they match target_aligned_alignment_coords as closely
     * as possible.
     *
     * \param translation Output. The x/y translation applied to alignment_coords and
     * coords to give (after \a rotation has also been applied) \a aligned_alignment_coords
     * and \a aligned_coords
     *
     * \param rotation Output. The rotation (theta) applied to alignment_coords and
     * coords to give (as long as translation was applied) aligned_alignment_coords.
     */
    void alignOptimally (const std::vector<cv::Point2d>& alignment_coords,
                         const std::vector<cv::Point2d>& target_aligned_alignment_coords,
                         cv::Point2d& translation,
                         double& rotation)
    {
        // Check that target frame has same number of coords, if not throw exception
        if (target_aligned_alignment_coords.size() != alignment_coords.size()) {
            std::stringstream ee;
            ee << "alignOptimally(2): Same number of alignment coords in both frames please! target_aligned..: "
               << target_aligned_alignment_coords.size() << ", alignment_coords: " << alignment_coords.size();
            std::cout << ee.str() << std::endl;
            return;
        }

        // Now we have a fit which has same number of bins as the previous frame,
        // this means we can compute SOS objective function

        // Store initial_vertices elements in order (x,y,theta). Note: Would like to
        // apply a limit on theta to avoid degenerate solutions; to avoid rotating by
        // anywhere close to 2pi, then rotation limit should be pi. So far this has not
        // been addressed.
        //
        // Use some proportion of the distance between the first landmark and its
        // equivalent on the previous slice to determine the values for the initial
        // vertices.
        cv::Point2d l1_offset = alignment_coords[0] - target_aligned_alignment_coords[0];
        double d = 0.3 * std::sqrt (l1_offset.x * l1_offset.x + l1_offset.y * l1_offset.y);
        std::cout << "distance used for initial vertices in xytheta space: " << d << std::endl;
        morph::vVector<morph::vVector<double>> initial_vertices;
        initial_vertices.push_back ({0.0, 0.0, 0.0});
        initial_vertices.push_back ({d, 0.0, 0.0});
        initial_vertices.push_back ({d, d, 0.0});
        initial_vertices.push_back ({d, 0.0, 0.1});

        morph::NM_Simplex<double> simp (initial_vertices);
        // Set a termination threshold for the SD of the vertices of the simplex
        simp.termination_threshold = 2.0 * std::numeric_limits<double>::epsilon();
        // Set an operation limit, in case the above threshold can't be reached
        simp.too_many_operations = 1000;

        while (simp.state != morph::NM_Simplex_State::ReadyToStop) {

            if (simp.state == morph::NM_Simplex_State::NeedToComputeThenOrder) {
                // 1. apply objective to each vertex
                for (unsigned int i = 0; i <= simp.n; ++i) {
                    simp.values[i] = this->computeSos3d (target_aligned_alignment_coords, alignment_coords, simp.vertices[i]);
                }
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToOrder) {
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeReflection) {
                double val = this->computeSos3d (target_aligned_alignment_coords, alignment_coords, simp.xr);
                simp.apply_reflection (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeExpansion) {
                double val = this->computeSos3d (target_aligned_alignment_coords, alignment_coords, simp.xe);
                simp.apply_expansion (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeContraction) {
                double val = this->computeSos3d (target_aligned_alignment_coords, alignment_coords, simp.xc);
                simp.apply_contraction (val);
            }
        }
        std::vector<double> vP = simp.best_vertex();
        double min_sos = simp.best_value();

        std::cout << "After " << simp.operation_count << " operations, best sos value: " << min_sos
                  << " and best x,y,theta: (" << vP[0] << "," << vP[1] << "," << vP[2] << ")" << std::endl;

        // Set translation and rotation
        translation.x = vP[0];
        translation.y = vP[1];
        rotation = vP[2];
    }
public:
    //! Common code to generate the frame name
    std::string getFrameName() const
    {
        std::stringstream ss;
        ss << "/Frame";
        ss.width(3);
        ss.fill('0');
        ss << (1+this->idx); // Count from 1 in the data file
        return ss.str();
    }

}; // FrameData
