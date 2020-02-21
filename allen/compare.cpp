/*!
 * Load the Allen ISH and corresponding "expression" images. For all pixels that Allen
 * reckons are "expressing", collect the RGB value in the ISH image.
 */

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using namespace cv;
#include <string>
using std::string;
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <set>
using std::set;
#include <cmath>
using std::abs;

// To allow comparison of Vec<uchar,3>s and thus creation of set<Vec<uchar,3>>
// This is a functor; a class that can be called like a function.
struct VecCompare
{
    bool operator() (const Vec<uchar,3>& lhs, const Vec<uchar,3>& rhs) const {
        if (lhs[0] < rhs[0]) {
            return true;
        } else if (lhs[1] < rhs[1]) {
            return true;
        } else if (lhs[2] < rhs[2]) {
            return true;
        }
        return false;
    }
};

int main()
{
    string fn = "e100041580_99_100873513";
    string fn_ex = fn + "_expr.jpg";
    string fn_ish = fn + ".jpg";

    Mat fr_ish = imread (fn_ish.c_str(), IMREAD_ANYCOLOR);
    Mat fr_ex = imread (fn_ex.c_str(), IMREAD_ANYCOLOR);

    if (fr_ex.rows != fr_ish.rows || fr_ex.cols != fr_ish.cols) {
        cout << "Size mismatch" << endl;
        return -1;
    }

#if 0
    cout << "Image size is " << fr_ex.rows << "x" << fr_ex.cols << " = " << (fr_ex.rows * fr_ex.cols) << " px" << endl; // 402447
    cout << "expression frame has " << fr_ex.channels() << " channels with depth " << fr_ex.depth() << endl;
    cout << "ISH frame has " << fr_ish.channels() << " channels with depth " << fr_ish.depth() << endl;
    if (fr_ish.depth() == CV_8U) {
        cout << "depth 0 is CV_8U\n";
    } else {
        cout << "depth 0 is NOT CV_8U\n";
    }
#endif

#if 0
    unsigned int pcount = 0;
    for (int r = 0; r < fr_ish.rows; ++r) {
        for (int c = 0; c < fr_ish.cols; ++c) {
            Point2i p(c,r);
            Vec<uchar,3> ish = fr_ish.at<Vec<uchar,3>>(p);
            cout << "ISH: " << static_cast<unsigned int>(ish[0])
                 << "," << static_cast<unsigned int>(ish[1])
                 << "," << static_cast<unsigned int>(ish[2]) << endl;
            pcount++;
        }
    }
    cout << "pcount = " << pcount << endl;
#endif

#if 1
    unsigned long long int pcount = 0;
    unsigned long long int ne_count = 0; // non-expressing count
    unsigned long long int se_count = 0; // slightly-expressing count
    unsigned long long int e_count = 0; // expressing count
    set<Vec<uchar,3>, VecCompare> expressing_ish;
    set<Vec<uchar,3>, VecCompare> nonexpressing_ish;
    set<Vec<uchar,3>, VecCompare> slightexpressing_ish;
    for (int r = 0; r < fr_ex.rows; ++r) {
        for (int c = 0; c < fr_ex.cols; ++c) {
            // Example pixel (c,r)
            Point2i p(c,r);
            Vec<uchar,3> pix_ex = fr_ex.at<Vec<uchar,3>>(p);
            unsigned int psum = static_cast<unsigned int>(pix_ex[0])
                + static_cast<unsigned int>(pix_ex[1])
                + static_cast<unsigned int>(pix_ex[2]);
            // FIXME: Convert pix_ex from colour map to scalar.
            Vec<uchar,3> ish = fr_ish.at<Vec<uchar,3>>(p);
            if (psum > 30) {
                e_count++;
                expressing_ish.insert (ish);
            } else if (psum == 0) {
                ne_count++;
                nonexpressing_ish.insert (ish);
            } else {
                se_count++;
                slightexpressing_ish.insert (ish);
            }
            ++pcount;
        }
    }
    cout << "pcount = " << pcount << ";\n";
    cout << "e_count = " << e_count << ";\n";
    cout << "ne_count = " << ne_count << ";\n";
    cout << "se_count = " << se_count << ";\n";

    // Expression/strongly expressing
    cout << "a = [\n";
    for (auto ish : expressing_ish) {
        cout << static_cast<unsigned int>(ish[0])
             << "," << static_cast<unsigned int>(ish[1])
             << "," << static_cast<unsigned int>(ish[2]) << endl;
    }
    cout << "];\n";

    // Non
    cout << "n = [\n";
    for (auto ishn : nonexpressing_ish) {
        cout << static_cast<unsigned int>(ishn[0])
             << "," << static_cast<unsigned int>(ishn[1])
             << "," << static_cast<unsigned int>(ishn[2]) << endl;
    }
    cout << "];\n";

    // Slight
    cout << "s = [\n";
    for (auto ishn : slightexpressing_ish) {
        cout << static_cast<unsigned int>(ishn[0])
             << "," << static_cast<unsigned int>(ishn[1])
             << "," << static_cast<unsigned int>(ishn[2]) << endl;
    }
    cout << "];\n";
#endif

#if 0
    namedWindow ("ISH", WINDOW_AUTOSIZE );// Create a window for display.
    imshow ("ISH", fr_ish);
    // Wait for a key, then exit
    waitKey();
#endif
    return 0;
}
