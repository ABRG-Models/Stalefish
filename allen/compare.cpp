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
using std::cout;
using std::endl;
#include <set>
using std::set;


int main()
{
    string fn = "e100041580_99_100873513";
    string fn_ex = fn + "_expr.jpg";
    string fn_ish = fn + ".jpg";

    Mat fr_ex = imread (fn_ex.c_str(), IMREAD_COLOR);
    Mat fr_ish = imread (fn_ish.c_str(), IMREAD_COLOR);

    if (fr_ex.rows != fr_ish.rows
        || fr_ex.cols != fr_ish.cols) {
        cout << "Size mismatch" << endl;
        return -1;
    }

#if 0
    cout << "Image size is " << (fr_ex.rows * fr_ex.cols) << " px" << endl;

    cout << "expression frame has " << fr_ex.channels() << " channels with depth " << fr_ex.depth() << endl;
    cout << "ISH frame has " << fr_ish.channels() << " channels with depth " << fr_ish.depth() << endl;
    if (fr_ish.depth() == CV_8U) {
        cout << "depth 0 is CV_8U\n";
    } else {
        cout << "depth 0 is NOT CV_8U\n";
    }
#endif

#if 0
    cout << "a = [\n";
    for (int r = 0; r < fr_ex.rows; ++r) {
        for (int c = 0; c < fr_ex.cols; ++c) {
            // Example pixel (r,c)
            Point2i p(r,c);
            Vec3b pix_ex = fr_ex.at<Vec3b>(p);
            // FIXME: Convert pix_ex from colour map to scalar.
            if (pix_ex[0] || pix_ex[1] || pix_ex[2]) {
                //cout << "Pixel Expr (" << r << "," << c << ") = " << pix_ex << endl;
                //cout << "Pixel ISH  (" << r << "," << c << ") = " << (fr_ish.at<Vec3b>(p)) << endl;
                Vec3b ish = fr_ish.at<Vec3b>(p);
                cout << (int)ish[0] << "," << (int)ish[1] << "," << (int)ish[2] << endl;
            }
        }
    }
    cout << "];\n";

    cout << "n = [\n";
    for (int r = 0; r < fr_ex.rows; ++r) {
        for (int c = 0; c < fr_ex.cols; ++c) {
            // Example pixel (r,c)
            Point2i p(r,c);
            Vec3b pix_ex = fr_ex.at<Vec3b>(p);
            if (pix_ex[0] || pix_ex[1] || pix_ex[2]) {
                // Do nothing this is an expressing pixel
            } else {
                // Non-expressing
                Vec3b ish = fr_ish.at<Vec3b>(p);
                int sum = ((int)ish[0] + (int)ish[1] + (int)ish[2]);
                if (sum && (sum < 690)) {
                    // Then it might not be white background...
                    cout << (int)ish[0] << "," << (int)ish[1] << "," << (int)ish[2] << endl;
                }
            }
        }
    }
    cout << "];\n";
#endif

#if 1
    namedWindow ("ISH", WINDOW_AUTOSIZE );// Create a window for display.
    imshow ("ISH", fr_ish);
    // Wait for a key, then exit
    waitKey();
#endif
    return 0;
}
