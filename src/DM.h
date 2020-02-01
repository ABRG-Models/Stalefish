#pragma once

#include <vector>
using std::vector;
#include <opencv2/opencv.hpp>
using cv::Mat;
#include "FrameData.h"

//! Singleton pattern data manager class to hold framedata
class DM
{
private:
    //! Private constructor/destructor
    DM() {};
    ~DM() {};
    //! A pointer returned to the single instance of this class
    static DM* pInstance;
    //! The frame data which is being managed
    vector<FrameData> vFrameData;
    //! The index into frames which is the current image
    int I = 0;
    //! The current image
    Mat img;

public:
    //! The instance public function. Short on purpose
    static DM* i (void) {
        if (DM::pInstance == 0) {
            DM::pInstance = new DM;
            DM::i()->init();
        }
        return DM::pInstance;
    }
    void init (void) {
        this->vFrameData.clear();
    }
    //! Add a frame to vFrameData
    void addFrame (Mat& frameImg) {
        cout << "addFrame called" << endl;
        FrameData fd(frameImg);
        cout << "FrameData made" << endl;
        this->vFrameData.push_back (fd);
        cout << "FrameData pushed back" << endl;
    }
    //! Return the size of vFrameData
    unsigned int getNumFrames (void) const {
        return this->vFrameData.size();
    }
    //! get current frame. Short name on purpose.
    FrameData* gcf (void) {
        if (!this->vFrameData.empty()) {
            return &(this->vFrameData[this->I]);
        }
        return (FrameData*)0;
    }
    //! Get the current frame number, counting from 1 like a human.
    int getFrameNum (void) const {
        return 1+this->I;
    }
    // Get a pointer to the persistent Mat img member attribute
    Mat* getImg (void) {
        return &(this->img);
    }
    //! Make the next frame current (or cycle back to the first)
    void nextFrame (void) {
        ++this->I %= this->vFrameData.size();
        DM::i()->gcf()->binA = DM::i()->binA;
        DM::i()->gcf()->binB = DM::i()->binB;
        DM::i()->gcf()->nBinsTarg = DM::i()->nBinsTarg;
    }
    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->vFrameData[this->I].frame.clone();
    }
    //! The application window name
    const string winName = "StaleFish";
    //! Saved/last cursor position
    int x = 0;
    int y = 0;
    //! Target number of bins; used by bins slider. Apply this to the framedata
    int nBinsTarg = 100;
    //! The bin lengths, set with a slider.
    int binA = 0;
    int binB = 40;
};

//! Globally initialise DM instance pointer to NULL
DM* DM::pInstance = 0;
