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
        }
        return DM::pInstance;
    }
    //! Add a frame to vFrameData
    void addFrame (Mat& frameImg) {
        this->vFrameData.push_back (FrameData (frameImg.clone()));
    }
    //! Return the size of vFrameData
    unsigned int getNumFrames (void) const {
        return this->vFrameData.size();
    }
    //! get current frame. Short name on purpose.
    FrameData* gcf (void) {
        if (!this->vFrameData.empty()) {
            return &(this->vFrameData[I]);
        }
        return (FrameData*)0;
    }
    //! Get the current frame number, counting from 1 like a human.
    int getFrameNum (void) const {
        return 1+I;
    }
    // Get a pointer to the persistent Mat img member attribute
    Mat* getImg (void) {
        return &(img);
    }
    //! Make the next frame current (or cycle back to the first)
    void nextFrame (void) {
        ++I %= this->vFrameData.size();
    }
    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->vFrameData[I].frame.clone();
    }
};

//! Globally initialise DM instance pointer to NULL
DM* DM::pInstance = 0;
