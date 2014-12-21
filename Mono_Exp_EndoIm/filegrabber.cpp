/* This file is part of the SceneLib2 Project.
* http://hanmekim.blogspot.com/2012/10/scenelib2-monoslam-open-source-library.html
* https://github.com/hanmekim/SceneLib2
*
* Copyright (c) 2012 Hanme Kim (hanme.kim@gmail.com)
*
* SceneLib2 is an open-source C++ library for SLAM originally designed and
* implemented by Andrew Davison and colleagues at the University of Oxford.
*
* I reimplemented his version with the following objectives;
*  1. Understand his MonoSLAM algorithm in code level.
*  2. Replace older libraries (i.e. VW34, GLOW, VNL, Pthread) with newer ones
*     (Pangolin, Eigen3, Boost).
*  3. Support USB camera instead of IEEE1394.
*  4. Make it more portable and convenient by using CMake and git repository.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy of
* this software and associated documentation files (the "Software"), to deal in
* the Software without restriction, including without limitation the rights to
* use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
* of the Software, and to permit persons to whom the Software is furnished to do
* so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#include "filegrabber.h"
//using namespace cv;
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>


FileGrabber::FileGrabber():init_frame_id_(0), max_frame_(0)
{
	initialised_ = false;
	frame_id_ = 0;
	frame_grabber_ = NULL;
}

FileGrabber::~FileGrabber()
{
	if(!files_vec_.empty()) {
		files_vec_.clear();
	}
	if(frame_grabber_ != NULL)
		delete frame_grabber_;

}

void FileGrabber::SetFrameRange(int init_frame_id, int max_frame)
{
	init_frame_id_ = init_frame_id;
	frame_id_ = init_frame_id_ - 1;
	max_frame_ = max_frame;
}


void FileGrabber::Init(const string &path, FrameGrabber *frame_grabber, int init_frame_id, int max_frame)
{
	ProcessFiles(path);
	frame_grabber_  = frame_grabber;
	initialised_    = true;
	SetFrameRange(init_frame_id, max_frame);
	fg_thread_ = boost::thread(boost::ref(*this));
	//fg_thread_.swap(boost::thread(boost::ref(*this)));
	//cout << "here is the first thread" << endl;
}

void FileGrabber::ProcessFiles(const boost::filesystem::path &directory)
{
	if (exists(directory)) {
		boost::filesystem::directory_iterator end;

		for (boost::filesystem::directory_iterator iter(directory); iter!=end; ++iter) {
			if (is_directory(*iter)) {
				ProcessFiles(*iter);
			} else {
				string  path_name = iter->path().string();
				/*cout << "path name is " << path_name <<endl;*/
				files_vec_.push_back(path_name);
			}
		}
		sort(files_vec_.begin(), files_vec_.end());
	}
}

void FileGrabber::operator ()()
{

	while (initialised_) {
		if (frame_grabber_->IsFrameBufferFull()==false) {
			if (files_vec_.size() > (unsigned int)frame_id_ && max_frame_ > frame_id_) {
				Frame_ frame;
				string  file_full_path = files_vec_.at(frame_id_);
				frame.frame_id = frame_id_;
				frame.data = GetImageFile(file_full_path);
				++frame_id_;
				frame_grabber_->SetFrame(frame);
				
			}
		}
		else {
			boost::this_thread::sleep(boost::posix_time::milliseconds(1));
		}
	}
}

cv::Mat FileGrabber::GetImageFile(const string &file_full_path)
{
	return  cv::imread(file_full_path, 0); // reads the image without color
}


