#include "framegrabber.h"

class UsbCamGrabber {
public:
	UsbCamGrabber();
	~UsbCamGrabber();
	void Init(const string &path, FrameGrabber *frame_grabber);


};
