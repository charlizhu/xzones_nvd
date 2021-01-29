# Step 1 is to use opencv just to capture frames.
# Refer to this source: http://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_gui/py_video_display/py_video_display.html

import numpy as np
import cv2

cap = cv2.VideoCapture(0)

print(cap.isOpened())

red = np.array([247, 230, 146])
yellow = np.array([107, 112, 2])

while(True):
    # Capture frame-by-frame
    ret, frame = cap.read()

    # Our operations on the frame come here
    colored = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)

    # See this stack exchange: https://stackoverflow.com/questions/56474278/how-to-mirror-live-webcam-video-when-using-cv2
    colored_flip = cv2.flip(colored, 1)
    # combined_window = np.hstack([gray, gray_flip])

    # Threshold the HSV image to get only blue colors
    mask = cv2.inRange(colored_flip, yellow, red)

    # Bitwise-AND mask and original image
    res = cv2.bitwise_and(colored_flip,colored_flip, mask= mask)

    # Refer to this source for pixel color detection: https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_colorspaces/py_colorspaces.html

    cv2.imshow('frame',colored_flip)
    # cv2.imshow('mask',mask)
    cv2.imshow('res',res)
    # Display the resulting frame
    # x = cv2.imshow('frame',colored_flip)
    # print(colored_flip)
    if cv2.waitKey(1) & 0xFF == ord('f'):
        break

# When everything done, release the capture
cap.release()
cv2.destroyAllWindows()