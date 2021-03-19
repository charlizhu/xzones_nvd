# NOTE: This script is currently incomplete, and seems to work so far, though I am more than certain that I missed
# many edge cases. With the current data it seems to work fine, but that might not always be the case.

# Current problems I can think of:
# 1. It currently works assuming that there is only 1 car in the vicinity.
# 2. Because of the above, it does not actually make a decision as to which TID is possibly the car.
# 3. I think that the script _should_ assume that the car being TID'ed the most is the correct car. It currently does not do that.
# 4. Eventually it should tell which car is most critical if two cars might be present, for example. Though it can be justified that this is not necessary because you can't fit two cars within the vicinity of interest.
# 5. The biggest problem is that the current script cannot switch from detecting one car to the next. i.e.: one car passes by, and another approaches.

class choosePoint(object):

    def __init__(self):
        self.start = timeit.default_timer()
        self.filename = "Converted" + "\\" + "moving_rearhook_normal_speed_3.txt"
        # self.writefile = open("parsed_" + self.filename, 'a')
        self.f = open(self.filename, "r")
        self.tidcurrent = []
        self.tidprevious = []
        self.tidtracked = []
        self.tidcounts = 0
        self.currenttid = None
        self.currenttidrepeats = None
        self.lx = []
        self.ly = []
        self.lxwant = []
        self.lywant = []
        self.vx = []
        self.vy = []
        self.xpoints = []
        self.ypoints = []
        self.angle = 2*math.pi/6
        self.sine = math.sin(self.angle)
        self.cosine = math.cos(self.angle)
        self.detectItems()

    def detectItems(self): # In reality this would be reading from the UART.

        for cars in self.f:
            datapoint = cars.split()
            for carIndex in range(0, len(datapoint), 14):
                print(type(datapoint[carIndex]))
                # This next snippet is taken from the TI MATLAB GUI script
                xLoc = float(datapoint[carIndex + 2]) + float(datapoint[carIndex + 3]) * 256
                yLoc = float(datapoint[carIndex + 4]) + float(datapoint[carIndex + 5]) * 256
                if xLoc > 32767:
                    xLoc = xLoc - 65536
                if yLoc > 32767:
                    yLoc = yLoc - 65536
                xp = xLoc / 128
                yp = yLoc / 128
                xLocRot = self.cosine * xp - self.sine * yp
                yLocRot = self.sine * xp + self.cosine * yp
                xVel = float(datapoint[carIndex + 6]) + float(datapoint[carIndex + 7]) * 256
                yVel = float(datapoint[carIndex + 8]) + float(datapoint[carIndex + 9]) * 256
                if xVel > 32767:
                    xVel = xVel - 65536
                if yVel > 32767:
                    yVel = yVel - 65536
                xv = xVel / 128
                yv = yVel / 128
                xVelRot = self.cosine * xv - self.sine * yv
                yVelRot = self.sine * xv + self.cosine * yv
                if datapoint[carIndex] == 152:
                    (self.lxwant).append(xLocRot)
                    (self.lywant).append(yLocRot)
                else:
                    (self.lx).append(xLocRot)
                    (self.ly).append(yLocRot)
                (self.vx).append(xVelRot)
                (self.vy).append(yVelRot)
                # print(self.lx,self.ly)
        plt.plot(self.lx,self.ly,'k.',alpha=0.2)
        plt.plot(self.lxwant, self.lywant, 'r')
        # plt.xlim([-5,20])
        plt.grid('on')
        plt.show()

        #     self.trackItems(datapoint)
        #     print(self.tid)
        #     print("hi", self.tidcounts)
        #     if self.tidcounts:
        #         print("hey", max(self.tidcounts))
        #         if max(self.tidcounts)>5:
        #             print("Logan's script here")
        #     # self.loganskalmanfilterscript() # TO DO --> use the global variables as inputs.
        # # Plotting is just for prototyping visualization.
        # # print(self.lx)
        # # print(self.ly)
        # plt.plot(self.lx,self.ly)
        # plt.axis('equal')
        # plt.show()


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import tracemalloc
    import timeit
    tracemalloc.start()
    choosePoint()