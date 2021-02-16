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
        self.f = open("bytes2.txt", "r")
        self.tidcurrent = []
        self.tidprevious = []
        self.tidtracked = []
        self.tidcounts = 0
        self.currenttid = None
        self.currenttidrepeats = None
        self.lx = []
        self.ly = []
        self.vx = []
        self.vy = []
        self.xpoints = []
        self.ypoints = []
        self.sine = math.sin(math.pi/6)
        self.cosine = math.cos(math.pi/6)
        self.detectItems()

    def detectItems(self): # In reality this would be reading from the UART.
        for cars in self.f:
            datapoint = cars.split()
            for carIndex in range(0, len(datapoint), 14):
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
                if yLocRot > 5 or xLocRot < -2: # to be modified, need more conditions here...
                    # print(xLocRot,yLocRot)
                    continue
                # ^^ from TI's script.
                (self.tidcurrent).append(datapoint[carIndex])
                (self.lx).append(xLocRot)
                (self.ly).append(yLocRot)
                (self.vx).append(xVelRot)
                (self.vy).append(yVelRot)
            if not self.tidcurrent:
                pass
            self.trackItems()
            if (self.currenttidrepeats is not None) and (self.currenttid is not None):
                # print(self.currenttid, (self.tidcurrent).index(self.currenttid))
                self.Kalman()
            self.lx = []
            self.ly = []
            self.vx = []
            self.vy = []
            self.tidcurrent = []
        plt.plot(self.xpoints,self.ypoints,'r-')
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

    def trackItems(self):
        if not self.tidprevious:
            self.tidprevious = self.tidcurrent
            pass
        else:
            for thisitem in self.tidcurrent:
                if thisitem in self.tidprevious:
                    if thisitem not in self.tidtracked:
                        (self.tidtracked).append(thisitem)
            for thatitem in self.tidtracked:
                if thatitem not in self.tidcurrent:
                    (self.tidtracked).pop(self.tidtracked.index(thatitem))
            if len(set(self.tidcurrent).intersection(set(self.tidprevious))) > 0:
                self.tidcounts = self.tidcounts + 1
            else:
                self.tidcounts = 0
            self.tidprevious = self.tidcurrent
        # print(self.tidtracked, self.tidcounts)

        if len(self.tidtracked) >= 1:
            self.currenttid = self.tidtracked[0]
            self.currenttidrepeats = self.tidcounts
        else:
            self.currenttid = None
            self.currenttidrepeats = None
            # # I'm sure there is a logic error below...
            # if yp < 5: # negating any points that are too far.
            #     if datapoint[carIndex] not in self.tid:
            #         (self.tid).append(datapoint[carIndex])
            #         (self.tidcounts).append(1)
            #
            #         # self.lx = []
            #         # self.ly = []
            #     else:
            #         indexnumber = (self.tid).index(datapoint[carIndex])
            #         (self.tidcounts)[indexnumber] = (self.tidcounts)[indexnumber] + 1
            #
            #     (self.lx).append(xLoc)
            #     (self.ly).append(yLoc)
    def Kalman(self):
        print("Logan's Kalman filter script")
        print("Position ",(self.lx)[(self.tidcurrent).index(self.currenttid)],(self.ly)[(self.tidcurrent).index(self.currenttid)])
        print("Velocity ",(self.vx)[(self.tidcurrent).index(self.currenttid)],(self.vy)[(self.tidcurrent).index(self.currenttid)])
        (self.xpoints).append((self.lx)[(self.tidcurrent).index(self.currenttid)])
        (self.ypoints).append((self.ly)[(self.tidcurrent).index(self.currenttid)])



if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    choosePoint()