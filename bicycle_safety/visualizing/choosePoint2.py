# NOTE: This script is currently incomplete, and seems to work so far, though I am more than certain that I missed
# many edge cases. With the current data it seems to work fine, but that might not always be the case.

class choosePoint(object):

    def __init__(self):
        self.f = open("bytes2.txt", "r")
        self.tid = []
        self.tidcounts = []
        self.currenttids = []
        self.lx = []
        self.ly = []
        self.sine = math.sin(math.pi/6)
        self.cosine = math.cos(math.pi/6)
        self.detectItems()

    def detectItems(self): # In reality this would be reading from the UART.
        for x in self.f:
            datapoint = x.split()
            self.trackItems(datapoint)
            print(self.tid)
            print(self.tidcounts)
            # self.loganskalmanfilterscript() # TO DO --> use the global variables as inputs.
        # Plotting is just for prototyping visualization.
        # print(self.lx)
        # print(self.ly)
        plt.plot(self.lx,self.ly)
        plt.axis('equal')
        plt.show()



    def trackItems(self,datapoint):
        for carIndex in range(0,len(datapoint),14):
            # This next snippet is taken from the TI MATLAB GUI script
            xLoc = float(datapoint[carIndex + 2]) + float(datapoint[carIndex + 3]) * 256
            yLoc = float(datapoint[carIndex + 4]) + float(datapoint[carIndex + 5]) * 256
            if xLoc > 32767:
                xLoc = xLoc - 65536
            if yLoc > 32767:
                yLoc = yLoc - 65536
            xp = xLoc / 128
            yp = yLoc / 128
            xLoc = self.cosine * xp - self.sine * yp
            yLoc = self.sine * xp + self.cosine * yp
            # ^^ from TI's script.
            if yp < 5: # negating any points that are too far.
                if datapoint[carIndex] not in self.tid:
                    (self.tid).append(datapoint[carIndex])
                    (self.tidcounts).append(1)
                    # self.lx = []
                    # self.ly = []
                else:
                    indexnumber = (self.tid).index(datapoint[carIndex])
                    (self.tidcounts)[indexnumber] = (self.tidcounts)[indexnumber] + 1
                (self.lx).append(xLoc)
                (self.ly).append(yLoc)


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    choosePoint()