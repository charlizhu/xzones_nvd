class choosePoint(object):

    def __init__(self):
        self.start = timeit.default_timer()
        self.filename = "Converted" + "\\" + "moving_rearhook_normal_speed_3.txt"
        self.writefile = open("screened_data\items.txt", 'a')
        self.f = open(self.filename, "r")
        self.numskips = 2
        self.commonitems = []
        self.elements = []
        self.tidcurrent = []
        self.tidprevious = []
        self.lx = []
        self.ly = []
        self.vx = []
        self.vy = []
        self.xpoints = []
        self.ypoints = []
        self.xvelocities = []
        self.yvelocities = []
        self.angle = -math.pi/6
        self.sine = math.sin(self.angle)
        self.cosine = math.cos(self.angle)
        self.convertFormat()

    def convertFormat(self): # In reality this would be reading from the UART.

        itemnumber = 1
        for cars in self.f:
            first_size, first_peak = tracemalloc.get_traced_memory()
            stop = timeit.default_timer()
            print(stop - self.start)
            self.start = stop
            print(first_size / 10 ** 6, first_peak / 10 ** 6)
            datapoint = cars.split()
            for carIndex in range(0, len(datapoint), 14):

                if itemnumber <= self.numskips:
                    (self.tidprevious).append([tid for tid in range(0, len(datapoint), 14)])
                    continue

                # Parts of this next snippet is taken from the TI MATLAB GUI script
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
                # End of MATLAB stuff
                (self.tidcurrent).append(datapoint[carIndex])
                (self.lx).append(xLocRot)
                (self.ly).append(yLocRot)
                (self.vx).append(xVelRot)
                (self.vy).append(yVelRot)
            itemnumber = itemnumber + 1

            items = (self.dataScreening())
            if (items == None):
                print("NONE")
            elif (not items):
                print("EMPTY")
            else:
                for x in items:
                    print((self.tidcurrent)[x])
                    (self.xpoints).append((self.lx)[x])
                    (self.ypoints).append((self.ly)[x])
                    (self.xvelocities).append((self.vx)[x])
                    (self.yvelocities).append((self.vy)[x])

            # Resetting and updating
            (self.tidprevious).pop(0)
            (self.tidprevious).append(self.tidcurrent)
            self.tidcurrent = []
            self.lx = []
            self.ly = []
            self.vx = []
            self.vy = []
        plt.plot(self.xpoints,self.ypoints,'r')
        plt.grid("on")

        for x in range(0,len(self.xpoints)):
            # print(type(' '))
            tofile = str(self.xpoints[x]) + " " + str(self.ypoints[x]) + " " + str(self.xvelocities[x]) + " " + str(self.yvelocities[x]) + "\n"
            # print(tofile)
            (self.writefile).write(tofile)
        plt.show()

    def dataScreening(self):

        trackedcars = self.trackItems()

        if not trackedcars:
            print("Nothing tracked")
        else:
            orderlist = []
            for item in trackedcars:
                order = (self.tidcurrent).index(item)
                if (((self.ly)[order] < 5) and ((self.lx)[order] > -10)) and ((self.vx)[order] > 0):
                # if (((self.ly)[order] < 5) and ((self.lx)[order] > 0)) and ((self.vx)[order] > 0): IF FRONT HOOK
                    orderlist.append(order)
            return orderlist

    def trackItems(self): # needs to be fixed
        if not self.tidprevious:
            print("empty")
            pass
        else:
            self.commonitems = (self.tidprevious)[0]
            for element in (self.tidprevious)[1:]:
                self.commonitems = [item for item in self.commonitems if item in element]
            return [item for item in self.commonitems if item in self.tidcurrent]


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import tracemalloc
    import timeit
    tracemalloc.start()
    choosePoint()