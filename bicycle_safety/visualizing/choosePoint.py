class choosePoint(object):

    def __init__(self):
        self.f = open("bytes1.txt", "r")
        self.tid = []
        self.tidcounts = []
        self.currenttids = []
        self.lx = []
        self.ly = []
        self.detectItems()

    def detectItems(self): # In reality this would be reading from the UART.
        for x in self.f:
            datapoint = x.split()
            for y in range(0,len(datapoint),14):
                # print(datapoint[y])
                (self.currenttids).append(datapoint[y])
            self.trackItems()
            self.currenttids = []
            # self.isCar1()
            car = self.isCar2()
            self.isValid(car,datapoint)
            # self.loganskalmanfilterscript() # TO DO
        # Plotting is just for prototyping visualization.
        print(self.lx)
        print(self.ly)
        plt.plot(self.lx,self.ly)
        plt.axis('equal')
        plt.show()

    def isValid(self,car,datapoint):
        for carIndex in range(0, len(datapoint), 14):
            if datapoint[carIndex] == car:
                print(car)
                # This next snippet is taken from the TI MATLAB GUI script
                xLoc = float(datapoint[carIndex+2]) + float(datapoint[carIndex+3])*256
                yLoc = float(datapoint[carIndex+4]) + float(datapoint[carIndex+5])*256
                if xLoc > 32767:
                    xLoc = xLoc - 65536
                if yLoc > 32767:
                    yLoc = yLoc - 65536
                xLoc = xLoc/128
                yLoc = yLoc/128
                print(xLoc,yLoc)
                if yLoc < 4.5: # Trying to exclude all points that are too far away, which are assumed to be erroneous datapoints.
                    (self.lx).append(xLoc)
                    (self.ly).append(yLoc)
                # print(self.lx)
                # print(self.ly)

    def trackItems(self):
        for x in self.currenttids:
            if x not in self.tid:
                (self.tid).append(x)
                (self.tidcounts).append(1)
            else:
                indexnumber = (self.tid).index(x)
                (self.tidcounts)[indexnumber] = (self.tidcounts)[indexnumber] + 1
        # print(self.tidcounts)
        # print(self.tid)

    # Assuming that the AWR correctly tracks cars if it is numbered more than n times.
    def isCar1(self):
        for i,j in enumerate(self.tidcounts):
            if i > 50: # 50 is arbitrary for now
                print((self.tid)[j])

    # Assuming that the AWR correctly tracks the car if it is TID'ed the most times.
    def isCar2(self):
        if max(self.tidcounts)>10: # 10 is arbitrary for now
            # print(type(self.tid))
            # print(self.tidcounts)
            return((self.tid)[(self.tidcounts).index(max(self.tidcounts))])


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    choosePoint()