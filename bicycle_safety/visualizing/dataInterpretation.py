class showData(object):

    def __init__(self):
        self.f = open("logansdata1.txt", "r")
        self.datapointnumber = 0;
        self.pointsofinterest = [];
        print(type(self.pointsofinterest))


    def printData(self):

        for x in self.f:
            # print("Loop number: ", self.datapointnumber)
            # print(x)
            self.datapointnumber = self.datapointnumber + 1
            temp = x.split()
            self.pointsofinterest.append([float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3])])

        self.plotData()


    def plotData(self):
        print(self.pointsofinterest)

        plt.figure()
        axes = plt.gca()
        axes.set_xlim([-5, 10])
        axes.set_ylim([-1, 7])
        largestX = (self.pointsofinterest)[0][0]
        count = 0
        curve = []
        for val,x in enumerate(self.pointsofinterest):
            plt.pause(0.01)
            # print(x[0],x[1])
            velocityvector = 1.5

            if val==0 or (self.pointsofinterest)[val][0]>largestX:
                if count<9:
                    curve.append([x[0],x[1]])
                largestX=(self.pointsofinterest)[val][0]
                length=np.sqrt(x[2]**2+x[3]**2)
                velocities=plt.plot([x[0],x[0]+x[2]*length*velocityvector],[x[1],x[1]+x[3]*length*velocityvector],'b',alpha=0.15)
                plt.plot(x[0],x[1],'ro')
                print(type(curve))
                # START DEBUGGING HERE
                if count >=9:
                    xvals = []
                    yvals = []
                    for loop in range(0,len(curve),3):
                        xvals.append(curve[loop][0])
                        yvals.append(curve[loop][1])
                    xvals = np.array(xvals)
                    yvals = np.array(yvals)
                    zvals = np.polyfit(xvals, yvals, 3)
                    xp = np.linspace(curve[0][0], curve[-1][0], 100)
                    p30 = np.poly1d(np.polyfit(xvals, yvals, 30))
                    p = np.poly1d(zvals)
                    plt.plot(xvals, yvals, '.', xp, p(xp), '-', xp, p30(xp), '--')
                    curve.pop(0)
                    curve.append([x[0],x[1]])
                count = count + 1
                velocities.remove
                plt.show


            else:
                plt.plot(x[0], x[1], 'ko', alpha=0.1)
                plt.show


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import time

    logansdata=showData()
    logansdata.printData()