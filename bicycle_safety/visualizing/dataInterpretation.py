class showData(object):

    def __init__(self):
        self.f = open("logansdata1.txt", "r")
        self.datapointnumber = 0;
        self.pointsofinterest = [];
        # print(type(self.pointsofinterest))


    def printData(self):

        for x in self.f:
            # print("Loop number: ", self.datapointnumber)
            # print(x)
            self.datapointnumber = self.datapointnumber + 1
            temp = x.split()
            self.pointsofinterest.append([float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3])])

        self.plotData()


    def plotData(self):
        # print(self.pointsofinterest)

        plt.figure()
        plt.pause(10)
        axes = plt.gca()
        axes.set_xlim([-5, 10])
        axes.set_ylim([-1, 14])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.grid()
        largestX = (self.pointsofinterest)[0][0]
        count = 0
        curve = []
        for val,x in enumerate(self.pointsofinterest):
            plt.pause(0.1)
            plt.cla()
            axes.set_xlim([-5, 10])
            axes.set_ylim([-1, 14])
            plt.gca().set_aspect('equal', adjustable='box')
            plt.grid()
            # print(x[0],x[1])
            velocityvector = 1.5

            if val==0 or (self.pointsofinterest)[val][0]>largestX:
                if count<9:
                    curve.append([x[0],x[1]])
                largestX=(self.pointsofinterest)[val][0]
                speed=np.sqrt(x[2]**2+x[3]**2)
                velocities=plt.plot([x[0],x[0]+x[2]*speed*velocityvector],[x[1],x[1]+x[3]*speed*velocityvector],'b',alpha=0.5)
                plt.plot(x[0],x[1],'ro')
                # print(curve)
                # START DEBUGGING HERE
                if count >=9:
                    xvals = []
                    yvals = []
                    for loop in range(0,len(curve)):
                        xvals.append(curve[loop][0])
                        yvals.append(curve[loop][1])
                    xvals = np.array(xvals)
                    yvals = np.array(yvals)
                    # zvals = np.polyfit(xvals, yvals, 3)
                    xp = np.linspace(curve[0][0], curve[-1][0], num=100,endpoint=True)
                    f2 = interp1d(xvals, yvals, kind='cubic')
                    plt.plot(xp, f2(xp), '--')
                    # print([[xp[-1],f2(xp)[-1]],[xp[-2],f2(xp)[-2]]]) # seems very irregular still...
                    # p30 = np.poly1d(np.polyfit(xvals, yvals, 30))
                    # p = np.poly1d(zvals)
                    # plt.plot(xvals, yvals, '.', xp, p(xp), '-', xp, p30(xp), '--')
                    curve.pop(0)
                    curve.append([x[0],x[1]])
                    # Fitting the spline
                    model = LinearRegression().fit(xp.reshape((-1,1)),f2(xp))
                    # print("Linear Regression results: ", model.intercept_,model.coef_)
                    plt.plot([x[0],x[0]+2.5],[model.coef_[0]*x[0]+model.intercept_,model.coef_[0]*(x[0]+2.5)+model.intercept_],'m',alpha=0.5)
                    front_pred = model.predict(np.linspace(0,10,20).reshape(-1,1))
                    rear_pred = model.predict(np.linspace(-10, 0, 20).reshape(-1, 1))
                    if x[0]<=0 and speed>0 and any(i <= 1 for i in rear_pred):
                        plt.title("REAR DANGER at loop #" + str(val))
                    elif x[0]>0 and speed>0 and any(i <= 1 for i in front_pred):
                        plt.title("FRONT DANGER at loop #" + str(val))
                    else:
                        plt.title("Safe")
                count = count + 1
                velocities.remove
                plt.plot(0,0,"k>",markersize=25)
                plt.ylabel('Distance to your left (meters)')
                plt.xlabel('Direction of travel (meters)')
                plt.show


            else:
                plt.plot(x[0], x[1], 'ko', alpha=0.1)
                plt.plot(0, 0, "k>",markersize=25)
                plt.ylabel('Distance to your left (meters)')
                plt.xlabel('Direction of travel (meters)')
                plt.title("Noisy or faulty data")
                plt.show




if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import time
    from scipy.interpolate import interp1d
    from sklearn.linear_model import LinearRegression
    import matplotlib.animation as animation

    logansdata=showData()
    logansdata.printData()