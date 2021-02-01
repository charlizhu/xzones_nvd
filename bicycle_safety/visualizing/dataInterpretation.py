# Creator: Charlie Zhu, Jan 2021
# Purpose: Virtual prototype to show concept of classification of the right-hook collision types.
# Notes: Still needs major re-work. Please rename variables as you see fit.

class showData(object):

    def __init__(self):
        # Taking file of data generated from the IWR1843 MATLAB GUI, which has been re-parsed in a text file.
        # Format of the text file: X pos, Y pos, X vel, Y vel, X acc, Y acc
        self.f = open("logansdata1.txt", "r")
        # There is perhaps a better way of doing this but these two global vars are used in parsing the file.
        self.datapointnumber = 0;
        self.pointsofinterest = [];

    def printData(self):

        for x in self.f:
            self.datapointnumber = self.datapointnumber + 1
            temp = x.split()
            self.pointsofinterest.append([float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3])])

        self.plotData()


    # Thought process: used some of the previous position points to generate a spline, do a linear fit on that.
    #                  If the trajectory goes within a certain "box" relative to the bike at (0,0)
    #                  and if the velocity is not zero, send warning. Should we used accel as well?
    def plotData(self):
        plt.figure()
        plt.axis('off')
        axes = plt.gca()
        plt.text(0.5, 0.55, 'In the following demo:', horizontalalignment='center',verticalalignment = 'center')
        plt.text(0.5, 0.5, 'BLUE LINE = Measured Velocity Direction', horizontalalignment='center', verticalalignment='center',color='b')
        plt.text(0.5, 0.45, 'MAGENTA = Predicted Travel Direction', horizontalalignment='center', verticalalignment='center',color='m')
        plt.pause(10) # This line is technically not needed. I just used it for conveniently filming the screen.
        # plt.text(0.5,0.5," ")
        plt.axis('on')
        plt.pause(0.5)
        # The axes limit are set just for the current file being read. Adjust as needed.
        axes.set_xlim([-5, 10])
        axes.set_ylim([-1, 14])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.grid()
        # largestX is to track how far forwards a car has gone. If the car moves backwards, that is perhaps unrealistic
        # on the road, and hence it is being treated as noise for now.
        largestX = (self.pointsofinterest)[0][0]
        # To store data used for splining.
        count = 0
        curve = []
        numPoints = 9 # number of points used for the splining. Change as needed.
        for val,x in enumerate(self.pointsofinterest):
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10 ** 6}MB; Peak was {peak / 10 ** 6}MB")
            plt.pause(0.1)
            plt.cla() # Update the plot.
            axes.set_xlim([-5, 10])
            axes.set_ylim([-1, 14])
            plt.gca().set_aspect('equal', adjustable='box')
            plt.grid()
            # velocityvector is just used for human-friendliness, so can be changed as needed.
            velocityvector = 1.5

            # For valid points.
            if val==0 or (self.pointsofinterest)[val][0]>largestX:
                if count<numPoints:
                    curve.append([x[0],x[1]])
                largestX=(self.pointsofinterest)[val][0]
                speed=np.sqrt(x[2]**2+x[3]**2)
                velocities=plt.plot([x[0],x[0]+x[2]*speed*velocityvector],[x[1],x[1]+x[3]*speed*velocityvector],'b',alpha=0.5)
                plt.plot(x[0],x[1],'ro',markersize=15)

                # Splining code below. Previously attempted a polynomial fit, but that easily becomes overfitted.
                # Documentation: https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
                if count >=numPoints:
                    xvals = []
                    yvals = []
                    for loop in range(0,len(curve)):
                        xvals.append(curve[loop][0])
                        yvals.append(curve[loop][1])
                    xvals = np.array(xvals)
                    yvals = np.array(yvals)
                    xp = np.linspace(curve[0][0], curve[-1][0], num=100,endpoint=True)
                    # f2 = interp1d(xvals, yvals, kind='cubic')
                    f2 = CubicSpline(xvals, yvals, extrapolate=True)
                    # cubicsplinedata = plt.plot(xp, f2(xp), '--') # UNCOMMENT TO SHOW

                    # To update the "curve" list.
                    curve.pop(0)
                    curve.append([x[0],x[1]])

                    # Fitting the spline
                    model = LinearRegression().fit(xp.reshape((-1,1)),f2(xp))
                    # Documentation: https://realpython.com/linear-regression-in-python/
                    predictedpath = plt.plot([x[0],x[0]+2.5],[model.coef_[0]*x[0]+model.intercept_,model.coef_[0]*(x[0]+2.5)+model.intercept_],'m',alpha=0.5)

                    # These values will need to change based on the requirements.
                    front_pred = model.predict(np.linspace(0,10,20).reshape(-1,1))
                    rear_pred = model.predict(np.linspace(-10, 0, 20).reshape(-1, 1))
                    if x[0]<=0 and speed>0 and any(i <= 1 for i in rear_pred):
                        plt.title("REAR CAR WITH POSSIBLE REAR DANGER")
                    elif x[0] > 0 and speed > 0 and any(i <= 1 for i in rear_pred):
                        plt.title("FRONTAL CAR WITH POSSIBLE REAR DANGER")
                    elif x[0]>0 and speed>0 and any(i <= 1 for i in front_pred):
                        plt.title("FRONTAL CAR WITH POSSIBLE FRONTAL DANGER")
                    elif x[0] <= 0 and speed > 0 and any(i <= 1 for i in front_pred):
                        plt.title("REAR CAR WITH POSSIBLE FRONTAL DANGER")
                    else:
                        plt.title("YOU ARE SAFE")
                    # plt.legend((velocities, predictedpath, cubicsplinedata), (
                    # 'Measured Velocity', 'Predicted Path', "Car's Previous Trajectory"))
                count = count + 1 # Housekeeping
                bike = plt.plot(0,0,"k>",markersize=25)
                plt.ylabel('Distance to your left (meters)')
                plt.xlabel('Direction of travel (meters)')

                plt.show

            # For all the "bad" points.
            else:
                plt.plot(x[0], x[1], 'ko', alpha=0.1,markersize=15)
                bike = plt.plot(0, 0, "k>",markersize=25)
                # plt.legend(bike,"Placement of Bike")
                plt.ylabel('Distance to your left (meters)')
                plt.xlabel('Direction of travel (meters)')
                plt.title("Noisy data")
                # print(val)
                plt.show

# Begin.
if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d, CubicSpline
    from sklearn.linear_model import LinearRegression
    import tracemalloc

    tracemalloc.start()


    # Using IWR1843 data "tm_demo_uart_stream_5_0-40_both_bundary_gating_2020novel"
    logansdata=showData()
    logansdata.printData()

    tracemalloc.stop()