class showData(object):

    def __init__(self):
        # Taking file of data generated from the IWR1843 MATLAB GUI, which has been re-parsed in a text file.
        # Format of the text file: X pos, Y pos, X vel, Y vel, X acc, Y acc
        self.f = open("AWR_Frontal_hook_data.txt", "r")
        # self.f = open("moving_AWR_rear_hook_data.txt", "r")

        # There is perhaps a better way of doing this but these two global vars are used in parsing the file.
        self.datapointnumber = 0
        self.pointsofinterest = []
        self.cov = np.array([4.2187, 1.3983, 0.33535, 0.28941])
        self.xy_data = []
        self.past_path = []
        self.frame_period = 0.08
        self.ProcessNoise = 1.0
        # for lowpass filter
        self.T = 22 * self.frame_period  # Sample Period
        self.fs = 1 / self.frame_period  # sample rate, Hz
        self.cutoff = 10  # desired cutoff frequency of the filter, Hz
        self.nyq = 0.5 * self.fs  # Nyquist Frequency
        self.order = 2  # sin wave can be approx represented as quadratic
        self.n = int(self.T * self.fs)  # total number of samples
        #for boundarybox
        self.bound_left = -9
        self.bound_up = 0.6
        self.bound_right = 8.5
        self.bound_rear = -2.6
        #for alarm
        self.danger = 'safe'
    def spline_filter(self, data, nsegs):
        """Detrend a possibly periodic timeseries by fitting a coarse piecewise
           smooth cubic spline

        Parameters
        ----------
        data : ndarray
            list of observed values
        nsegs : number
            number of spline segments

        Returns
        -------
        filtered : ndarray

        """
        index = np.arange(len(data))
        nknots = max(2, nsegs + 1)
        knots = np.linspace(index[0], index[-1], nknots + 2)[1:-2]
        #print("length of data is:" + str(len(data)) + " length of index is: " + str(len(index)))
        #print(data)
        return LSQUnivariateSpline(index, data, knots[-1:1], ext='extrapolate')
    def coord_rotation2D(self, degree, data_pair):
        newX = data_pair[0]*math.cos(math.radians(degree)) - data_pair[1]*math.sin(math.radians(degree))
        newY = data_pair[0] * math.sin(math.radians(degree)) + data_pair[1] * math.cos(math.radians(degree))
        newVX = data_pair[2] * math.cos(math.radians(degree)) - data_pair[3] * math.sin(math.radians(degree))
        newVY = data_pair[2] * math.sin(math.radians(degree)) + data_pair[3] * math.cos(math.radians(degree))
        return [newX, newY, newVX, newVY]
    def spline_filter_xy(self, datax, datay, nsegs):
        """Detrend a possibly periodic timeseries by fitting a coarse piecewise
           smooth cubic spline

        Parameters
        ----------
        data : ndarray
            list of observed values
        nsegs : number
            number of spline segments

        Returns
        -------
        filtered : ndarray

        """
        index = np.arange(len(datax))
        nknots = max(2, nsegs + 1)
        knots = np.linspace(index[0], index[-1], nknots + 2)[1:-2]
        #print("length of data is:" + str(len(data)) + " length of index is: " + str(len(index)))
        #print(data)
        return LSQUnivariateSpline(datax, datay, knots[-1:1], ext='extrapolate')

    def butter_lowpass_filter(self, data, cutoff, fs, order):
        normal_cutoff = cutoff / (fs*0.5) #nyq freq
        # Get the filter coefficients
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data)
        return y

    def printData(self):

        for x in self.f:
            self.datapointnumber = self.datapointnumber + 1
            temp = x.split()
            self.pointsofinterest.append([float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3])])
        self.plotData()

    def point_average(self, curve, ave_point, num_split):
        for i in range(0, len(curve)//num_split):
            temp_vel = []

            ave_point.append(np.mean([curve[i], curve[i+1], curve[i+2], curve[i+3]], axis=0))
    def point_average_points(self, curve, num_split, velocities):
        ave_point = []
        ave_point.append(np.mean([curve[-1], curve[-2], curve[-3]], axis=0))
        return ave_point
    # Thought process: used some of the previous position points to generate a spline, do a linear fit on that.
    #                  If the trajectory goes within a certain "box" relative to the bike at (0,0)
    #                  and if the velocity is not zero, send warning. Should we used accel as well?
    def plotData(self):
        #memory usage, not used here
        #start_mem = hpy().heap().size
        plt.figure()
        plt.axis('off')
        axes = plt.gca()
        plt.text(0.5, 0.55, 'In the following demo:', horizontalalignment='center',verticalalignment = 'center')
        plt.text(0.5, 0.5, 'BLUE LINE = Measured Velocity Direction', horizontalalignment='center', verticalalignment='center',color='b')
        plt.text(0.5, 0.45, 'MAGENTA = Predicted Travel Direction', horizontalalignment='center', verticalalignment='center',color='m')
        plt.pause(1) # This line is technically not needed. I just used it for conveniently filming the screen.
        # plt.text(0.5,0.5," ")
        plt.axis('on')
        #plt.pause(0.5)
        # The axes limit are set just for the current file being read. Adjust as needed.

        plt.gca().set_aspect('equal', adjustable='box')
        plt.grid()
        # largestX is to track how far forwards a car has gone. If the car moves backwards, that is perhaps unrealistic
        # on the road, and hence it is being treated as noise for now.
        #get memory usage

        #largestX = (self.pointsofinterest)[0][0]
        #largestY = (self.pointsofinterest)[0][1]
        # To store data used for splining.
        count = 0
        curve = []
        velocity = []
        numPoints = 12
        num_split = 3
        start_flag = 1
        #validpoint = True
        #pastvel = [np.sign(1),np.sign(1)]
        for val,x in enumerate(self.pointsofinterest):

            #current, peak = tracemalloc.get_traced_memory()
            #print(f"Current memory usage is {current / 10 ** 6}MB; Peak was {peak / 10 ** 6}MB")
            plt.pause(0.05)
            plt.cla() # Update the plot.
            axes.set_xlim([-10, 10])
            axes.set_ylim([-5, 5])
            plt.gca().set_aspect('equal', adjustable='box')
            plt.grid()
            # velocityvector is just used for human-friendliness, so can be changed as needed.
            velocityvector = 1.5

            # For valid points.
            #point_difference = [np.sign((self.pointsofinterest)[val][0] - largestX), np.sign((self.pointsofinterest)[val][0] - largestY)]

            if 1: #val == 0 or (point_difference[0] == pastvel[0] and point_difference[1] == pastvel[1]):
                if count < numPoints:
                    newXY = self.coord_rotation2D(-30, [x[0],x[1],x[2], x[3]])
                    curve.append([newXY[0],newXY[1], newXY[2],newXY[3]])
                    #velocity.append([x[2], x[3]])
                    #pastvel = [np.sign((self.pointsofinterest)[val][2]), np.sign((self.pointsofinterest)[val][3])]
                    if start_flag == 1:
                        kf = Kalman(curve[0][0], curve[0][1], curve[-1][2], curve[-1][3], self.ProcessNoise)
                        start_flag -= 1
                    elif start_flag == 0:
                        accelx = (curve[-1][2] - curve[-2][2]) / self.frame_period
                        accely = (curve[-1][3] - curve[-2][3]) / self.frame_period
                        kf.predict(self.frame_period, 3, [accelx, accely, accelx, accely])
                        kf.update(curve[-1], self.cov, curve)
                        #start_flag -= 1
                    self.xy_data.append([kf._x[0], kf._x[1]])
                    # print(self.xy_data)
                #largestX = (self.
                # pointsofinterest)[val][0]
                #largestY = (self.pointsofinterest)[val][1]

                speed=np.sqrt(x[2]**2+x[3]**2)
                #velocities = plt.plot([x[0],x[0]+x[2]*speed*velocityvector],[x[1],x[1]+x[3]*speed*velocityvector],'b',alpha=0.5)

                # Splining code below. Previously attempted a polynomial fit, but that easily becomes overfitted.
                # Documentation: https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
                if count >=numPoints:
                    start_time = time.time()
                    #xvals = []
                    #yvals = []

                    curve.pop(0)
                    self.xy_data.pop(0)
                    newXY = self.coord_rotation2D(-30, [x[0],x[1],x[2], x[3]])
                    curve.append([newXY[0],newXY[1], newXY[2],newXY[3]])
                    #if count > numPoints:
                    #    kf.update(curve[1], self.cov)
                    self.past_path.append(kf._x)
                    #may use velocity for weighted average, right now it is not implemented
                    #velocity.pop(0)
                    #velocity.append([x[2], x[3]])
                    ave_point = []
                    all_ave_data = []

                    self.point_average(curve, ave_point, num_split)
                    #print(ave_point)
                    #self.point_average(curve, all_ave_data, num_split, velocity)
                    plt.plot(ave_point[-1][0], ave_point[-1][1], 'ro', markersize=10)
                    accelx = (ave_point[- 1][2] - ave_point[- 2][2])/self.frame_period
                    accely = (ave_point[- 1][3] - ave_point[- 2][3])/self.frame_period
                    kf.predict(self.frame_period, 5, [accelx, accely, accelx, accely])
                    predicted_path = kf.update(ave_point[-1], self.cov, curve)
                    warn = (WarningAlgo())
                    self.danger = warn.scan_data(predicted_path)
                    self.xy_data.append([kf._x[0],kf._x[1]])
                    # for i in range(len(predicted_path)):
                    #     if (predicted_path[i][0] <= self.bound_right and predicted_path[i][0] >= self.bound_left) \
                    #             and (predicted_path[i][1] <= self.bound_up/2
                    #                  and predicted_path[i][1] >= -1*self.bound_up/2):
                    #         self.danger = True
                    #         break
                    #     else:
                    #         self.danger = False
                    '''
                    for loop in range(0, len(self.xy_data)):
                        xvals.append(self.xy_data[loop][0])
                        yvals.append(self.xy_data[loop][1])
                        # xvals.append(ave_point[loop][0])
                        # yvals.append(ave_point[loop][1])
                    xvals = np.array(xvals)
                    yvals = np.array(yvals)
                    '''
                    # doing the least squared univariate spline
                    '''
                    if (count >= 15):
                        fx = self.spline_filter(xvals, 5)
                        fy = self.spline_filter(yvals, 5)
                        xp = np.linspace(xvals[-1], xvals[-1] + 50, 100, endpoint=True)
                        plt.plot(fx(xp), fy(xp), 'r', lw=2)
                    '''
                    # fxy = self.spline_filter_xy(xvals, yvals, 2)


                    #print([i[0] for i in predicted_path])
                    #print("predicted x: " + str(predicted_path[98][0]) + "predicted y: " + str(predicted_path[98][1]))
                    #print("actual x: " + str(curve[-1][0]) + "actual y: " + str(curve[-1][1]))

                    #time taken & memory used
                    #end_mem = hpy().heap().size
                    end_time = time.time()
                    #print("time taken is: " + str(end_time - start_time) + "; memory usage is: " + str(end_mem-start_mem))
                    plt.plot([i[0] for i in predicted_path], [i[1] for i in predicted_path], 'g', lw=3)
                    plt.plot([i[0] for i in self.past_path], [i[1] for i in self.past_path], 'b', lw=2)
                    plt.gca().add_patch(Rectangle((self.bound_left, -1*self.bound_up),
                                                  self.bound_right - self.bound_left, 2*self.bound_up,
                                                  edgecolor='red',
                                                  facecolor='none',
                                                  lw=2))
                    plt.plot([self.bound_rear,self.bound_rear],[self.bound_up,-self.bound_up],'r')

                    if (self.danger == "lateral"):
                        plt.title('REAR SIDE DANGER', horizontalalignment='center', fontsize = 'large', color = 'red')
                    elif (self.danger == "rear"):
                        plt.title('REAR HOOK DANGER', horizontalalignment='center', fontsize='large', color='red')
                    elif (self.danger == "front"):
                        plt.title('FRONT HOOK DANGER', horizontalalignment='center', fontsize='large', color='red')
                    else:
                        plt.title('SAFE', horizontalalignment='center', fontsize = 'large', color = 'green')

                    '''
                    #filtering results
                    filx = self.butter_lowpass_filter(xvals, self.cutoff, self.fs, 6)
                    fily = self.butter_lowpass_filter(yvals, self.cutoff, self.fs, 6)
                    '''
                count = count + 1 # Housekeeping
                plt.plot(0,0,"k>",markersize=15) #where the bike is
                plt.ylabel('Distance to your left (meters)')
                plt.xlabel('Direction of travel (meters)')

                plt.show

            # For all the "bad" points.
            else:
                plt.plot(x[0], x[1], 'ko', alpha=0.1,markersize=15)
                bike = plt.plot(0, 0, "k>",markersize=15)
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
    from scipy.interpolate import interp1d, CubicSpline, LSQBivariateSpline, LSQUnivariateSpline, UnivariateSpline
    from sklearn.linear_model import LinearRegression
    import tracemalloc
    # import bspline
    # import bspline.splinelab as splinelab
    from scipy.signal import butter, filtfilt
    import math
    from sklearn.preprocessing import PolynomialFeatures
    from Kalman_filter import Kalman
    from matplotlib.patches import Rectangle
    import time
    from WarningAlgo import WarningAlgo
    # from guppy import hpy

    # Using IWR1843 data "tm_demo_uart_stream_5_0-40_both_bundary_gating_2020novel"
    logansdata=showData()
    logansdata.printData()

    tracemalloc.stop()