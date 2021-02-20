# Creator: Charlie Zhu, Jan 2021
# Purpose: Virtual prototype to show concept of classification of the right-hook collision types.
# Notes: Still needs major re-work. Please rename variables as you see fit.

class showData(object):

    def __init__(self):
        # Taking file of data generated from the IWR1843 MATLAB GUI, which has been re-parsed in a text file.
        # Format of the text file: X pos, Y pos, X vel, Y vel, X acc, Y acc
        self.f = open("AWR_Frontal_hook_data.txt", "r")
        # There is perhaps a better way of doing this but these two global vars are used in parsing the file.
        self.datapointnumber = 0;
        self.pointsofinterest = [];
        #for lowpass filter
        self.T = 22*0.03  # Sample Period
        self.fs = 1/0.03 # sample rate, Hz
        self.cutoff = 10  # desired cutoff frequency of the filter, Hz
        self.nyq = 0.5 * self.fs  # Nyquist Frequency
        self.order = 2  # sin wave can be approx represented as quadratic
        self.n = int(self.T * self.fs)  # total number of samples
        #for boundarybox
        self.bound_left = -7.5
        self.bound_up = 3
        self.bound_right = 7.5
        #for alarm
        self.danger = False
        #for spline
        self.extrapolated_length = 30
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
        return [newX, newY]
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
    # Thought process: used some of the previous position points to generate a spline, do a linear fit on that.
    #                  If the trajectory goes within a certain "box" relative to the bike at (0,0)
    #                  and if the velocity is not zero, send warning. Should we used accel as well?
    def plotData(self):
        '''
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
        '''
        # largestX is to track how far forwards a car has gone. If the car moves backwards, that is perhaps unrealistic
        # on the road, and hence it is being treated as noise for now.
        largestX = (self.pointsofinterest)[0][0]
        # To store data used for splining.
        count = 0
        curve = []
        #velocity = []
        numPoints = 22 # number of points used for the splining. Change as needed.
        num_split = 3
        for val,x in enumerate(self.pointsofinterest):

            #current, peak = tracemalloc.get_traced_memory()
            #print(f"Current memory usage is {current / 10 ** 6}MB; Peak was {peak / 10 ** 6}MB")
            '''
            plt.pause(0.05)
            plt.cla() # Update the plot.
            axes.set_xlim([-10, 10])
            axes.set_ylim([-10, 20])
            plt.gca().set_aspect('equal', adjustable='box')
            plt.grid()
            '''
            # velocityvector is just used for human-friendliness, so can be changed as needed.
            #velocityvector = 1.5

            # For valid points.
            if 1:# val == 0 or (self.pointsofinterest)[val][0] > largestX:
                start_mem = hpy().heap().size
                if count < numPoints:
                    newXY = self.coord_rotation2D(-30, [x[0],x[1]])
                    curve.append([newXY[0],newXY[1]])
                    #velocity.append([x[2], x[3]])
                #largestX = (self.pointsofinterest)[val][0]
                #speed=np.sqrt(x[2]**2+x[3]**2)
                #velocities = plt.plot([x[0],x[0]+x[2]*speed*velocityvector],[x[1],x[1]+x[3]*speed*velocityvector],'b',alpha=0.5)

                # Splining code below. Previously attempted a polynomial fit, but that easily becomes overfitted.
                # Documentation: https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
                if count >=numPoints:
                    start_time = time.time()
                    xvals = []
                    yvals = []

                    curve.pop(0)
                    newXY = self.coord_rotation2D(-30, [x[0], x[1]])
                    curve.append([newXY[0],newXY[1]])
                    #may use velocity for weighted average, right now it is not implemented
                    #velocity.pop(0)
                    #velocity.append([x[2], x[3]])
                    ave_point = []
                    all_ave_data = []

                    self.point_average(curve, ave_point, num_split)
                    self.point_average(curve, all_ave_data, num_split)
                    #plt.plot(curve[-1][0], curve[-1][1], 'ro', markersize=10)
                    #plt.gca().add_patch(Rectangle((self.bound_left, -1*self.bound_up/2),
                                                  #self.bound_right - self.bound_left, self.bound_up,
                                                  #edgecolor='red',
                                                  #facecolor='none',
                                                  #lw=2))
                    for loop in range(0,len(curve)):
                        xvals.append(curve[loop][0])
                        yvals.append(curve[loop][1])
                        #xvals.append(ave_point[loop][0])
                        #yvals.append(ave_point[loop][1])
                    xvals = np.array(xvals)
                    yvals = np.array(yvals)

                    #filtering results
                    filx = self.butter_lowpass_filter(xvals, self.cutoff, self.fs, 6)
                    fily = self.butter_lowpass_filter(yvals, self.cutoff, self.fs, 6)
                    #doing the least squared univariate spline
                    fx = self.spline_filter(filx, 15)
                    fy = self.spline_filter(fily,15)
                    for i in range(0,self.extrapolated_length):
                        if (fx(i+xvals[-1]) <= self.bound_right and fx(i+xvals[-1]) >= self.bound_left) \
                                and (fy(i+xvals[-1]) <= self.bound_up / 2
                                     and fy(i+xvals[-1]) >= -1 * self.bound_up / 2):
                            self.danger = True
                            break
                        else:
                            self.danger = False

                    xp = np.linspace(xvals[0], xvals[-1]+self.extrapolated_length, 100, endpoint=True)
                    #fxy = self.spline_filter_xy(xvals, yvals, 2)
                    #plt.plot(fx(xp), fy(xp), 'g', lw=3)
                    #if (self.danger):
                    #    plt.text(0, 22, 'DANGER', horizontalalignment='center', fontsize='large', color='red')
                    #else:
                    #    plt.text(0, 22, 'SAFE  ', horizontalalignment='center', fontsize='large', color='green')
                    #print (fx(xvals[-1]))
                    # using polyfit
                    '''
                    poly = PolynomialFeatures(degree=4)
                    X_poly = poly.fit_transform(X)
                    poly.fit(X_poly, yvals)
                    '''
                    #plt.plot(xvals, yvals, 'b', lw=4)
                    #plt.plot(xp, fxy(xp), 'g', lw=3)

                    # time taken & memory used
                    end_mem = hpy().heap().size
                    end_time = time.time()
                    print("time taken is: " + str(end_time - start_time) + "; memory usage is: " + str(end_mem - start_mem))
                count = count + 1 # Housekeeping
                #bike = plt.plot(0,0,"k>",markersize=25)
                #plt.ylabel('Distance to your left (meters)')
                #plt.xlabel('Direction of travel (meters)')

                #plt.show

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
    from scipy.interpolate import interp1d, CubicSpline, LSQBivariateSpline, LSQUnivariateSpline, UnivariateSpline
    from sklearn.linear_model import LinearRegression
    import tracemalloc
    import bspline
    import bspline.splinelab as splinelab
    from scipy.signal import butter, filtfilt
    import math
    from sklearn.preprocessing import PolynomialFeatures
    from matplotlib.patches import Rectangle
    import time
    from guppy import hpy

    #tracemalloc.start()


    # Using IWR1843 data "tm_demo_uart_stream_5_0-40_both_bundary_gating_2020novel"
    start_mem_all = hpy().heap().size
    logansdata=showData()
    logansdata.printData()
    print("all memory used is: " + str(hpy().heap().size - start_mem_all))

    tracemalloc.stop()