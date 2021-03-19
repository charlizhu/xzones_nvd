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
        self.filename = "rear_righthook.txt"
        self.writefile = open("temp.txt", 'a+')
        self.f = open(self.filename, "r")
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
        # previous pos and vel
        self.lx_p = []
        self.ly_p = []
        self.vx_p = []
        self.vy_p = []

        self.xpoints = []
        self.xpoints.append({
            "tid": 1,
            "pos": [],
        })
        self.ypoints = []
        self.ypoints.append({
            "tid": 1,
            "pos": [],
        })
        self.previous_tracked = []
        self.threadhold_error =0.8
        self.time_stamp = 0.08
        self.new_tid_num = 0
        self.tid_write = [None] * 25
        self.detectItems()

        #for colors only
        self.colors = None
        self.by_hsv = None
        self.sorted_names = None
        self.color_start_index = 15
        self.prev_points = 3


    def add_new_tid(self):
        self.new_tid_num += 1
        return self.new_tid_num
    def coord_rotation2D(self, degree, data_pair):
        newX = data_pair[0]*math.cos(math.radians(degree)) - data_pair[1]*math.sin(math.radians(degree))
        newY = data_pair[0] * math.sin(math.radians(degree)) + data_pair[1] * math.cos(math.radians(degree))
        newVX = data_pair[2] * math.cos(math.radians(degree)) - data_pair[3] * math.sin(math.radians(degree))
        newVY = data_pair[2] * math.sin(math.radians(degree)) + data_pair[3] * math.cos(math.radians(degree))
        return [newX, newY, newVX, newVY]

    def parsing_objects_line (self, datapoint, carIndex):
        xLoc = float(datapoint[carIndex + 2]) + float(datapoint[carIndex + 3]) * 256
        yLoc = float(datapoint[carIndex + 4]) + float(datapoint[carIndex + 5]) * 256
        if xLoc > 32767:
            xLoc = xLoc - 65536
        if yLoc > 32767:
            yLoc = yLoc - 65536
        xp = xLoc / 128
        yp = yLoc / 128
        xVel = float(datapoint[carIndex + 6]) + float(datapoint[carIndex + 7]) * 256
        yVel = float(datapoint[carIndex + 8]) + float(datapoint[carIndex + 9]) * 256
        if xVel > 32767:
            xVel = xVel - 65536
        if yVel > 32767:
            yVel = yVel - 65536
        xv = xVel / 128
        yv = yVel / 128
        return xp, yp, xv, yv

    def screening(self, RotXY):
        if RotXY[1] < -5 or RotXY[0] < -5 or RotXY[1] > 15 or RotXY[3] < -3 or RotXY[0] > 8:
            return -1
        else:
            return 1

    def color_array(self):
        self.colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        self.by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                        for name, color in self.colors.items())
        self.sorted_names = [name for hsv, name in self.by_hsv]

    def detectItems(self): # In reality this would be reading from the UART.
        self.color_array()
        for cars in self.f:
            first_size, first_peak = tracemalloc.get_traced_memory()
            stop = timeit.default_timer()
            #print(stop - self.start)
            self.start = stop
            #print(first_size / 10 ** 6, first_peak / 10 ** 6)
            datapoint = cars.split()
            #for each "object" in one line of file
            #need to wait until 5-6 iterations has passed for variance
            #print("datapoint insides: " + str(datapoint))
            for carIndex in range(0, len(datapoint), 14):

                # This next snippet is taken from the TI MATLAB GUI script
                #parsing
                xp, yp, xv, yv = self.parsing_objects_line(datapoint, carIndex)
                #rotation matrix
                RotXY = self.coord_rotation2D(-30, [xp, yp, xv, yv])
                #print(RotXY[0],RotXY[1])
                if not(self.screening(RotXY) == 1): # to be modified, need more conditions here...
                    # print(xLocRot,yLocRot)
                    continue
                # if two appear at the same time then pass
                # appending to lists for comparison for tracking
                self.tidcurrent.append(datapoint[carIndex])
                if len(self.tidcurrent) != len(set(self.tidcurrent)):
                    self.tidcurrent.pop()
                #print("tid currently appended", self.tidcurrent)
                print("rotated xy is at FIRST:", RotXY)

                self.lx.append(RotXY[0])
                self.ly.append(RotXY[1])
                self.vx.append(RotXY[2])
                self.vy.append(RotXY[3])
                #print("list lx!!!!!!!!!!!!!!!!! ", self.lx)
            #if len(self.lx) < self.prev_points:
            #        continue
            if (not self.tidcurrent):
                pass
            # print("datapint is here:" + str(datapoint))
            self.trackItems()
            if (self.currenttidrepeats is not None) and (self.currenttid is not None):
                # print(self.currenttid, (self.tidcurrent).index(self.currenttid))
                self.write_to_file()
            else:
                self.writefile.write("000" + '\n')
            self.lx = []
            self.ly = []
            self.vx = []
            self.vy = []
            self.tidcurrent = []
        print(self.xpoints)
        # plot all the objects
        plt.figure()
        for i in range(len(self.xpoints)):
            tid_colour = (self.xpoints[i]['tid'])

            plt.plot(self.xpoints[i]['pos'], self.ypoints[i]['pos'],
                         color=self.colors[self.sorted_names[tid_colour + 15]])
            plt.annotate(str(tid_colour), xy = (self.xpoints[i]['pos'][0], self.ypoints[i]['pos'][0]))
            #print(str(self.xpoints[i][j]['pos']) + " " + str(self.ypoints[i][j]['pos']))
        plt.xlim((-5,10))
        plt.grid('on')
        plt.show()


    '''
    percentage error of predicted position: 
    (current_pos - (previous pos + velocity*time_passed))/current_pos
    '''
    def percentage_error(self, current_list, index1, previous_list, index2, velocity, time_mult):
        return abs(abs(current_list[index1] - (previous_list[index2]+self.time_stamp*time_mult*velocity[index2]))
                   / current_list[index1])

    def tid_write_switch(self, index_track):
        self.tid_write.pop(index_track)
        self.tid_write.append(None)

    def variance_screening(self, x_val, y_val, num):
        return

    def movement_cap(self):
        return
    # tracking the object with tolerance of the previous pos
    def tracking_with_tolerance(self):
        #print ("current tid before toloerance runs: " + str(self.tidcurrent))
        # if the current tid is empty, then nothing is tracked, and so clear the list
        if not self.tidcurrent:
            self.tidtracked.clear()
            self.tid_write = [None] * 25
        for tid in self.tidcurrent:
            if (tid in self.tidprevious[-1]) and (len(self.tidprevious[-1]) <= len(self.tidcurrent)):
                if tid not in self.tidtracked:
                    # adding variance condition later
                    self.tidtracked.append(tid)
                    # this means it is a newly tracked object, so append (tid_write)
                    self.tid_write[self.tidtracked.index(tid)] = self.add_new_tid()
            else:
                #print("not in previous ,and tracked list is currently: " + str(self.tidtracked))
                #print("previous tid is: " + str(self.tidprevious))

                for tid_prev_list in range(0, len(self.tidprevious)):
                    tid_count = 0
                    prev_index = -tid_prev_list - 1
                    for tid_prev in self.tidprevious[prev_index]:
                        # get the percentage error  of predicted pos vs actual pos
                        print("prev_index######################", prev_index)
                        x_err = self.percentage_error(self.lx, self.tidcurrent.index(tid), self.lx_p[prev_index],
                                                      self.tidprevious[prev_index].index(tid_prev), self.vx_p[prev_index], -prev_index)
                        y_err = self.percentage_error(self.ly, self.tidcurrent.index(tid), self.ly_p[prev_index],
                                                      self.tidprevious[prev_index].index(tid_prev), self.vy_p[prev_index], -prev_index)
                        # if predicted is within the error threshold, it means it is a vehicle to be tracked, add to the
                        # tracked list
                        print("tid is at: " + str(tid) + " previous tid_iter is currently at: " + str(tid_prev) + " x err and y err are: " + str(x_err) + ", " + str(y_err) )
                        if x_err <= self.threadhold_error and y_err <= self.threadhold_error:
                            print("less than error")
                            if 1: #tid not in self.tidtracked:  #maybe some logic error here?
                                if tid not in self.tidtracked:
                                    self.tidtracked.append(tid)
                                # print("current tid added: " + str(tid))

                                # if tru(tid prev is in tracked), then it is a previously tracked object (tid_write)
                                print("tis track here in bool", self.tidtracked)
                                if self.tidtracked and (tid_prev != tid and (tid_prev in self.tidtracked)):
                                    print("prev tracked POPPED: " + str(tid_prev))
                                    index_track = self.tidtracked.index(tid_prev)
                                    self.tidtracked.pop(index_track)
                                    # need to check this part
                                    self.tid_write[self.tidtracked.index(tid)] = self.tid_write[index_track]
                                    if index_track != self.tidtracked.index(tid):
                                        self.tid_write_switch(index_track)
                                # if not then it is a new object tracked(tid_write)
                                elif tid_prev not in self.tidtracked:  # if tid_prev not in self.tidtracked:
                                    self.tid_write[self.tidtracked.index(tid)] = self.add_new_tid()
                            tid_count += 1
                        else:
                            # also pop whatever is in new tracked list (tid_write)
                            if (tid_prev not in self.tidcurrent) and (tid_prev in self.tidtracked) \
                                    and (-prev_index >= len(self.tidprevious)):
                                #print("prev tracked is: " + str(self.tidprevious))
                                print("POPPED because of out of torlrance: ", tid_prev)
                                index_track = self.tidtracked.index(tid_prev)
                                self.tidtracked.pop(index_track)
                                self.tid_write_switch(index_track)
                    # if predicted is not within the threshold, then pop the component if recorded in the tracked list
                    if tid_count == 0:
                        #print("tid count is 0")
                        # also pop whatever is in new tracked list  (tid_write)
                        if self.tidtracked and (tid in self.tidtracked):
                            index_track = self.tidtracked.index(tid)
                            print("POPPED because tid-count is 0: ", tid)
                            self.tidtracked.pop(index_track)
                            self.tid_write_switch(index_track)
                    else:
                        # get rid of excess points and exit
                        index_pop = [index for (index, item) in enumerate(self.tidtracked)
                                     if item not in self.tidcurrent]
                        for pop in range(len(index_pop)-1, -1, -1):
                            self.tidtracked.pop(index_pop[pop])
                            self.tid_write_switch(pop)
                        break
        return

    def update_previous(self, start):
        if start:
            self.tidprevious.append(self.tidcurrent)
            (self.lx_p).append(self.lx)
            (self.ly_p).append(self.ly)
            (self.vx_p).append(self.vx)
            (self.vy_p).append(self.vy)
        else:
            self.tidprevious.append(self.tidcurrent)
            self.tidprevious.pop(0)
            (self.lx_p).append(self.lx)
            self.lx_p.pop(0)
            (self.ly_p).append(self.ly)
            self.ly_p.pop(0)
            (self.vx_p).append(self.vx)
            self.vx_p.pop(0)
            (self.vy_p).append(self.vy)
            self.vy_p.pop(0)
    '''
    compares previous and current timestamp of objects to screen out shits
    '''
    def trackItems(self):
        # previous_tracked
        if len(self.tidprevious) < 3:
            # update previous points
            self.update_previous(True)
            pass
        else:
            self.tracking_with_tolerance()
            # update previous points
            self.update_previous(False)
        # print(self.tidtracked, self.tidcounts)

        if len(self.tidtracked) >= 1:
            # print("tidtracked not empty, it is: " + str(self.tidtracked))
            self.currenttid = self.tidtracked[0]
            self.currenttidrepeats = self.tidcounts
        else:
            self.currenttid = None
            self.currenttidrepeats = None
        print("tracked points:" + str(self.tidtracked))
    def write_to_file(self):
        print("Logan's Kalman filter script")
        print("current point is: " + str(self.tidcurrent))
        print("tid to write is: " + str(self.tid_write))
        print("tid track right now after print write is: " + str(self.tidtracked))
        if not self.tidcurrent:
            return
        print("Position ", self.lx, self.ly)
        print("Velocity ", self.vx, self.vy)
        # write all tracked data points to a file
        for i in range(0, len(self.tid_write)):
            id_cnt = 0
            if self.tid_write[i] != None:
                self.writefile.write(str(len(self.tidtracked)) + ' ' + str(self.tid_write[i]) + ' ' +
                                     str(self.lx[self.tidcurrent.index(self.tidtracked[i])]) + ' ' +
                                     str(self.ly[self.tidcurrent.index(self.tidtracked[i])]) + ' ' +
                                     str(self.vx[self.tidcurrent.index(self.tidtracked[i])]) + ' ' +
                                     str(self.vy[self.tidcurrent.index(self.tidtracked[i])]) + ' ')
                for id in range(len(self.xpoints)):
                    if self.xpoints[id]['tid'] == self.tid_write[i]:
                        self.xpoints[id]['pos'].append(self.lx[self.tidcurrent.index(self.tidtracked[i])])
                        self.ypoints[id]['pos'].append(self.ly[self.tidcurrent.index(self.tidtracked[i])])
                        id_cnt += 1
                if id_cnt == 0:
                    self.xpoints.append({
                        "tid": self.tid_write[i],
                        "pos": [self.lx[self.tidcurrent.index(self.tidtracked[i])]],
                    })
                    self.ypoints.append({
                        "tid": self.tid_write[i],
                        "pos": [self.ly[self.tidcurrent.index(self.tidtracked[i])]],
                    })

        self.writefile.write('\n')



if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import tracemalloc
    import timeit
    from matplotlib import colors as mcolors
    tracemalloc.start()
    choosePoint()