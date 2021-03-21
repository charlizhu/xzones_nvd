class WarningAlgo(object):

    def __init__(self):

        self.rear_bounds = [-2.6,0]
        self.front_bounds = [0,8.5]
        self.side_bounds = [0,0.6]
        self.rear_sides = [-9,0]
        self.turning_tolerance = 0.035
        self.backwards_tolerance = 0.5 # percentage of the predicted path's points are going "behind" the car.
        self.xdir_tolerance = 0.4 # tolerance in the x direction
        self.noise_tolerance = 10
        # self.scan_data()

    def scan_data(self,predicted_path):
        # First case is where the car is approaching too close to the side.
        # This is just to show that the algorithm can have other scenarios too
        # such as the rear-end head-on collision, for after the MVP.

        pp = predicted_path[1:]
        cl = predicted_path[0] # The first item in predicted_path is the current_location (i.e.: cl).
        if (self.side_bounds[0] < cl[1] < self.side_bounds[1]) \
                and (self.rear_sides[0] < cl[0] < self.rear_sides[1]) \
                and ((pp[pp.index(min([row[1] for row in pp]))][1] - cl[1]) < self.turning_tolerance): # this line might need fixing
            return "lateral"
        else:
            count = 0
            # Now it checks for side hook conditions.
            for point in pp:
                # if sum(i[0] < (cl[0]-self.xdir_tolerance) for i in pp)/len(pp) > self.backwards_tolerance:
                #     return "skip"
                # Rear side hook.
                if (self.rear_bounds[0] < point[0] < self.rear_bounds[1]):
                    if (point[1] < self.side_bounds[1]):
                        # print(self.cl, point, count)
                        return "rear"
                # Front side hook.
                if (self.front_bounds[0] < point[0] < self.front_bounds[1]):
                    if (point[1] < self.side_bounds[1]):
                        # print(self.cl, point, count)
                        return "front"
                count = count + 1
            # If no warnings were previously raised, then the cyclist is safe.
            return "safe"

    def checkNoise(self,predicted_path,previous_path):
        if abs(predicted_path[-1][1] - previous_path[-1][1]) > self.noise_tolerance:
            return 'noisy'
        return 'clean'

# Modify the Main function if the algorithm is to directly interface with the Kalman script.
if __name__ == "__main__":
    import tracemalloc
    import timeit

    tracemalloc.start()
    start_time = timeit.default_timer()

    numPoints = 6
    count = 0
    predicted_path = []

    # Inputs: car's predicted path, and current location.
    # IF I UNDERSTAND CORRECTLY: predicted_path is a n-by-2 matrix of [x,y] positions.
    myFile = open("sim_data\one_car_rear_sidehook.txt",'r') # change as needed
    for points in myFile:
        current_line = points.split()
        predicted_path.append([float(x) for x in current_line])
        if count >= numPoints:
            x = WarningAlgo(predicted_path)
            print(x.scan_data())
            predicted_path.pop(0)
        count = count + 1
        memory_size, memory_peak = tracemalloc.get_traced_memory()
        run_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        print("Memory Usage (MB): " + str(memory_peak/(10**6)) + " Run Time (sec): " + str(run_time))