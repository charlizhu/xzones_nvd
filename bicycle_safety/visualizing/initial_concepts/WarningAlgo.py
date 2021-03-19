class warning_algo(object):

    def __init__(self, predicted_path):
        self.pp = predicted_path
        self.cl = predicted_path[0] # The first item in predicted_path is the current_location (i.e.: cl).
        self.rear_bounds = [-2.6,0]
        self.front_bounds = [0,8.5]
        self.side_bounds = [0,0.6]
        self.rear_sides = [-9,0]
        self.turning_tolerance = 0.035
        self.scan_data()

    def scan_data(self):
        # First case is where the car is approaching too close to the side.
        # This is just to show that the algorithm can have other scenarios too
        # such as the rear-end head-on collision, for after the MVP.

        if (self.side_bounds[0] < self.cl[1] < self.side_bounds[1]) \
                and (self.rear_sides[0] < self.cl[0] < self.rear_sides[1]) \
                and ((self.pp[-1][-1] - self.cl[1]) < self.turning_tolerance):
            return "lateral"
        else:
            count = 0
            # Now it checks for side hook conditions.
            for point in self.pp:
                # Rear side hook.
                if (self.rear_bounds[0] < point[0] < self.rear_bounds[1]):
                    if (point[1] < self.side_bounds[1]):
                        print(self.cl, point, count)
                        return "rear"
                # Front side hook.
                if (self.front_bounds[0] < point[0] < self.front_bounds[1]):
                    if (point[1] < self.side_bounds[1]):
                        print(self.cl, point, count)
                        return "front"
                count = count + 1
            # If no warnings were previously raised, then the cyclist is safe.
            return "safe"

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
            x = warning_algo(predicted_path)
            print(x.scan_data())
            predicted_path.pop(0)
        count = count + 1
        memory_size, memory_peak = tracemalloc.get_traced_memory()
        run_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        print("Memory Usage (MB): " + str(memory_peak/(10**6)) + " Run Time (sec): " + str(run_time))