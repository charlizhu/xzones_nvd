import numpy as np

class warning(object):

    def __init__(self, predicted_path, current_location):
        self.pp = predicted_path
        self.cl = current_location
        self.rear_bounds = [-2.6,0]
        self.front_bounds = [0,8.5]
        self.side_bounds = [0,0.6]
        self.rear_sides = [-9,0]
        self.scan_data()

    def scan_data(self):
        if (self.cl[1] < self.side_bounds[1] and self.cl[1 > self.side_bounds[0]]) \
                and (self.cl[0] > self.rear_sides[0] and self.cl[0] < self.rear_sides[1]):
            # print("WARNING: Too Close Laterally")
            return "lateral"
        else:
            # I was thinking if using any() would be less memory-intensive but I couldn't figure that out...
            for point in self.pp:
                if (point[0] > self.rear_bounds[0]) and (point[0] < self.rear_bounds[1]):
                    if point[1] < self.side_bounds[1]:
                        # print("WARNING: Rear Hook")
                        return "rear"
                if (point[0] > self.front_bounds[0]) and (point[0] < self.front_bounds[1]):
                    if point[1] < self.side_bounds[1]:
                        # print("WARNING: Front Hook")
                        return "front"
            # print("SAFE")
            return "safe"



# MAIN is just used for testing right now. Get rid if no longer needed.
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
    myFile = open("sim_data\one_car_safe.txt",'r') # change as needed
    for points in myFile:
        current_line = points.split()
        predicted_path.append([float(x) for x in current_line])
        if count >= numPoints:
            current_location = predicted_path[0]
            x = warning(predicted_path,current_location)
            print(x.scan_data())
            predicted_path.pop(0)
        count = count + 1
        memory_size, memory_peak = tracemalloc.get_traced_memory()
        run_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        print("Memory Usage (MB): " + str(memory_peak/(10**6)) + " Run Time (sec): " + str(run_time))