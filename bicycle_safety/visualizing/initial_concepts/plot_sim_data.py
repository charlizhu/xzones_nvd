import matplotlib.pyplot as plt

myFile = open("sim_data\one_car_rear_sidehook.txt",'r') # change as needed
for points in myFile:
    current_line = points.split()
    plt.plot(float(current_line[0]),float(current_line[1]),'r.')

plt.ylim([-3,2])

plt.grid('on')
plt.show()