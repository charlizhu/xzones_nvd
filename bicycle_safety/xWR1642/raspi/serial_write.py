import time
import serial
import tracemalloc
# import int

tracemalloc.start()
start_time = time.time()
myfile = open('data8.txt','a')

count = 0

while True:
    print(count)
    mydata = serial.Serial(port='/dev/ttyACM1')
    x = mydata.readline()
    myfile.write(str(x) + '\n')
    count = count + 1
    current,peak = tracemalloc.get_traced_memory()
    print(current/10**6,peak/10**6)
    current_time = time.time()
    print(current_time-start_time)
    start_time=current_time
    if count == 10:
        myfile.close()
        break
    continue
#     y = int.from_bytes(x,'little',signed=True)
#     print(y)
