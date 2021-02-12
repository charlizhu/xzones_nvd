import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline
from sklearn.linear_model import LinearRegression
import tracemalloc

def countObjects():
    f = open("data3_2.txt", "r")

    allObjects = []
    everyID = []
    record = []
    count = 0

    for x in f:
        count = count + 1
        temp = x.split()
        # print(temp)
        if temp[0] not in allObjects:
            allObjects.append(temp[0])
        everyID.append(temp[0])

    print(count)

    print(len(allObjects))

    for x in allObjects:
        record.append([x,everyID.count(x)])

    record.sort(key=lambda x: x[1])
    print((record))

def traceObjs():
    f = open("data3_2.txt", "r")

    obj2x = []
    obj2y = []
    obj11x = []
    obj11y = []
    obj3x = []
    obj3y = []
    obj7x = []
    obj7y = []
    otherx = []
    othery = []
    for x in f:
        temp = x.split()

        if temp[0] == '2':
            obj2x.append(float(temp[1]))
            obj2y.append(float(temp[2]))

        elif temp[0] == '3':
            obj3x.append(float(temp[1]))
            obj3y.append(float(temp[2]))

        elif temp[0] == '11':
            obj11x.append(float(temp[1]))
            obj11y.append(float(temp[2]))

        elif temp[0] == '7':
            obj7x.append(float(temp[1]))
            obj7y.append(float(temp[2]))

        else:
            otherx.append(float(temp[1]))
            othery.append(float(temp[2]))
    #
    plt.plot(obj2x,obj2y,'r-')
    # plt.plot(obj3x,obj3y,'g-')
    plt.plot(obj11x, obj11y, 'b-')
    plt.plot(obj7x, obj7y, 'm-')
    # plt.plot(obj1x, obj1y, 'o-')
    # plt.plot(obj12x, obj12y, 'y-')
    # plt.plot(otherx,othery,'k.',alpha=0.05)
    # plt.axis('equal')
    plt.grid('on')
    plt.show()

def correctness2():
    f = open("data5_2.txt", "r")

    count = 0

    car2 = []
    car11 = []
    car7 = []
    numObjs = []
    countarray = []
    confused = []

    for x in f:
        print(x)
        temp = x.split()
        print(temp)
        confused.append(temp)

    for temp in confused:

        numObjs.append(len(temp))

        if '2' in temp:
            car2.append(2)
        elif '2' not in temp:
            car2.append(1)


        if '11' in temp:
            car11.append(0)
        elif '11' not in temp:
            car11.append(-1)

        if '7' in temp:
            car7.append(-2)
        elif '7' not in temp:
            car7.append(-3)
        countarray.append(count)
        count = count + 1

    # print(len(car6))
    # print(len(car13))
    # print(len(car17))
    # print(len(numObjs))
    # print(len(countarray))

    f, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(countarray, car2,'r')
    ax1.plot(countarray, car11,'b')
    ax1.plot(countarray, car7,'m')
    ax1.set_title('Locking on to "maybe cars" objects')
    ax2.plot(countarray,numObjs)
    ax2.set_title("Number of objects ID'd in total")
    ax1.axes.get_xaxis().set_visible(False)
    ax2.set_xlabel("Data point number")
    ax1.set_yticks([-3,-2,-1,0,1,2])
    ax1.set_yticklabels((['Not locked on', 'Locked on', 'Not locked on', 'Locked on', 'Not locked on', 'Locked on']))

    plt.show()

def correctness():
    f = open("connect2.txt", "r")

    count = 0

    car = []
    notcar = []
    noise = []
    countarray = []
    confused = []

    for x in f:
        print(x)
        temp = x.split()
        print(temp)
        confused.append(temp)

    print(confused)

    for temp in confused:

        if '27' in temp:
            car.append(2)
            noise.append(len(temp)-1)
        elif '27' not in temp:
            car.append(1)
            noise.append(len(temp))

        if '0' in temp:
            notcar.append(-1)
        elif '0' not in temp:
            notcar.append(-2)
        countarray.append(count)
        count = count + 1

    print(len(countarray))
    print(len(car))
    print(len(notcar))
    print(len(noise))

    f, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(countarray, car)
    ax1.plot(countarray, notcar)
    ax1.set_title('Locking on to "main" objects')
    ax2.plot(countarray,noise)
    ax2.set_title('Amount of noise detected (includes non-car)')
    ax1.axes.get_xaxis().set_visible(False)
    ax2.set_xlabel("Data point number")
    ax1.set_yticks([-2,-1,0,1,2])
    ax1.set_yticklabels((['Not locked on to non-car', 'Locked on to non-car', ' ', 'Not locked on to car', 'Locked on to car']))

    plt.show()
# countObjects()
traceObjs()
# correctness2()
