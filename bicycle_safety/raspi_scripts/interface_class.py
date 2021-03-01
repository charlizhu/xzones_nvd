#!/usr/bin/python

class interface(object):
    def __init__(self, port_user, port_data, baud_usr, baud_data):
        self._PortUsr = port_user
        self._PortData = port_data
        self._BaudUsr = baud_usr
        self._BaudData = baud_data
        self.user_ser = serial.Serial(self._PortUsr, self._BaudUsr, timeout=2)
        self.data_ser = serial.Serial(self._PortData, self._BaudData, timeout=1)
    def load_config(self):
        cliCFG = ['advFrameCfg', 'sensorStart']
        ticks = 0
        '''
        while 1:
            ticks += 1
            self.user_ser.write('x'.encode())
            cc = self.user_ser.read(100)
            print(cc)
            cc = cc.decode("utf-8")
            cc = cc.replace(chr(10), '')
            cc = cc.replace(chr(13), '')
            
            if not (not (cc)):
                break
            elif ticks == 10:
                print("cannot connect to the radar device")
                return -1
            time.sleep(0.1)
        '''
        print('Sending configuration to Radar sensor')
        for i in range(0, len(cliCFG)):
            print('Command: ' + cliCFG[i])
            self.user_ser.write(cliCFG[i].encode() + '\n'.encode())
            cc = self.user_ser.readline()
            print(cc)
            if (cc.decode("utf-8").find('Done') == 1):
                print(cc)
                break
            elif (cc.decode("utf-8")).find('not recognized as a CLI command') != -1:
                print(cc)
                return -2
            time.sleep(0.1)
        return 1
    def data_readline(self):
        return self.data_ser.readline()
    def data_read(self, bytes):
        return self.data_ser.read(bytes)


if __name__ == "__main__":
    import time
    import serial
    import tracemalloc

    interface_pi = interface('/dev/ttyACM0','/dev/ttyACM1', 115200, 921600)
    interface_pi.load_config()
    time.sleep(2)
    
    tracemalloc.start()
    start_time = time.time()
    myfile = open('data8.txt','a')

    count = 0
    #mydata = serial.Serial('/dev/ttyACM1', 921600)
    while True:
        print(count)
        if interface_pi.data_ser.in_waiting:
            x = interface_pi.data_read(interface_pi.data_ser.in_waiting)
            myfile.write(str(x))
            count = count + 1
            current,peak = tracemalloc.get_traced_memory()
            print(current/10**6,peak/10**6)
            current_time = time.time()
            print(current_time-start_time)
            start_time=current_time
