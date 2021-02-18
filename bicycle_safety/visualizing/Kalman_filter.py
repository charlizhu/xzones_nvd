import numpy as np


class Kalman(object):
    def __init__(self, initial_xx: float,
                 initial_xy: float,
                 initial_vx: float,
                 initial_vy: float,
                 accel_variance: float):
        # mean of state GRV
        self.DIM = 4
        self.offsetX = 0
        self.offsetV = 2
        self._x = np.zeros(self.DIM)

        self._x[self.offsetX] = initial_xx
        self._x[self.offsetX + 1] = initial_xy
        self._x[self.offsetV] = initial_vx
        self._x[self.offsetV + 1] = initial_vy

        self._accel_variance = accel_variance

        # covariance of state GRV
        self._P = np.eye(self.DIM)
        # A and B and others for predict ahead
        self._DT = 0
        self._dt = 0
        self._accel = np.zeros(self.DIM)
        self._A = np.eye(self.DIM)
        self._B = np.eye(self.DIM)

    def predict_ahead(self, state, A, B, accel, DT, dt):
        time_step = dt
        predicted_x = []
        for i in range(int(DT/dt)):
            A[self.offsetX, self.offsetV] = time_step
            A[self.offsetX + 1, self.offsetV + 1] = time_step

            B[self.offsetX][self.offsetX] = 0.5 * time_step ** 2
            B[self.offsetX + 1][self.offsetX + 1] = 0.5 * time_step ** 2
            B[self.offsetV][self.offsetV] = time_step
            B[self.offsetV + 1][self.offsetV + 1] = time_step
            result = A.dot(state) + B.dot(accel)
            predicted_x.append(result.tolist())
            time_step += dt
        return predicted_x

    def predict(self, dt: float, interpDT: float, accel):
        # x = A x + B*u(t) + wt
        # P = A P A^T + Q (ie. G Gt a)
        #accel = [ax, ay, ax, ay]
        A = np.eye(self.DIM)
        A[self.offsetX, self.offsetV] = dt
        A[self.offsetX+1, self.offsetV+1] = dt


        G = np.zeros((self.DIM, 1))
        G[self.offsetX] = 0.5 * dt ** 2
        G[self.offsetX+1] = 0.5 * dt ** 2
        G[self.offsetV] = dt
        G[self.offsetV+1] = dt

        B = np.eye(self.DIM)
        B[self.offsetX][self.offsetX] = 0.5 * dt ** 2
        B[self.offsetX + 1][self.offsetX + 1] = 0.5 * dt ** 2
        B[self.offsetV][self.offsetV] = dt
        B[self.offsetV + 1][self.offsetV + 1] = dt

        new_x = A.dot(self._x) + B.dot(accel)
        new_P = A.dot(self._P).dot(A.T) + G.dot(G.T) * self._accel_variance ** 2

        self._P = new_P
        self._x = new_x
        #prediction = self.predict_ahead(self._x, A, B, accel, interpDT - dt, dt)
        self._DT = interpDT
        self._dt = dt
        self._accel = np.copy(accel)
        self._A = np.copy(A)
        self._B = np.copy(B)

        #return prediction

    def update(self, meas_value, meas_variance, measurements):
        # y = z - H x
        # S = H P Ht + R
        # K = P Ht S^-1
        # x = x + K y
        # P = (I - K H) * P

        H = np.eye(self.DIM)
        z = meas_value
        '''
        R = np.matrix([[np.std(measurements[0, (n - i):n]) ** 2, 0.0],
                       [0.0, np.std(measurements[1, (n - i):n]) ** 2]])
        '''
        '''
        R = np.eye(self.DIM)
        R[self.offsetX][self.offsetX] = meas_variance[self.offsetX] ** 2
        R[self.offsetX + 1][self.offsetX + 1] = meas_variance[self.offsetX + 1] ** 2
        R[self.offsetV][self.offsetV] = meas_variance[self.offsetV] ** 2
        R[self.offsetV + 1][self.offsetV + 1] = meas_variance[self.offsetV + 1] ** 2
        '''
        R = np.array([[np.std([i[0] for i in measurements]) ** 2, 0.0, 0.0, 0.0],
                       [0.0, np.std([i[1] for i in measurements]) ** 2, 0.0, 0.0],
                       [0.0, 0.0, np.std([i[2] for i in measurements]) ** 2, 0.0],
                       [0.0, 0.0, 0.0, np.std([i[3] for i in measurements]) ** 2]])
        y = z - H.dot(self._x)
        S = H.dot(self._P).dot(H.T) + R
        try:
            K = self._P.dot(H.T).dot(np.linalg.inv(S))
        except:
            #numpy.linalg.pinv
            K = self._P.dot(H.T).dot(np.linalg.pinv(S))
            print("pseudo inverse")


        new_x = self._x + K.dot(y)
        new_P = (np.eye(4) - K.dot(H)).dot(self._P)

        self._P = new_P
        self._x = new_x
        prediction = self.predict_ahead(self._x, self._A, self._B, self._accel, self._DT - self._dt, self._dt)
        return prediction
