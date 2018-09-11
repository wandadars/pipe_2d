import sys
import numpy as np


class PipeGrid(object):
    def __init__(self, length, height, width, dh_wall, dh_centerline, num_i, num_j, num_k):
        self.length = length
        self.height = height
        self.width = width
        self.dh_wall = dh_wall
        self.dh_centerline = dh_centerline
        self.num_i = num_i
        self.num_j = num_j
        self.num_k = num_k

    def compute_delta(self):
        tolerance = 1.0e-10
        
        b = self.compute_b_parameter()
        converged = False
        delta = np.sqrt(6.0 * (b - 1))  #initial guess
        while converged == False:
            numerator = -(np.sinh(delta) / delta - b)
            denominator = -np.sinh(delta) / (delta**2) + np.cosh(delta) / delta
            d_delta = numerator / denominator 
            delta += d_delta
            if abs(d_delta) / delta < tolerance:
                converged = True
                print('Converged value of delta: {0:10.6e}'.format(delta))

        return delta

    def compute_b_parameter(self):
        return 1.0 / ((self.num_j - 1)//2 * np.sqrt(self.dh_wall * self.dh_centerline))
    
    def compute_y_distribution(self):
        delta = self.compute_delta() 
        a = np.sqrt(self.dh_centerline / self.dh_wall)

        s = np.zeros(self.num_j)
        s[0] = 0.0
        for j in range(1, (self.num_j + 1)//2):
            numerator = np.tanh(delta*(j / ((self.num_j-1)/2.0) - 0.5))
            denominator = np.tanh(delta / 2.0)
            u = 0.5 * ( 1.0 + numerator / denominator) 

            value = 0.5 * self.height * u /(a + u*(1.0 - a))
            s[j] = value 
            print('j, s: ' + str(j) + ', ' + str(s[j]))

        for j in range(self.num_j - 1, (self.num_j -1)//2, -1):
            s[j] = self.height - s[num_j-1-j]
            print('j, s: ' + str(j) + ', ' + str(s[j]))

        return s

    def write_grid(self):
        coord_format = '{0:<17.10e}\n'
        f = open('./file.grd', 'w+')
        f.write('1\n')
        f.write('{0:d} {1:d} {2:d}\n'.format(self.num_i, self.num_j, self.num_k))
        for k in range(self.num_k):
            for j in range(self.num_j):
                for i in range(self.num_i):
                    x = self.length * (float(i)/float(self.num_i-1))
                    f.write(coord_format.format(x))

        s = self.compute_y_distribution()
        for k in range(self.num_k):
            for j in range(self.num_j):
                for i in range(self.num_i):
                    y = s[j]
                    f.write(coord_format.format(y))
        
        for k in range(self.num_k):
            for j in range(self.num_j):
                for i in range(self.num_i):
                    z = self.width * float(k) / float(self.num_k-1)
                    f.write(coord_format.format(z))

        f.close()


class HalfPipeGrid(PipeGrid):
    def __init__(self, length, height, width, dh_wall, dh_centerline, num_i, num_j, num_k):
        super(HalfPipeGrid, self).__init__(length, height, width, dh_wall, dh_centerline, num_i, num_j, num_k)

    def compute_b_parameter(self):
        return 1.0 / ((self.num_j - 1) * np.sqrt(self.dh_wall * self.dh_centerline))
    
    def compute_y_distribution(self):
        """
        Compute the y distribution using Thompson's method detailed on page 306 
        of his grid generation book.
        """
        delta = self.compute_delta() 
        a = np.sqrt(self.dh_centerline / self.dh_wall)

        s = np.zeros(self.num_j)
        s[0] = 0.0
        for j in range(1, self.num_j ):
            numerator = np.tanh(delta*(j / (self.num_j-1) - 0.5))
            denominator = np.tanh(delta / 2.0)
            u = 0.5 * ( 1.0 + numerator / denominator) 
            value = self.height * u /(a + u*(1.0 - a))
            s[j] = value 
            print('j, s: ' + str(j) + ', ' + str(s[j]))

        return s


if len(sys.argv) < 9:
    print('Usage: pipe-2D nI nJ L H dhAxis dhCenterline half/full > file.grd')
    sys.exit()


length = float(sys.argv[3])
height = float(sys.argv[4])
width = float(sys.argv[5])
dh_wall = float(sys.argv[6])
dh_centerline = float(sys.argv[7])
num_i = int(sys.argv[1])
num_j = int(sys.argv[2])
num_k = 2

grid_type = sys.argv[8]

if num_j  % 2 < 1:
    print('Number of points in J direction(nJ) must be odd')
    sys.exit()

if grid_type.strip() == 'half':
    grid = HalfPipeGrid(length, height, width, dh_wall, dh_centerline, num_i, num_j, num_k)
elif grid_type.strip() == 'full':
    grid = PipeGrid(length, height, width, dh_wall, dh_centerline, num_i, num_j, num_k)

grid.write_grid()

