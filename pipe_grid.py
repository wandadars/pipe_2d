import sys
import numpy as np


class HalfPipeGrid(object):
    def __init__(self):
        self.length = None
        self.height = None 
        self.width = None 
        self.dh_wall = None 
        self.dh_centerline = None
        self.num_i = None 
        self.num_j = None 
        self.num_k = None 
        
        self.y_points = None

        self.parse_input()

    def parse_input(self):
        if len(sys.argv) < 10:
            self.print_usage()
            sys.exit()

        self.length = float(sys.argv[4])
        
        self.y_min = 0.0
        self.y_max = float(sys.argv[5])
        self.width = float(sys.argv[6])
        self.dh_wall = float(sys.argv[7])
        self.dh_centerline = float(sys.argv[8])
        self.num_i = int(sys.argv[1])
        self.num_j = int(sys.argv[2])
        self.num_k = int(sys.argv[3])

    def print_usage(self):
        print('Usage: pipe_grid nX nY nZ L H W dhWall dhCenterline half')

    def compute_delta(self):
        tolerance = 1.0e-10
        
        b = 1.0 / ((self.num_j - 1) * np.sqrt(self.dh_wall * self.dh_centerline)) 
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

    def compute_arclength_coordinates(self):
        """
        Compute the y distribution using Thompson's method detailed on page 306 
        of his grid generation book.
        """
        a = np.sqrt(self.dh_centerline / self.dh_wall)
        delta = self.compute_delta() 
        
        s = np.zeros(self.num_j)
        s[0] = 0.0 
        for j in range(1, len(s) ):
            epsilon = j
            I = self.num_j - 1
            numerator = np.tanh(delta*((epsilon / I) - 0.5))
            denominator = np.tanh(delta / 2.0)
            u = 0.5 * ( 1.0 + numerator / denominator) 
            s[j] = u / (a + u * (1.0 - a)) 
            print('j, s: {0:d}, {1:<19.10e}'.format(j, s[j]))
        return s

    def compute_y_distribution(self):
        arclength_points = self.compute_arclength_coordinates()
        self.y_points = []
        for i, s in enumerate(arclength_points):
            self.y_points.append(self.y_min + (self.y_max - self.y_min)*s)

    def write_grid(self):
        coord_format = '{0:<17.10e}\n'
        self.compute_y_distribution()
        f = open('./file.grd', 'w+')
        f.write('1\n')
        f.write('{0:d} {1:d} {2:d}\n'.format(self.num_i, len(self.y_points), self.num_k))
        for k in range(self.num_k):
            for j in range(len(self.y_points)):
                for i in range(self.num_i):
                    x = self.compute_x_coordinate(i, j, k)
                    f.write(coord_format.format(x))

        for k in range(self.num_k):
            for j in range(len(self.y_points)):
                for i in range(self.num_i):
                    y = self.compute_y_coordinate(i, j, k)
                    f.write(coord_format.format(y))
        
        for k in range(self.num_k):
            for j in range(len(self.y_points)):
                for i in range(self.num_i):
                    z = self.compute_z_coordinate(i, j, k)
                    f.write(coord_format.format(z))
        f.close()

    def compute_x_coordinate(self, i, j, k):
        x = self.length * (float(i)/float(self.num_i-1))
        return x

    def compute_y_coordinate(self, i, j, k):
        return self.y_points[j]

    def compute_z_coordinate(self, i, j, k):
        z = self.width * (float(k)/float(self.num_k-1))
        return z


class FullPipeGrid(HalfPipeGrid):
    def __init__(self):
        super(FullPipeGrid, self).__init__()

    def print_usage(self):
        print('Usage: pipe_grid nX nY nZ L H W dhWall dhCenterline full')

    def compute_y_distribution(self):
        arclength_points = self.compute_arclength_coordinates()
        self.y_points = []
        for s in arclength_points:
            y_coord = self.y_min + 0.5*(self.y_max - self.y_min)*s
            self.y_points.append(y_coord)

        for s in reversed(arclength_points[:-1]):
            y_coord = self.y_max - (self.y_min + 0.5*(self.y_max - self.y_min)*s)
            self.y_points.append(y_coord)


class PipeGrid3D(HalfPipeGrid):
    def __init__(self):
        self.length = None
        self.r_inner = None 
        self.r_outer = None 
        self.dh_wall = None 
        self.dh_centerline = None
        self.num_i = None  #Axial distribution
        self.num_j = None  #Theta distribution
        self.num_k = None  #Radial distribution
    
        self.y_points = None

        self.parse_input()
        self.print_inputs()

    def parse_input(self):
        if len(sys.argv) < 10:
            self.print_usage()
            sys.exit()

        self.r_inner = float(sys.argv[4])
        self.r_outer = float(sys.argv[5])
        self.length = float(sys.argv[6])
        self.dh_wall = float(sys.argv[7])
        self.dh_centerline = float(sys.argv[8])
        self.num_i = int(sys.argv[1])
        self.num_j = int(sys.argv[2])
        self.num_k = int(sys.argv[3])

    def print_inputs(self):
        print('NI: ' + str(self.num_i))
        print('NJ: ' + str(self.num_j))
        print('NK: ' + str(self.num_k))
        print('R_inner: ' + str(self.r_inner))
        print('R_outer: ' + str(self.r_outer)) 
        print('Length: ' + str(self.length))
        print('Wall Spacing: ' + str(self.dh_wall))
        print('Center Spacing: ' + str(self.dh_centerline)) 

    def print_usage(self):
        print('Usage: pipe_grid nX nY nZ R_inner R_outer L dhWall dhCenterline 3d')

    def compute_y_distribution(self):
        arclength_points = self.compute_arclength_coordinates()
        self.y_points = []
        for i, s in enumerate(arclength_points):
            self.y_points.append(self.r_inner + (self.r_outer - self.r_inner)*s)

    def compute_x_coordinate(self, i, j, k):
        x = self.length * (float(i)/float(self.num_i-1))
        return x

    def compute_y_coordinate(self, i, j, k):
        r = self.y_points[j]
        theta = 2.0 * np.pi * (float(k) / float(self.num_k-1))
        y = r*np.cos(theta)
        return y 

    def compute_z_coordinate(self, i, j, k):
        r = self.y_points[j]
        theta = 2.0 * np.pi * (float(k) / float(self.num_k-1))
        z = r*np.sin(theta)
        return z


#Main
grid_type = sys.argv[-1]
print('Detected Grid Type: ' + grid_type)
if grid_type.strip() == 'half':
    grid = HalfPipeGrid()
elif grid_type.strip() == 'full':
    grid = FullPipeGrid()
elif grid_type.strip() == '3d':
    grid = PipeGrid3D()
else:
    print('No Grid Type detected')
    sys.exit()

grid.write_grid()


