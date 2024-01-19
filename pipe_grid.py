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

        self.parse_input()

    def parse_input(self):
        if len(sys.argv) < 10:
            self.print_usage()
            sys.exit()

        self.num_i = int(sys.argv[1])
        self.num_j = int(sys.argv[2])
        self.num_k = int(sys.argv[3])
        self.length = float(sys.argv[4])
        self.y_max = float(sys.argv[5])
        self.width = float(sys.argv[6])
        self.dh_wall = float(sys.argv[7])
        self.dh_centerline = float(sys.argv[8])
        
        self.y_min = 0.0 #This value is hardcoded for simplicity, but can be added as an argument if needed
      
    def print_usage(self):
        print('Usage: pipe_grid nX nY nZ Lx Ly Lz dhWall dhCenterline half')

    def compute_delta(self, wall_spacing, centerline_spacing, num_points, tolerance=1.0e-10):
        """
        Computes the delta value used in the generation of the y-coordinate grid spacing.

        This method iteratively calculates the delta value, which is used to determine the
        non-linear distribution of grid points in the y-direction. The calculation is based
        on a specified tolerance and continues until convergence.

        Returns:
            float: The converged delta value used for grid spacing calculations.
            float: The value of the centerline_spacing input. This value can be updated in this method.

        Raises:
            ValueError: If input parameters lead to an invalid computation.
        """

        # Validate inputs
        if num_points <= 1 or wall_spacing <= 0 or centerline_spacing <= 0:
            raise ValueError(
                "Invalid input parameters. Ensure num_j > 1, dh_wall > 0, and dh_centerline > 0.")

        converged = False
        while converged == False:
            b = 1.0 / ((num_points - 1) * np.sqrt(wall_spacing * centerline_spacing))
            # Check for valid square root computation
            if b < 1:
                print('Invalid combination of num_j, dh_wall, and dh_centerline leading to negative square root. Adjusting the spacing of the final cell to resolve the issue.')
                self.dh_centerline = 0.5 * centerline_spacing
                centerline_spacing *= 0.5
                print('New value for dh_centerline is: ' + str(centerline_spacing))
                b = 1.0 / ((num_points - 1) * np.sqrt(wall_spacing * centerline_spacing))
            else:
                converged = True

        converged = False
        delta = np.sqrt(6.0 * (b - 1))  # initial guess
        print('Initial guess for Delta: ' + str(delta))
        while converged == False:
            numerator = -(np.sinh(delta) / delta - b)
            denominator = -np.sinh(delta) / (delta**2) + np.cosh(delta) / delta
            d_delta = numerator / denominator
            delta += d_delta
            if abs(d_delta) / delta < tolerance:
                converged = True
                print('Converged value of delta: {0:10.6e}'.format(delta))

        return delta, centerline_spacing

    def compute_arclength_coordinates(self, wall_spacing, centerline_spacing, num_points):
        """
        Computes the y-coordinate distribution using Thompson's method detailed on page 306
        of his grid generation book.

        This method applies Thompson's spacing method to compute the distribution of grid points
        along the y-axis. The method involves calculating a scaling factor 's' for each grid point
        based on the computed delta value and the ratio of the height at the wall to the height at
        the centerline.

        Returns:
            numpy.ndarray: An array of scaled values 's' for each y-coordinate grid point.
        """

        delta, centerline_spacing = self.compute_delta(wall_spacing, centerline_spacing, num_points)

        a = np.sqrt(centerline_spacing / wall_spacing)
        s = []
        for j in range(0, self.num_j):
            epsilon = j
            I = num_points - 1
            numerator = np.tanh(delta * ((epsilon / I) - 0.5))
            denominator = np.tanh(delta / 2.0)
            u = 0.5 * (1.0 + numerator / denominator)

            s.append(u / (a + u * (1.0 - a)))
            print('j, s: {0:d}, {1:<19.10e}'.format(j, s[-1]))
        return s

    def compute_y_distribution(self):
        """
        Computes the actual y-coordinate values for the grid points.

        This method uses the scaled values obtained from 'compute_arclength_coordinates' to
        calculate the actual y-coordinates of the grid points. The y-coordinates are interpolated
        between the minimum (y_min) and maximum (y_max) values of the grid.

        Returns:
            y_points: A list of computed y-coordinate values for the grid.
        """
        arclength_points = self.compute_arclength_coordinates(self.dh_wall, self.dh_centerline, self.num_j)
        y_points = []
        for i, s in enumerate(arclength_points):
            y_points.append(self.y_min + (self.y_max - self.y_min) * s)

        return y_points

    def compute_x_distribution(self):
        x = []
        for i in range(self.num_i):
            x.append(self.length * (float(i) / float(self.num_i - 1)))
        return x
        
    def compute_z_distribution(self):
        z = []
        for i in range(self.num_k):
            z.append(self.width * (float(i) / float(self.num_k - 1)))
        return z
    
    def compute_coordinates(self):
        """
        Computes and stores the x, y, and z coordinates for all grid points.

        This method precomputes the grid coordinates for each point in the grid and stores them
        in a structured format for easy access. The coordinates are stored in a 3D matrix
        where each element is a tuple (x, y, z) representing the coordinates of a grid point.
        """
        # Precompute coordinate distributions
        x_distribution = self.compute_x_distribution()
        y_distribution = self.compute_y_distribution()
        z_distribution = self.compute_z_distribution()
        
        # Initialize a 3D matrix to store coordinates
        coordinates = np.empty((self.num_i, self.num_j, self.num_k), dtype=tuple)

        # Compute and store coordinates
        for i in range(self.num_i):
            for j in range(self.num_j):
                for k in range(self.num_k):
                    x = x_distribution[i]
                    y = y_distribution[j]
                    z = z_distribution[k]
                    coordinates[i, j, k] = (x, y, z)
        return coordinates

    def write_grid(self):
        """
        Writes the computed grid to a file.

        This method outputs the precomputed grid coordinates to a file in a specified format.
        Each line in the file corresponds to the x, y, and z coordinates of a grid point.
        """
        coord_format = '{0:<17.10e}\n'
        coordinates = self.compute_coordinates()
        f = open('./file.grd', 'w+')
        f.write('1\n')
        f.write(
            '{0:d} {1:d} {2:d}\n'.format(
                self.num_i,
                self.num_j,
                self.num_k))
        for k in range(self.num_k):
            for j in range(self.num_j):
                for i in range(self.num_i):
                    x = coordinates[i, j, k][0]
                    f.write(coord_format.format(x))

        for k in range(self.num_k):
            for j in range(self.num_j):
                for i in range(self.num_i):
                    y = coordinates[i, j, k][1]
                    f.write(coord_format.format(y))

        for k in range(self.num_k):
            for j in range(self.num_j):
                for i in range(self.num_i):
                    z = coordinates[i, j, k][2]
                    f.write(coord_format.format(z))
        f.close()


class FullPipeGrid(HalfPipeGrid):
    def __init__(self):
        super().__init__()
       
    def print_usage(self):
        print('Usage: pipe_grid nX nY nZ Lx Ly Lz dhWall dhCenterline full')

    def compute_y_distribution(self):
        """
        Computes the actual y-coordinate values for the grid points for a full pipe grid.

        This method extends the scaled values obtained from 'compute_arclength_coordinates'
        to create a symmetrical distribution of y-coordinates for a full pipe. It mirrors
        the y-coordinates about the centerline.

        Returns:
            y_points: A list of computed y-coordinate values for the full grid.
        """
        half_pipe_y_points = super().compute_y_distribution()
        
        # Adjust the number of y points for full pipe
        self.num_j = 2 * (self.num_j - 1) + 1

        # Mirror the y-coordinates, scaling non-wall points by 0.5
        full_pipe_y_points = [self.y_min]  # Start with the first point at original scale
        for y in half_pipe_y_points[1:]:
            # Scale down non-wall points by 0.5
            scaled_y = self.y_min + 0.5 * (y - self.y_min)
            full_pipe_y_points.append(scaled_y)

        for y in reversed(half_pipe_y_points[0:-1]):
            # Mirror the scaled points and adjust to top half of the pipe
            mirrored_y = self.y_max - 0.5 * (y - self.y_min)
            full_pipe_y_points.append(mirrored_y)

        return full_pipe_y_points


class PipeGrid3D(HalfPipeGrid):
    def __init__(self):
        self.length = None
        self.y_min = None
        self.y_max = None
        self.dh_wall = None
        self.dh_centerline = None
        self.num_i = None  # Axial distribution
        self.num_j = None  # Radial distribution
        self.num_k = None  # Theta  distribution

        self.parse_input()

    def parse_input(self):
        if len(sys.argv) < 10:
            self.print_usage()
            sys.exit()

        self.num_i = int(sys.argv[1])
        self.num_j = int(sys.argv[2])
        self.num_k = int(sys.argv[3])
        self.y_min = float(sys.argv[4])
        self.y_max = float(sys.argv[5])
        self.length = float(sys.argv[6])
        self.dh_wall = float(sys.argv[7])
        self.dh_centerline = float(sys.argv[8])
        
        print('NI (Axial): ' + str(self.num_i))
        print('NJ (Radial): ' + str(self.num_j))
        print('NK (Theta): ' + str(self.num_k))
        print('R_inner: ' + str(self.y_min))
        print('R_outer: ' + str(self.y_max))
        print('Length: ' + str(self.length))
        print('Wall Spacing: ' + str(self.dh_wall))
        print('Center Spacing: ' + str(self.dh_centerline))       

    def print_usage(self):
        print('Usage: pipe_grid nX nY nZ R_inner R_outer Lx dhWall dhCenterline 3d')

    def compute_coordinates(self):
        """
        Computes and stores the x, y, and z coordinates for all grid points in a 3D rotated pipe grid.

        This method overrides the compute_coordinates method to account for the 3D structure of the pipe.
        """
        # Precompute radial (y) distribution and x distribution
        radial_distribution = self.compute_y_distribution()
        x_distribution = self.compute_x_distribution()

        # Initialize a 3D matrix to store coordinates
        coordinates = np.empty((self.num_i, self.num_j, self.num_k), dtype=tuple)

        # Compute and store coordinates
        for i in range(self.num_i):
            x = x_distribution[i]  # Axial distribution
            for j in range(self.num_j):
                r = radial_distribution[j]  # Radial distribution
                for k in range(self.num_k):
                    theta = 2.0 * np.pi * (float(k) / float(self.num_k - 1))  # Theta distribution
                    y = r * np.cos(theta)
                    z = r * np.sin(theta)
                    coordinates[i, j, k] = (x, y, z)

        return coordinates
        

def create_grid(grid_type):
    if grid_type == 'half':
        return HalfPipeGrid()
    elif grid_type == 'full':
        return FullPipeGrid()
    elif grid_type == '3d':
        return PipeGrid3D()
    else:
        raise ValueError(f"Invalid grid type: {grid_type}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py [arguments] [grid_type]")
        sys.exit(1)
        
    grid_type = sys.argv[-1]
    try:
        grid = create_grid(grid_type)
        grid.write_grid()
    except ValueError as e:
        print(e)
        sys.exit(1)

if __name__ == "__main__":
    main()
