This utility is a python script that generates simple meshes in the Plot3D format.

Sample usage to generate a grid that has a boundary layer at y = 0:
python3 pipe_grid.py 9 9 2 1 1 5.0e-4 3.0e-6 0.1 half 

grid has:
7 points in the x direction
9 points in the y direction(boundary layer direction)
2 points in the z direction(single cell thick)
A length of 1
A height of 1
A width of 5.0e-4
An initial cell spacing of 3.0e-6
A final cell spacing(at y=1) of 0.1
"half" signifies that there is no boundary layer at y=1


Sample usage to generate a grid that has a boundary layer at y = 0 and y = 1:
python3 pipe_grid.py 11 11 11 4.0 1.0 0.01 1.0e-3 0.1 full 

grid has:
11 points in the x direction
21 points in the y direction(boundary layer direction), 2*11 -1(for overlapped middle point) 
11 points in the z direction
A length of 4
A width of 0.01(z direction)
A height of 1
An initial cell spacing of 1.0e-3
A final cell spacing(at y=0.5*H) of 0.1 (at the centerline of pipe)
"full" signifies that there is a boundary layer at y=H
This grid will be symmetric about y = 0.5H


Sample usage to generate a 3D cylindrical grid that has a boundary layer at the wall:
pipe_grid.py 11 11 11 0.0 1.0 4.0 1.0e-3 0.05 3d 
11 points in the x direction
11 points in the y direction(boundary layer direction)
11 points in the theta direction
Inner radius of 0.0
Outer radius of 1.0
A length of 4.0
An initial cell spacing of 1.0e-3
A final cell spacing of 0.05 (at the centerline of pipe)
"3d" signifies that the output mesh will be a full cylinder mesh 


To generate a Loci VOG formatted grid use the following command:
plot3d2vog -m file     (the utility assumes the .grd suffix)


