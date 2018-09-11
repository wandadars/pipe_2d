Sample usage to generate a grid that has a boundary layer at y = 0:
python3 pipe-2D.py 7 9 1 1 5.0e-4 3.0e-6 0.1 half

grid has:
7 points in the x direction
9 points in the y direction(boundary layer direction)
A length of 1
A height of 1
A width of 5.0e-4
An initial cell spacing of 3.0e-6
A final cell spacing(at y=1) of 0.1
"half" signifies that there is no boundary layer at y=1



Sample usage to generate a grid that has a boundary layer at y = 0 and y = 1:
python3 pipe-2D.py 7 9 1 1 5.0e-4 3.0e-6 0.1 full

grid has:
7 points in the x direction
9 points in the y direction(boundary layer direction), 9 points spans the entire vertical space
A length of 1
A height of 1
A width of 5.0e-4
An initial cell spacing of 3.0e-6
A final cell spacing(at y=1) of 0.1
"full" signifies that there is a boundary layer at y=1
This grid will be symmetric about y = 0.5
