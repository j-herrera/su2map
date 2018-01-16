# su2map
Maps a [SU2](https://github.com/su2code/SU2) CFD solution on VTK format onto a new mesh, creating a suitable restart file. The tool works for any cell shape in both 2D and 3D.

The syntax is: _python su2map.py flow.vtk new_mesh.su2 solution_flow.dat_