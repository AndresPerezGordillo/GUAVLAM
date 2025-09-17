README
GUAVLAM 1.0

General Unsteady Aerodynamics Vortex-Lattice Method (GUAVLAM) is a code designed for unsteady aerodynamic analysis of lifting surfaces. Surfaces can be fixed (wings on an airplane) or rotary (rotor blades on a helicopter). Bluff bodies can be included, which serve as obstructions to the flow.

How do I get set up?
When you download this repository, the bulk of the source code is in the 'source/' directory. The code has currently been tested using the open95 Fortran compiler form the Open64 compiler suite. To compile the code:

Make sure you have the 'main.f90' file in your directory, as well as the Makefile. Copies of these can be found in the 'main/' directory and the 'NewPlaca/' directory.
Modify the Makefile to include the desired Fortran compilar and flags. The compiler should have support for Fortran 90.
Run the command 'make all'.
Run the executable ./guavlam
An example case is found in the 'NewPlaca/' directory. You can control parameters of the simulation using the 'CaseControl' file

Who do I talk to?
Juan D. Colmenares juaneco2710(at)gmail.com
