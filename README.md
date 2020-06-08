# pion_picnic
This Repository is an extension of my previous Pion-Code to calculate the mass of the Pion. The extension uses the FULL Dirac-Structure of the Pion (with 4 Dressing Functions E-H). Furthermore, this time "form" was used to genereated the desired/required Dirac-Traces. 

How to use:
- cd to build_main.
- execute cmake with link to all files in folder above.
- make executeable 'maintest2' and move it to folder above (where all the links to libraries etc work).
- cd to folder above.
- execute maintest2 to calculate the pion mass.
```bash
cd build_main
cmake ..
make && mv maintest2 ..
cd ..
./maintest2
```
- be happy

On first execution, make sure the parameter
```C++
#define iterate_quark
```
is turned ON (not commented out). This will generate a file needed for the iteration/extrapolation of the pion.
Once this file is generated the parameter above can be turned off, unless you want to experiment with the numerics. Then you can change the paramters
```C++
#define absciss_points 200
#define ang_absciss_points 100
#define ang_absciss_points_2 16
#define BSE_absciss_points 64
#define BSE_ang_absciss_points 16
#define max_step 2000
#define max_iter 200
```
and to your liking.
CHEERS
