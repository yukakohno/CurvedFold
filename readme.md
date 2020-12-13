### General information

This is a software to model and design a curved folding model with a single crease and a "rotationally symmetric curved folding" composed of mutiple pieces with single creases.
The parameters of the curved fold are optimized to fit the surface to the 3D target points.<br>

This code is build by Visual Studio 2019 with libraries below installed.<br>
- openCV 2.4.9
- fltk-1.3.5
- freeglut 3.0.0 MSVC Package (https://www.transmissionzero.co.uk/software/freeglut-devel/)

For the basic information of this software, please refer to:<br>
Yuka Watanabe, Jun Mitani, "Modelling the Folding Motions of a Curved Fold", in Origami7: Proceedings of the 7th International Meeting on Origami in Science, Oxford, England, September 5-7, 2018, pages 1135-1150.<br>
Yuka Watanabe, Jun Mitani, "Visualization of Folding Motion of Rotationally Symmetric Curved Folding", Computer-Aided Design & Applications, 17(3), 2019, 513-524.<br>

### How to optimize the parameters

#### Curved folding with a single crease

See readme.md in branch "SingleCrease2."

#### Curved folding with multiple pieces with single creases


###### 1. Load Curve Parameters (P.txt)<br>

[load] button -> choose P1.txt<br>
Adjust rotation, scale, translation of the 3D object.

<img src="readme_images/01_LoadP1.bmp" width="100%"><br>

###### 2. Make multiple pieces <br>

Choose DIVISION NUMBER in [D3] tab to dupulicate the 3D object and to make "rotationally symmetric curved folding."

<img src="readme_images/02_DivisionNumber.bmp" width="100%"><br>

###### 3. Choose an intermediate frame <br>

Choose an intermediate frame by choosing a number between 0 and -20 of the value slider FOLDING MOTION.

<img src="readme_images/03_IntermediateState.bmp" width="100%"><br>

###### 4. Optimize the pose of the piece.<br>

Optimize and set the pose of the pieces by
 [opt] button -> [set] button -> [opt] button in matrix row.<br>

<img src="readme_images/04_OptimizeMatrix.bmp" width="100%"><br>

Below DISPLAY in [D0] tab,<br>
Check the checkbox [ONE] to show one piece with texture and others as wire frame models.<br>
Check the checkbox [stch] to show the pairs of boundary points "stitched" and their mid points.

<img src="readme_images/05_dispStitch.bmp" width="100%"><br>

###### 5. Set target points <br>

Set the mid points of the pairs of boundary points as the target points by
 [opt] button in param row.<br>
Target point are shown in blue and the original points, or the boundary points on one side, is shown in red.

<img src="readme_images/06_SetTarget.bmp" width="100%"><br>

###### 6. Optimize the parameters, the torsions and the folding angles<br>

On [D4] tab, check the checkbox [rot] and press [Rulings and Angles] button.<br>
If the optimization process does not start, press [Random Rulings] button beforehand.<br>
In the optimization process, the torsions and the folding angels of the curved crease are optimized to minimize the total distances between the target points and the original points.

<img src="readme_images/07_OptimizeRulingsAngles.bmp" width="100%"><br>

###### 7. Set the optimized parameters <br>

Go back to [D3] tab and set the optimized parameters by [set] button in param row.<br>

<img src="readme_images/08_SetParam.bmp" width="100%"><br>

###### 8. Check the result <br>

Check the result.  Repeat step 4-7 if necessary.<br>

<img src="readme_images/09_Result.bmp" width="100%"><br>
