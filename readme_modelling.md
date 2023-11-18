How to use the GUI interface

### 2. Modelling a Curved Folded Shapes and the Folding Motion

#### 2.1. Display Option Settings<br>
Check the checkbox[CPT] to show the control points on the crease curve.

<img src="readme_images/201_DisplaySettings.bmp" width="100%"><br>

#### 2.2. Load Initial Configuration<br>

Load "CurvedFold/GUI/input/P.txt".<br>
This step has been done when the software is launched.<br>
To reset the model shape, reload the same file.
<img src="readme_images/202_LoadInit.bmp" width="100%"><br>

#### 2.3. Choose the Control Mode, Parameter, and the Control Point to Modify the Parameter<br>

<img src="readme_images/203_ChooseMode.bmp" width="100%"><br>

#### 2.4. Adjust the Parameter<br>

Choose the control point and change its parameter alternately to adjust the model shape as desired.
<img src="readme_images/204_ChangeParam.bmp" width="100%"><br>

#### 2.5. Change the Position and the Angle of the 2D Crease Curve<br>
- Rotation: Left mouse drag on 2D pane.
- Translation: Right mouse drag on 2D pane.

<img src="readme_images/205_2D_PosAngle.bmp" width="100%"><br>

#### 2.6. Change the Position and the Pose of the 3D Model<br>
- Rotation : shift + right mouse drag on 3D pane
- Translation : ctrl + right mouse drag on 3D pane

<img src="readme_images/206_3D_Pose.bmp" width="100%"><br>

#### 2.7. Folding Motions<br>

<img src="readme_images/207_FoldingMotion1.bmp" width="100%"><br>
<img src="readme_images/207_FoldingMotion2.bmp" width="100%"><br>

#### 2.8. Export the Configuration<br>
By clicking on the button [log] under "SAVE" section, following files are saved in "CurvedFold/GUI/output/HHMMSS/*".
- CP.obj : 3D OBJ Polygon model
- CP1.svg : 2D Crease pattern + rulings
- P.txt, m2m3.txt : Input file of this software. 
- model.bmp : screenshot of 3D pane
- X.csv : parameters on the crease curve
- quad.csv, d180.csv : errors of the 3D model

<img src="readme_images/208_Save.bmp" width="100%"><br>

P.txt contains the control mode and the curvature, the torsion, the folding angle, and the 2D curvature of the control points.<br>
m2m3.txt describes the 3D and 2D poses of the starting point of the curved crease, in the form of homogeneous coordinates.

| P.txt | m2m3.txt ||
| ---- | ---- | ---- |
| 10 1 # file type, mode: B, kv,tr,fa,k2d &nbsp; &nbsp; &nbsp; <br>7 # plot count<br>&kappa;<sub>0</sub> &nbsp; &tau;<sub>0</sub> &nbsp; &alpha;<sub>0</sub> &nbsp; &kappa;<sub>2D 0</sub><br>&kappa;<sub>1</sub> &nbsp; &tau;<sub>1</sub> &nbsp; &alpha;<sub>1</sub> &nbsp; &kappa;<sub>2D 1</sub><br>&kappa;<sub>2</sub> &nbsp; &tau;<sub>2</sub> &nbsp; &alpha;<sub>2</sub> &nbsp; &kappa;<sub>2D 2</sub><br>&kappa;<sub>3</sub> &nbsp; &tau;<sub>3</sub> &nbsp; &alpha;<sub>3</sub> &nbsp; &kappa;<sub>2D 3</sub><br>&kappa;<sub>4</sub> &nbsp; &tau;<sub>4</sub> &nbsp; &alpha;<sub>4</sub> &nbsp; &kappa;<sub>2D 4</sub><br>&kappa;<sub>5</sub> &nbsp; &tau;<sub>5</sub> &nbsp; &alpha;<sub>5</sub> &nbsp; &kappa;<sub>2D 5</sub><br>&kappa;<sub>6</sub> &nbsp; &tau;<sub>6</sub> &nbsp; &alpha;<sub>6</sub> &nbsp; &kappa;<sub>2D 6</sub> | m2 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;<br> X<sub>x</sub> &nbsp; X<sub>y</sub> &nbsp; 0.0<br>Y<sub>x</sub> &nbsp; Y<sub>y</sub> &nbsp; 0.0<br>P<sub>x</sub> &nbsp; P<sub>y</sub> &nbsp; 1.0<br>m3<br>X<sub>x</sub> &nbsp; X<sub>y</sub> &nbsp; X<sub>z</sub> &nbsp; 0.0<br>Y<sub>x</sub> &nbsp; Y<sub>y</sub> &nbsp; Y<sub>z</sub> &nbsp; 0.0<br>Z<sub>x</sub> &nbsp; Z<sub>y</sub> &nbsp; Z<sub>z</sub> &nbsp; 0.0<br>P<sub>x</sub> &nbsp; P<sub>y</sub> &nbsp; P<sub>z</sub> &nbsp; 1.0<br> | <img src="readme_images/209_MsCPs.bmp" width="300"> |
