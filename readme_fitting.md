How to use the GUI interface

### 3. Fitting Curved Folded Model to the User Specified Points

#### 3.1. Load Curve Parameters (P.txt)<br>

[load] button -> choose P.txt 

<img src="readme_images/02_LoadP.bmp" width="100%"><br>

All the input files are in "CurvedFold/GUI/input/".<br>
In P.txt, the folding angle and the 2D curvature of the control points on the 2D crease curve are defined on 3rd and 4th columns.<br>
m2m3.txt is automatically loaded with P.txt.  It describes the 3D and 2D poses of the starting point of the curved crease in the form of homogeneous coordinates.

| P.txt | m2m3.txt ||
| ---- | ---- | ---- |
| 10 1 # file type, mode: B, kv,tr,fa,k2d &nbsp; &nbsp; &nbsp; <br>7 # plot count<br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>0</sub> &nbsp; &kappa;<sub>2D 0</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>1</sub> &nbsp; &kappa;<sub>2D 1</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>2</sub> &nbsp; &kappa;<sub>2D 2</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>3</sub> &nbsp; &kappa;<sub>2D 3</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>4</sub> &nbsp; &kappa;<sub>2D 4</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>5</sub> &nbsp; &kappa;<sub>2D 5</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>6</sub> &nbsp; &kappa;<sub>2D 6</sub> | m2 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;<br> X<sub>x</sub> &nbsp; X<sub>y</sub> &nbsp; 0.0<br>Y<sub>x</sub> &nbsp; Y<sub>y</sub> &nbsp; 0.0<br>P<sub>x</sub> &nbsp; P<sub>y</sub> &nbsp; 1.0<br>m3<br>X<sub>x</sub> &nbsp; X<sub>y</sub> &nbsp; X<sub>z</sub> &nbsp; 0.0<br>Y<sub>x</sub> &nbsp; Y<sub>y</sub> &nbsp; Y<sub>z</sub> &nbsp; 0.0<br>Z<sub>x</sub> &nbsp; Z<sub>y</sub> &nbsp; Z<sub>z</sub> &nbsp; 0.0<br>P<sub>x</sub> &nbsp; P<sub>y</sub> &nbsp; P<sub>z</sub> &nbsp; 1.0<br> | <img src="readme_images/209_MsCPs.bmp" width="300"> |

#### 3.2. Load Target Points (target*.txt)<br>

Click [tgt_pt] button and choose target*.txt.<br>
3D target points and the corresponding points on the curved fold are loaded.<br>
The color shows the distance between 2D and 3D points.<br>
Red/blue indicates larger distance and green for smaller distance.

<img src="readme_images/04_LoadTarget.bmp" width="100%"><br>

target*.txt includes pairs of 3D target point positions, 1st to 3rd columns, and 2D control point positions, the 4th and 5th columns.

| target*.txt |  |
| ---- | :----: |
| tp0<sub>x</sub> &nbsp; tp0<sub>y</sub> &nbsp; tp0<sub>z</sub> &nbsp; cp0<sub>x</sub> &nbsp; cp0<sub>y</sub><br>tp1<sub>x</sub> &nbsp; tp1<sub>y</sub> &nbsp; tp1<sub>z</sub> &nbsp; cp1<sub>x</sub> &nbsp; cp1<sub>y</sub><br>tp2<sub>x</sub> &nbsp; tp2<sub>y</sub> &nbsp; tp2<sub>z</sub> &nbsp; cp2<sub>x</sub> &nbsp; cp2<sub>y</sub><br>tp3<sub>x</sub> &nbsp; tp3<sub>y</sub> &nbsp; tp3<sub>z</sub> &nbsp; cp3<sub>x</sub> &nbsp; cp3<sub>y</sub><br>tp4<sub>x</sub> &nbsp; tp4<sub>y</sub> &nbsp; tp4<sub>z</sub> &nbsp; cp4<sub>x</sub> &nbsp; cp4<sub>y</sub><br>... | <img src="readme_images/00_CPsTs.bmp" width="300"> |

#### 3.3. Load Target Mask (tmask*.txt) (optional)

Click [tmask] button and choose tmask*.txt.  Target points on mask=0 are deleted.
<img src="readme_images/06_LoadTargetMask.bmp" width="100%">

| tmask*.txt<br> |
| ---- |
| 1 1 1 1 1 0 0 0 0<br> 1 1 1 1 1 0 0 0 0<br> 1 1 1 1 1 0 0 0 0<br> 1 1 1 1 1 0 0 0 0<br> 1 1 1 1 1 1 1 1 1<br> 0 0 0 0 1 1 1 1 1<br> 0 0 0 0 1 1 1 1 1<br> 0 0 0 0 1 1 1 1 1<br> 0 0 0 0 1 1 1 1 1<br> |

#### 3.4. Optimize the Shape and the Pose of the Curved Fold Model

Click [Rulings and Angles] button under "OPTIMIZATION" section on the tab [D4]<br>
The shape and the pose of the curved fold model is optimized to fit the target points.
<img src="readme_images/07_Optimize.bmp" width="100%"><br>

