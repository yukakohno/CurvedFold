How to use the GUI interface

### 3. Fitting Curved Folded Model to the User Specified Points

#### 3.1. Load Curve Parameters (P.txt)<br>

[load] button -> choose P.txt 

<img src="readme_images/02_LoadP.bmp" width="100%"><br>

All the input files are in "CurvedFold/GUI/input/".<br>
In P.txt, the folding angle and the 2D curvature of the control points on the 2D crease curve are defined on 3rd and 4th columns.

| P.txt |  |
| ---- | :----: |
| 10 1 # file type, mode: B, kv,tr,fa,k2d<br>7 # plot count<br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>0</sub> &nbsp; &kappa;<sub>2D 0</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>1</sub> &nbsp; &kappa;<sub>2D 1</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>2</sub> &nbsp; &kappa;<sub>2D 2</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>3</sub> &nbsp; &kappa;<sub>2D 3</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>4</sub> &nbsp; &kappa;<sub>2D 4</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>5</sub> &nbsp; &kappa;<sub>2D 5</sub><br>0.0 &nbsp; 0.0 &nbsp; &alpha;<sub>6</sub> &nbsp; &kappa;<sub>2D 6</sub> | <img src="readme_images/00_CPs.bmp" width="300"> |


<br>

As P.txt is loaded, m2m3.txt is automatically loaded.  It describes the poses, in the form of homogeneous coordinates, of the starting point of the curved crease in the 2D paper and the 3D space.

| m2m3.txt |  |
| ---- | :----: |
| m2 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;<br> X<sub>x</sub> &nbsp; X<sub>y</sub> &nbsp; 0.0<br>Y<sub>x</sub> &nbsp; Y<sub>y</sub> &nbsp; 0.0<br>P<sub>x</sub> &nbsp; P<sub>y</sub> &nbsp; 1.0<br>m3<br>X<sub>x</sub> &nbsp; X<sub>y</sub> &nbsp; X<sub>z</sub> &nbsp; 0.0<br>Y<sub>x</sub> &nbsp; Y<sub>y</sub> &nbsp; Y<sub>z</sub> &nbsp; 0.0<br>Z<sub>x</sub> &nbsp; Z<sub>y</sub> &nbsp; Z<sub>z</sub> &nbsp; 0.0<br>P<sub>x</sub> &nbsp; P<sub>y</sub> &nbsp; P<sub>z</sub> &nbsp; 1.0<br> | <img src="readme_images/00_Ms.bmp" width="300"> |


#### 3.2. Load Target Points (target*.txt)<br>

[tgt_pt] button -> choose target*.txt 

<img src="readme_images/04_LoadTarget.bmp" width="100%"><br>

3D target points and the corresponding points on the curved fold are loaded.
The color shows the distance between 2D and 3D points.
Red/blue indicates larger distance and green for smaller distance.

target*.txt includes pairs of 3D target point positions, 1st to 3rd columns, and 2D control point positions, the 4th and 5th columns.

| target*.txt |  |
| ---- | :----: |
| tp0<sub>x</sub> &nbsp; tp0<sub>y</sub> &nbsp; tp0<sub>z</sub> &nbsp; cp0<sub>x</sub> &nbsp; cp0<sub>y</sub><br>tp1<sub>x</sub> &nbsp; tp1<sub>y</sub> &nbsp; tp1<sub>z</sub> &nbsp; cp1<sub>x</sub> &nbsp; cp1<sub>y</sub><br>tp2<sub>x</sub> &nbsp; tp2<sub>y</sub> &nbsp; tp2<sub>z</sub> &nbsp; cp2<sub>x</sub> &nbsp; cp2<sub>y</sub><br>tp3<sub>x</sub> &nbsp; tp3<sub>y</sub> &nbsp; tp3<sub>z</sub> &nbsp; cp3<sub>x</sub> &nbsp; cp3<sub>y</sub><br>tp4<sub>x</sub> &nbsp; tp4<sub>y</sub> &nbsp; tp4<sub>z</sub> &nbsp; cp4<sub>x</sub> &nbsp; cp4<sub>y</sub><br>... | <img src="readme_images/00_CPsTs.bmp" width="300"> |



#### 3.3. Load Target Mask (tmask*.txt) (optional)

[tmask] button -> choose tmask*.txt 

<img src="readme_images/06_LoadTargetMask.bmp" width="100%">

Target points on mask=0 are deleted.

> 1 1 1 1 1 0 0 0 0<br>
> 1 1 1 1 1 0 0 0 0<br>
> 1 1 1 1 1 0 0 0 0<br>
> 1 1 1 1 1 0 0 0 0<br>
> 1 1 1 1 1 1 1 1 1<br>
> 0 0 0 0 1 1 1 1 1<br>
> 0 0 0 0 1 1 1 1 1<br>
> 0 0 0 0 1 1 1 1 1<br>
> 0 0 0 0 1 1 1 1 1<br>

#### 3.4. Load Rulings (ruling*.txt) (optional)

[ruling] button -> choose ruling*.txt 

<img src="readme_images/08_LoadRulings.bmp" width="100%">

Rulings are loaded.

In ruling*.txt, the ruling angles of the control points on the left and the right sides of the 2D crease curve are defined.

> 2.022672	1.991598	1.849320	1.562967	1.298554	1.156659	1.127758<br>
> 1.118921	1.149995	1.292273	1.578625	1.843039	1.984934	2.013835<br>


#### 3.5. Switch Rulings (optional)

The rulings may be switched by [switch] button on [D4] tab.

<img src="readme_images/09_switch_ruling00.bmp" width="49%"> <img src="readme_images/09_switch_ruling01.bmp" width="49%">
<img src="readme_images/09_switch_ruling02.bmp" width="49%"> <img src="readme_images/09_switch_ruling03.bmp" width="49%">
<img src="readme_images/09_switch_ruling04.bmp" width="49%"> <img src="readme_images/09_switch_ruling05.bmp" width="49%">
<img src="readme_images/09_switch_ruling06.bmp" width="49%">

#### 3.6. Start/Stop Folding Animation (w/wo pose optimization)

The folding motion starts/stops by [Start][Stop] buttons.
The folding angle is increased or decreased with fixed rulings.
If you check the checkbox [mat], the pose is set to minimize the total distance to the target points, or the distances between the target points and the corresponding points on the curved fold.

<img src="readme_images/10_start_stop_00.bmp" width="100%">
<!--img src="readme_images/10_start_stop_01.bmp" height="35%">  <img src="readme_images/10_start_stop_02.bmp" height="35%"-->

#### 3.7. Optimize Folding Angle with Fixed Rulings

[Angle with Fixed Rulings] button sets the folding angle to minimize the total distance to the target points.

<img src="readme_images/11_Angle_w_Fixed_Rulings.bmp" width="100%">

#### 3.8. Random Rulings

[Random Rulings] button makes a small changes int the rulings.
If the the checkbox [optimze] is checked, the rulings change is accepted only if the total distance to the target points are decreased by the change.

<img src="readme_images/12_Random_Rulings_00.bmp" width="100%">
<!--img src="readme_images/12_Random_Rulings_01.bmp" width="75%"-->
<!--img src="readme_images/12_Random_Rulings_02.bmp" width="75%"-->

#### 3.9. Optimize Rulings and Angles

[Rulings and Angles] button optimizes the rulings and the folding angles to minimized the distance to the target points, so that the curved folded surface approximates the target points.

<img src="readme_images/13_Rulings_and_Angles.bmp" width="100%">
