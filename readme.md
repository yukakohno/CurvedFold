### General information

This source code is used for the experiments of CAD21 paper:
- Yuka Watanabe, Jun Mitani, "Fitting Single Crease Curved-Fold Model to the User Specified Points", Computer-Aided Design & Applications, 19(2), 2022, 387-404.

The source code and sample data for some applications listed in the paper are provided in other branches:
- "Intuitive User Interface for a Single Crease Model Manipulation" (Section 4.1 and Figure 2) : "SingleCrease2".<br>
https://github.com/yukakohno/CurvedFold/tree/SingleCrease2
- "Stitch the Boundaries of Two Surface Patches" (Section 4.2 and Figure 5) : "MultiPieces2".<br>
https://github.com/yukakohno/CurvedFold/tree/MultiPieces2

This code is build by Visual Studio 2019 with libraries below installed.
- openCV 2.4.9
- fltk-1.3.5
- freeglut 3.0.0 MSVC Package (https://www.transmissionzero.co.uk/software/freeglut-devel/)


### Some GUI widgets used in this work:

#### tab [D0]

--- LOAD ---<br>
[load]    load initial curved-fold model.<br>
[tgt_pt]  load target points from file input.<br>
[tmask]   load mask area of the target points aligned in 9x9.<br>
[rulings] load ruling angels.


#### tab [D4]

--- RULING 2 CURVE---<br>
scroll bars: change folding angles under given rulings.<br>
[switch]  switch ruling angles (among given 6 sets of predefined rulings)<br>
[Rul->TA] calculate torsions and folding angels under given rulings.<br>
[Opt Mat] optimize the pose of the curved fold model<br>
[Start],[Stop] Start/stop the folding animzation.<br>
[*]mat  if checked, optimize the pose while the folding animation is played<br>
[*]rot  if checked and []mat is unchecked, optimize the pose with fixed origin while the folding animation is played<br>

--- TARGET POINTS ---<br>
[add / modify / delete] acitivate the edit mode of surface control points<br>
idx [scroll bar] choose the surface control points to edit<br>

--- OPTIMIZATION ---<br>
[Angle with Fixed Rulings] optimize the folding angles under given rulings.<br>
[Random Rulings]           set random change in the rulings, optimize the folding angles and the pose according the checkbox settings.<br>
[*]optimization            if checked, adopt [Random Rulings] result only if the gap decreases.<br>
[Rulings and Angles]       optimize the curved-fold model based on the rulings and folding angles.<br>
[Random Torsion and Angles] try random torsions and folding angles.<br>
[Torsion and Angles]       optimize the curved-fold model based on the torsions and the folding angles.


#### tab[D5]

--- PARAM/TARGETS LIST ---<br>
[param list (torsion/angle)]  generate parameter lists based on random torsion and foldin angels<br>
[target list (torsion/angle)] generate target points lists based on random torsion and foldin angels<br>
[start][stop]                 show/stop the curved fold models of the generated parameters<br>
[param list (rulings)]        generate parameter lists based on random rulings<br>
[target list (rulings)]       generate parameter lists based on random rulings<br>
[start][stop]                 show/stop the curved fold models of the generated parameters<br>
[target masks]                choose an index of the mask(surface control points alignments)<br>

--- BATCH PROC ---<br>
input[]      index of input shape to start batch process (*1)<br>
tgt[]        index of target shape to start batch process (*1)<br>
mask[]       index of the mask (surface control points alignments) to start batch process (*1)<br>
[batch proc] start batch process (*2)<br>
[stop]       stop batch process (*1)<br>

(*1) valid only if build with<br>
  #define PROC_BATCH_IDLE 1 // 0:batch, 1:idle (ControlPanel_BatchProc.cpp line 67)<br>
(*2) the optimization method is defined by<br>
  #define PROC_OPT_MODE 2 // 0: ruling, 1: torsion & fold angle, 2: hybrid (ControlPanel_BatchProc.cpp line 68)
