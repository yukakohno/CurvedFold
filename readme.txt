This code is build by Visual Studio 2019 with libraries below installed.
-openCV 2.4.9
-fltk-1.3.5
-freeglut 3.0.0 MSVC Package (https://www.transmissionzero.co.uk/software/freeglut-devel/)

This source code is used for the experiments of CAD21 paper:
Yuka Watanabe, Jun Mitani, "Fitting Single Crease Curved-Fold Model to the User Specified Points", Computer-Aided Design & Applications (to appear).

Some GUI widgets related to this work is listed below:

tab[D0]

--- LOAD ---
[load]    load initial curved-fold model.
[tgt_pt]  load target points from file input.
[tmask]   load mask area of the target points aligned in 9x9.
[rulings] load ruling angels.


tab[D4]

--- RULING 2 CURVE---
scroll bars: change folding angles under given rulings.
[switch]  switch ruling angles (among given 6 sets of predefined rulings)
[Rul->TA] calculate torsions and folding angels under given rulings.
[Opt Mat] optimize the pose of the curved fold model
[Start],[Stop] Start/stop the folding animzation.
[*]mat  if checked, optimize the pose while the folding animation is played
[*]rot  if checked and []mat is unchecked, optimize the pose with fixed origin while the folding animation is played

--- TARGET POINTS ---
[add / modify / delete] acitivate the edit mode of surface control points
idx [scroll bar] choose the surface control points to edit

--- OPTIMIZATION ---
[Angle with Fixed Rulings] optimize the folding angles under given rulings.
[Random Rulings]           set random change in the rulings, optimize the folding angles and the pose according the checkbox settings.
[*]optimization            if checked, adopt [Random Rulings] result only if the gap decreases.
[Rulings and Angles]       optimize the curved-fold model based on the rulings and folding angles.
[Random Torsion and Angles] try random torsions and folding angles.
[Torsion and Angles]       optimize the curved-fold model based on the torsions and the folding angles.


tab[D5]

--- PARAM/TARGETS LIST ---
[param list (torsion/angle)]  generate parameter lists based on random torsion and foldin angels
[target list (torsion/angle)] generate target points lists based on random torsion and foldin angels
[start][stop]                 show/stop the curved fold models of the generated parameters
[param list (rulings)]        generate parameter lists based on random rulings
[target list (rulings)]       generate parameter lists based on random rulings
[start][stop]                 show/stop the curved fold models of the generated parameters
[target masks]                choose an index of the mask(surface control points alignments)

--- BATCH PROC ---
input[]      index of input shape to start batch process (*1)
tgt[]        index of target shape to start batch process (*1)
mask[]       index of the mask (surface control points alignments) to start batch process (*1)
[batch proc] start batch process (*2)
[stop]       stop batch process (*1)

(*1) valid only if build with
  #define PROC_BATCH_IDLE 1 // 0:batch, 1:idle (ControlPanel_BatchProc.cpp line 67)
(*2) the optimization method is defined by
  #define PROC_OPT_MODE 2 // 0: ruling, 1: torsion & fold angle, 2: hybrid (ControlPanel_BatchProc.cpp line 68)
