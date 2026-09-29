This is a simple handbook perameterised design method for creating 
the intial design for an clasic airhsip, capable of performing stable flight.

The purpos of this work is to demonstate how even huristic perameterised 
methods can help imporve intuition over a problem.


Inputs:
Hull
1. Desired length & Diametier
3. Skin weight

Fin
LeadingEdgeTapperAngle_rad = deg2rad(15);
tailingEdgeHight_m = 1;
xPosTailEdge_m = 11.75;
dihedral_rad = deg2rad(50);
numbFins = 3;


Intuition learning inputs
1. Lifting gas (Helium, hydrogen)
2. Altitude
3. Tempterature
4. Fin mass
5. Add balist to balance ship (Yes No)


It returns:
1. Blimp geometry
2. Fin sizing
3. Lift
4. Inirtia matrices
5. Center of mass
6. static moment
7. Max flight Altitude



Note: As this model is focused on dynamic and static stability, we have not 
considered propolsion systems or payload. We also do not modle airodynamic 
drag or lift.