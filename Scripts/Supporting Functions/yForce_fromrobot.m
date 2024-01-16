function Fpush = yForce_fromrobot(Fx,Fz,el_ang)

Fpush = -cosd(270 - el_ang)*Fx + sind(270 - el_ang)*Fz;