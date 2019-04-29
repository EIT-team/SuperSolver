# SuperSolver
Ye olde EIT super solver, for those times you have to run in it matlab, or you want some intermediate variables which you cant get out of PEITS.

Now runs as a function, 

```matlab
% initialise all the structures
[Mesh,Fem,Fwd,Inv,Sol] = supersolver_init(vtx,tri,mat_ref,elec_pos,gnd_pos,prt_full);
% setup system and boundary conditions
[Mesh,Fem,Fwd,Inv,Sol] = supersolver_setup(Mesh,Fem,Fwd,Inv,Sol);
% solve forward and get boundary voltages
[Mesh,Fem,Fwd,Inv,Sol,Data] = supersolver_runfwd(Mesh,Fem,Fwd,Inv,Sol);
% get jacobian
[Mesh,Fem,Fwd,Inv,Sol,J] = supersolver_solve(Mesh,Fem,Fwd,Inv,Sol);

```


See `Reconstruction_script` in the examples folder to generate a EIT Forward solution in a cylindrical tank. 
