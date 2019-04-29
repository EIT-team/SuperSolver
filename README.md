# SuperSolver
Ye olde EIT super solver, for those times you have to run in it matlab, or you want some intermediate variables which you cant get out of PEITS.

Now runs as a function, 

```
% initialise all the structures
[Mesh,Fem,Fwd,Inv,Sol] = supersolver_init(vtx,tri,mat_ref,elec_pos,gnd_pos,prt_full);

```


See `Reconstruction_script` in the examples folder to generate a EIT Forward solution in a cylindrical tank. 
