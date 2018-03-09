# Cloth simulation features
* Cloth<->object, cloth<->other cloth, and cloth<->self intersection
* Stretch and compress limits to get the cloth to behave more realistically (not like rubber)
* Bend springs to keep the cloth locally flat to achieve realistic looking folds
* Forward and Backward Euler integration methods
* simulation run in c++, frames printed out in .poly file format so they can fed into houdini for visualization

**Backward Euler Sim**<br />
![](BE.gif)


**Forward Euler Sim**<br />
![](FE.gif)



**Notes on collisions**<br />
* Cloth for collision purposes is simply represented as spheres
