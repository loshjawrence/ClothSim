# Cloth simulation features
* Cloth<->object, cloth<->other cloth, and cloth<->self intersection
* Stretch and compress limits to get the cloth to behave more realistically (not like rubber)
* Bend springs to keep the cloth locally flat to achieve realistic looking folds
* Forward and Backward Euler integration methods
* simulation run in c++, frames printed out in .poly file format so they can fed into houdini for visualization

**Cloth Setup**<br />
![](clothSetup.png)

**Backward Euler Sim**<br />
![](BE.gif)

**Forward Euler Sim**<br />
![](FE.gif)

**Backward Euler Sim (30x30 cloth)**<br />
![](BE3.gif)

**Graphs**<br />
![](graph1.png)
![](graph2.png)


**Notes on collisions**<br />
* For collision, the cloth particles are represented as spheres. Collision is then checking every sphere against every other sphere (for both self collision and other cloth collision).
