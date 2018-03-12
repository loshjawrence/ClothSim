# Cloth simulation features
* Cloth<->object, cloth<->other cloth, and cloth<->self intersection
* Stretch and compress limits to get the cloth to behave more realistically (not like rubber)
* Bend springs to keep the cloth locally flat to achieve realistic looking folds
* Forward and Backward Euler integration methods
* simulation run in c++, frames printed out in .poly file format so they can fed into houdini for visualization

# Cloth Setup
* NxM grid of vertices. 
* Structural springs exist between all adjacent neighbors (up, down, left, right, diagonals).
* Structural springs ensure the cloth exhibits cloth like behavior along the longitude, latitude, and shear directions.
* Bend springs exist between between every other vertex on up, down, left, right direction. 
* Bend springs ensure the cloth is locally flat by activating during compression by pushing away adjacent quad faces of the cloth.
* Stretch and compression limits are set on the cloth at an amount of 10% 
* This is done to make the cloth behave more like yarn material rather than a rubber material.
![](clothSetup.png)

# Notes on collisions
* For collision, the cloth particles are represented as spheres. 
* Collision is then checking every sphere against every other sphere (for both self collision and other cloth collision).
* If within range, if within range the particles are pushed apart so they are no longer touching.
* Velocity adjusted by subtracting off the component perpendicular to the plane of collision.
* This is similar to Granm-Schmidt orthogonalization
![](GramSchmidt.png)
![](collisionCloth.png)


**Backward Euler Sim (30x30 cloth)**<br />
![](BE3.gif)

**Backward Euler Sim**<br />
![](BE.gif)

**Forward Euler Sim**<br />
![](FE.gif)

**Backward Euler Sim (30x30 cloth)**<br />
![](BE3.gif)

# Thoughts on FE vs BE
* FE is still much faster even when taking many more steps than BE.
* BE code is more complex than FE.
* BE must begin to include capsule capsule collision tests at larger time steps.
* BE must do more complex collision testing and adjustments when taking large time steps.
**Graphs**<br />
![](graph1.png)
![](graph2.png)

