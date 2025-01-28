
## Hi

This is some practice writing a graph slam backend from scratch in C to gain a solid
understanding of how it works. I will try to keep it to no libraries, and do all
the computation in house so to speak. The end goal would be to use this backend
in a full slam system and see how it does.

## Some more details

This is a 2d graph slam algorithm that holds nodes with values x, y, and theta. As
a user you must supply the parameters for the nodes. Observation constraints take 
the form of xj <- xi, that is: how would xj look relative to xi. In homogeneous 
coordinates, this is the transformation Z which has the proporties Z-1 * xj = xi
and Z * xi = xj. The graph is expected to be built up as the user traverses the
environment. When a loop closure occurs, or after multiple have occured, the 
optimization can be triggered by the user, with an arbitrary node set as static.
This will modify all the poses held by the pose graph, and the user can extract
them to update thier own representation when the process is done.
