
# Graph Slam Backend

## Description

This is a graph slam backend written from scratch in C, mainly to help me gain
understanding of how it works. No libraries (other then stdlib) are used, and
the project currently contains a full working backend for a 2d (3 dof) system.
I will note that in the description I assume that the user has basic knowledge
of robotics, and focus primarily on the description of the system. The source
code is heavily documented, and all sources used are provided below.

## Build && Install

Build the project and ensure that it is passing all tests.

```
$ mkdir build
$ make 
$ make test
```

Clone the repository and execute the following command:
```
$ sudo make install
```
Note: The install command assumes that you are using a unix-like file system
and will fail on Windows or MacOS.   
   
There is also a CMakeLists.txt for building this project as a package
to be used in c++ cmake projects. This can be used in the following 
way:

```
$ mkdir build && cd build
$ cmake ..
$ sudo make install
```

This allows the find_package call to be used to include the 
project instead of target_include_directories:
```

/* using cmake build & install */

find_package(GraphSlamBackend REQUIRED)
...
target_link_libraries(my_lib GraphSlamBackend::graph_slam_backend)

/* using make build & install */

target_link_libraries(my_lib graph_slam_backend)

```

It is kind of redundant to have two build systems, but I spent
a while googling on how to do cmake packages, so I have left it
in here partly for future reference for myself.

## Usage

Provided you have completed the previous steps, you can include
the backend in your project using:
```

/* in c projects */

#include <graph_slam_backend/graph_slam.h>

/* in c++ projects */

extern "C" {
#include <graph_slam_backend/graph_slam.h>
}

```

### Setup

Getting the system setup is pretty simple. Call
```
    pg * pose_graph = construct_pose_graph();
```
to construct a pose graph. This graph will hold
all the poses and update them as the graph is
optimized. In order to take advantage of this, you
will have to loop through the internal node list
of this structure when you want to work with the
poses in it. Doing this looks something like this:
```

/* note: structure of type(pv) */
typedef /* ... */ {
    float x; /* the x position in cartesian coordinates of the node */
    float y; /* the y position ''                                   */
    float t; /* the rotation at the node in radians, theta          */
} pv;

    for (int i = 0; i < pose_graph->node_count; ++i) {

        pv pv_value = pose_graph->node[i].pos;

        /* convert pv_value into my representation of position */
        (/* my representation */).pos = pv_value;
    }

```

### Nodes

As is standard in a graph, our pose graph will be composed of nodes
and edges. Each node can be thought of as a variable, which is constrained
by the edges in the graph (covered next). A good visual analogy of this
is that each node is bungee-corded or rubber-banded to a set of rods
surounding it. The node will move to the position that equally distributes
the pulling of the bands. Here the bands are the constraints, and the rods
are actually other nodes, but it is better in the analogy to think of them
as fixed. We add more nodes to our system because we will benefit from
what the final state of these nodes means for our system. Often times
these nodes will end up holding other information, like the actual position
of the robot. If so, we need to know where the robot is when the node came
to rest. Having other nodes helps us know where the robot has been and 
allows us to potentially build a map using that information. Nodes can be
generated from any observation or measurement about the position of the robot,
but they need at least one edge to constrain them. Otherwise it is basically like
trying to get the node to hover in the air with no band linking it. This node
would be unconstrained, and would result in a failure of the system to find
an optimal positioning of the node, given that all positions would have the
same pull of zero! This is a very loose treatment of the graph slam reasoning,
a more detailed treatment is given in the sources at the bottom of this document,
which are recomended if you want an understanding of why this system works.   
   
Moving forward. The addition of a node will look something like this:
```

    /* collect some measurement of the environment */
    /* ... */ measurement = /* measurement function */

    /* somehow convert the measurement to a pv */
    pv pos = /* conversion */(measurement);

    pgn_handle handle = add_node(pos,graph);

    /* store handle with related observations 
     * for constructing the path of the robot
     * and adding edges in the future        */
```

Note again the conversion of the collected data to a
position vector (pv). There will be a lot of that kind of stuff.
Note also that a pgn_handle is returned. This value is an integer
that holds the position in the node array shown in the previous
section that our added node occupies. This value should be kept
track of by the user, as it is needed for adding edges to the
graph.

### Edges

These are the bungee cords/rubber bands mentioned in the previous
section. The more of these we add, the more constrained our system
will be. Note again that all nodes need at least one edge, and more
edges are usually better for results than less. Be careful though,
edges should hold positional data you are more certain of, as if you
add an edge that is incorrect, it will drag your node in the wrong
direction. There are a few important pieces to an edge. The first
is the handles of the two nodes involved. The first node is the
reference node. It basically says that we are viewing the second
node from it's perspective. It is important to keep that in mind,
as it affects the meaning of the remaining components. If you get these
values mixed up, your pose graph will optimize backwards, or diverge,
which is not good. The next element is the measurement. This is a 
position vector (pv) and should be a relative measurement based on
the first node. In other terms, this should be the position of the
second node provided the first one is defined as the origin of the
coordinate system. In other terms, this is the change in rotation,
x position, and y position of the second node from the first node.
Note that if the second and first nodes we swapped, this value
would be inverted. The final element in an edge is the information
matrix. This is a symmetric 3x3 matrix and defined as the inverse
of the covariance matrix. Wikipedia is a good resource if the terms
are unfamiliar, but I will summerize it simply here. The covariance 
describes how much the measurement can be trusted, and encodes this
as values between zero and 1, with higher values indicating more
certainty. As the information matrix is the inverse, it does not have
to be bounded by the \[0,1] range, and can have arbitrary values,
as long as the matrix remains symmetric and positive semi-definite. So
no negative or zero values on the diagonal. This value should be
estimated directly when an correlation or edge is detected between
two nodes, or should be collected by computing the inverse of a 
covariance that was estimated directly. Once these elements have been
assembled, an edge can be added to the pose graph. The syntax is:
```

/* Definition of sm33: */
typedef /* ... */ {
    float a00;
    float a10, a11;
    float a20, a21, a22;
} sm33


    pgn_handle A;
    pgn_handle B;

    /* assume an observation has detected
     * that A is related to B somehow. This
     * would typically be done via ICP, scan
     * matching, photogrammetry, odometry, or
     * some other similar system.            */

    pv movement = /* compute reference position
                   * B given A using observation */

    sm33 information = /* compute the information matrix of this
                        * observation using the probability the
                        * data is correct at each vector element */

    add_edge(A,B,&movement,&information,graph);

```
Once we have collected enough edges, and usually after we have returned
to a previously visited location and produced edges that relate nodes
that are many timesteps away from ones at our current position, we can
move to the optimization phase.

### Optimization

Optimization is the most important part of the graph slam algorithm. It
works by moving the nodes around until they meet the constraints introduced
by the edges as best as possible. This can be achived in many ways (even
through random guessing!) but here it is achived using the sparse cholesky
decomposition, given the observation that the information matrix is constructed
in a way that always results in a positive semi-definite matrix, and that the
number of edges for each node is often much less than the total number of
nodes, i.e. the graph has very few connections. To use this method, a sparse
matrix with a set of constaints is constructed, and is then solved repeatedly
until the solutions are not any different from those given in the previous
iteration, at which point the function returns. The reason you are being told
all of this is to understand that the success of the optimization function depends
on the correct construction of the information matrix. If you do not provide a
sufficient number of edges (at least one per node!), have massive values in
your information matrix or have a non positive seim-definite information matrix,
or perform other similarly poor actions, the call to optimize will diverge, hang
or cause a crash of the program. There is little that can really be done by
the function to prevent this, as it does not have the power to tell if an
information matrix is correctly constructed. It can, however, check for divergence
and will print a warning message and exit with the value -1 if it detects
a divergence occuring. This should be seen as a message that your setup
of edges is incorrect and needs fixing. It could also be an error in my code. I
am not perfect. However, this system has worked on the inputs I have tested, and
I have confidence in it. I would suggest looking at your code before digging for
errors in the solver. The call to optimize is rather simple, provided everything
else is correct:
```
    int res = optimize(graph,NULL);

    if (res == -1)
        /* system diverged */
    else if (res == 0)
        /* success */
```

After a successful optimization, the node positions can be recovered using
the syntax discussed in the Setup section. These optimized nodes can be used
to update the current position of the robot, or to compute a new map given
the observations. Note that this call may take a while to return, as it is
doing a lot of work. If this becomes a problem in practice, refer to the
below section for information on limiting the runtime of the optimization
and on collecting information about the call.

#### Modifying Behavior

If you have constraints on the time the optimization can run, or want
more information about how long the process takes, you can pass a gsod
(graph_slam_optimization_data) structure in as an optional second argument.
This allows you to set additional exit conditions, change the existing ones,
and collect more detailed information about why the optimization terminated.
Here are some examples of usage:
```
    /* change the cutoff value for exit */
    /* see source code for details      */
    gsod settings;
    bzero(&settings,sizeof(gsod));
    settings.cutoff = 1-e2;

    optimize(graph,&settings);
```
```
    /* set the limit of steps the optimization
     * can perform. This may be useful if the
     * system is time constrained and you want
     * to perform the optimization in stages */
    gsod settings;
    bzero(&settings,sizeof(gsod));
    settings.step_limit = 10;

    optimize(graph,&settings);

    if (settings.end_state = OES_step_limit) {
        /* the call to optimize exited due to hitting
         * the step limit, not because of success or
         * divergence.                              */
    }

```
```
    /* set the limit of time the optimization
     * can perform. This may be useful if the
     * system is time constrained and you want
     * to perform the optimization in stages. 
     * Note that this value is in clocks, so
     * in this example we are setting the time
     * limit to 3 seconds                     */
    gsod settings;
    bzero(&settings,sizeof(gsod));
    settings.t_limit = 3 * CLOCKS_PER_SEC;

    optimize(graph,&settings);

    if (settings.end_state = OES_time_limit) {
        /* the call to optimize exited due to hitting
         * the time limit, not because of success or
         * divergence.                              */
    }

```
```
    /* another way to check for divergence */
    gsod settings;
    bzero(&settings,sizeof(gsod));

    optimize(graph,&settings);

    if (settings.end_state = OES_failure) {
        /* the system has diverged */
    }
    else 
        /* settings.end_state == OES_success */

```
```
    /* print detailed information about the
     * optimization function                */
    gsod settings;
    bzero(&settings,sizeof(gsod));

    optimize(graph,&settings);

    print_graph_slam_optimization_data(settings);

```

### Fixing Nodes

If your robot has access to some kind of GPS or local positioning system,
you can optionally fix nodes instead of simply assigning a high weight
to their edges via the information matrix.

#### A Warning

Because we are fixing nodes in a nonlinear system and applying a linear
solver, we introduce the possibility that we will never find a local
minimum for our system. Because of this possiblity, it is recommended
to do some careful testing with this feature, as the optimize function
may never return... You can also mitigate this by setting a step or time
limit as described above. The nodes you get back should be useable, they
will simply be osillating between two local minima.

#### Usage

Usage of this feature is pretty simple:
```
    /* the node we want to fix */
    pgn_handle A;

    fix_node(A,graph);

```

There is also this if you want to use it:
```
    /* remove the fixed status from a node */
    unfix_node(A,graph);
```

## Closing Thoughts

I personally think graph slam is a nice algorithm, and am
impressed with the people who were able to construct such a useful
mathematical model out of the robotics slam problem. 
   
I would like to leave a shoutout to Cyrill Stachniss and Sebastian
Thrun, who are great researchers and teachers whose material I benefitted
from greatly during this project.   

I would also like to thank Ump45(404) and Rinne Ohara for the indirect
motivation on this project. \:wa:waq

## TODO

 * Integrate the system with a full front end and run a full test   
 * clean up the test file   
 * add more explaination in the solver files   

## Sources

```
A great reference on many modern robotics algorithms.

Thrun, Sebastian. “Probabilistic robotics.” Commun. ACM 45 (2002): 52-57.
```
```
Useful paper and introduction to graph slam. Most of the variable names and mathematics are derived from this.

Grisetti, Giorgio & Kümmerle, Rainer & Stachniss, Cyrill & Burgard, Wolfram. (2010). A tutorial on graph-based SLAM. IEEE Transactions on Intelligent Transportation Systems Magazine. 2. 31-43. 10.1109/MITS.2010.939925. 
```
```
Source used in the construction of cholesky decomposition solver.

Stewart, G. W.. “Building an Old-Fashioned Sparse Solver.” (2003).
```
```
Useful lecture on the graph slam algorithm. Method of fixing nodes was derived from here.

Stachniss, C. (2020) Graph-based SLAM using Pose Graphs, YouTube. Available at: https://www.youtube.com/watch?v=uHbRKvD8TWg (Accessed: 28 January 2025). 
```
```
Helpful for double-checking my approach of fixing nodes.

Eustice, Ryan M., Hanumant Singh and John J. Leonard. “Exactly Sparse Delayed-State Filters for View-Based SLAM.” IEEE Transactions on Robotics 22 (2006): 1100-1114.
```
