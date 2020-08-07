# Cpplap
Originally cpplap (c++ lines and planes) has been developed for a small [raytracing project](https://github.com/hsandmeyer/dmd_traycing).
However it might become useful for anybody who needs a simple, one header library to compute interactions between vectors, lines and planes in 3D-space.
For instance to compute the intersection point between a line and a plane, use
```cpp
cpplap::Plane<double> plane(3, 1, 5, 3); //3x + 2y + 5z = 3
cpplap::Line<double> line({1, 2, 3}, {3, 1, 2}); //{1, 2, 3} + t * {3, 1, 2}
cpplap::Vect<double> intersection = plane.intersectionPointWith(line);
```
To compute the distance vector between a plane and a point use
```cpp
cpp::Vect<double> dist = plane.distanceTo({2, 2, 2});
```
For a list of all possible computations, have a look at [example.cpp](test/example.cpp)
