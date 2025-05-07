# Computational Geometry

This is an implementation (and some halfway implementations) of some of the
data structures and algorithms from "Computational Geometry", by Berg et al.

It includes basic data structures, points, vectors, segments, lines, rays, and
polygons. A general data structures for planar subdivision, double connected
edge list (DCEL), and a half implementation of a trapezoidal map.

The library come with geometrical constant operations, such as intersection,
left of, right of, and it implemente the next algorithms

 - convex hull (Graham scan, Jarvis march, and quick hull),
 - triangulation (includes make monotone polygons, and triangulate those),
 - voronoi (incremental voronoi).

## Examples

You can found the project that inspired this implementation under
[examples](/examples).

### Known issues

If two or more edges intersect on the same point, it will surely hit different
points of intersection, or worse, miss the intersection. This need to be
addressed with better implementations of geometric primitives and careful
post-processing.

---

© 2025 José Joaquín Zubieta Rico

