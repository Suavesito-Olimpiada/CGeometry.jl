# Example

## Art gallery

In `art_gallery.jl`, the problem of cameras in an art gallery is solved. The big
colored (red, green or blue) markers at the end of the gifs are the positions
of the cameras inside the polygon.

First, a radially ordered polygon is created for test. It is partitioned in
$y$-monotone polygons, which are then triangulated.

Two examples of this can be found in `triangulation40.gif` and
`triangulation100.gif`. A third examples is shown in
`triangulationvoronoi.gif`, where a voronoi planar subdivision is taken and
triangulated (the cells are monotone in every direction).

### Make monotone and triangulation

![Triangulation (40)](./triangulation40.gif)

![Triangulation (100)](./triangulation100.gif)

### Triangulation

![Triangulation voronoi](./triangulationvoronoi.gif)

## Convex hull

Three algorithms covered in the book were implemented, with several examples
each. These were obtained with code similar to what is inside of `convexhull.jl`.

### Jarvis March

![Jarvis march (16)](./jarvismarch16.gif)

![Jarvis march, in grid (16)](./jarvismarch16_squared.gif)

![Jarvis march (64)](./jarvismarch64.gif)

### Graham scan

![Graham scan 64](./grahamscan64.gif)

![Graham scan, in circle (64)](./grahamscan64_round.gif)

### Quick hull

![Quick hull, in grid (16)](./quickhull16_squared.gif)

![Quick hull (64)](./quickhull64.gif)

![Quick hull, in circle (64)](./quickhull64_round.gif)

## Voronoi

The incremental voronoi algorithm is implemented, which is quadratic in
complexity, but easier to program. Several examples to show it in the workings,
all implemented in `voronoi.jl`.

![Voronoi (16)](./voronoi16.gif)

![Voronoi, in circle (16)](./voronoi16_round.gif)

![Voronoi, in grid (16)](./voronoi16_square.gif)

![Voronoi (32)](./voronoi32.gif)
