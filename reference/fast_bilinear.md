# A fast bilinear interpolation function

It assumes equal grid intervals, thus can find the correct position in
\\O(1)\\ time.

## Usage

``` r
fast_bilinear(x, y, z, x0, y0)
```

## Arguments

- x:

  a vector containing the `x` coordinates of the rectangular data grid.

- y:

  a vector containing the `y` coordinates of the rectangular data grid.

- z:

  a matrix containing the `z[i,j]` data values for the grid points
  (`x[i]`,`y[j]`).

- x0:

  `x` value used to interpolate at.

- y0:

  `y` value used to interpolate at.

## Value

A list which contains only one element, `z`.

## Details

The following is from the documentation of
[`akima::bilinear()`](https://rdrr.io/pkg/akima/man/bilinear.html).

This is an implementation of a bilinear interpolating function. For a
point (x0,y0) contained in a rectangle (x1,y1),(x2,y1), (x2,y2),(x1,y2)
and x1\<x2, y1\<y2, the first step is to get z() at locations (x0,y1)
and (x0,y2) as convex linear combinations
z(x0,y\*)=a*z(x1,y*)+(1-a)*z(x2,y*) where a=(x2-x1)/(x0-x1) for
y\*=y1,y2. In a second step z(x0,y0) is calculated as convex linear
combination between z(x0,y1) and z(x0,y2) as
z(x0,y1)=b\*z(x0,y1)+(1-b)\*z(x0,y2) where b=(y2-y1)/(y0-y1).

Finally, z(x0,y0) is a convex linear combination of the z values at the
corners of the containing rectangle with weights according to the
distance from (x0,y0) to these corners.
