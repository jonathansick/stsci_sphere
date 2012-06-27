User documentation
==================

.. currentmodule:: sphere

The `sphere` library is a pure Python package for handling spherical
polygons that represent arbitrary regions of the sky.

Requirements
------------

- Python 2.7

- Numpy 1.4 or later

- PyFITS

- PyWCS

- STWCS

Coordinate representation
-------------------------

Coordinates in world space are traditionally represented by right
ascension and declination (*ra* and *dec*), or longitude and latitude.
While these representations are convenient, they have discontinuities
at the poles, making operations on them trickier at arbitrary
locations on the sky sphere.  Therefore, all internal operations of
this library are done in 3D vector space, where coordinates are
represented as (*x*, *y*, *z*) vectors.  The `sphere.vector` module
contains functions to convert between (*ra*, *dec*) and (*x*, *y*,
*z*) representations.

While any (*x*, *y*, *z*) triple represents a vector and therefore a
location on the sky sphere, a distinction must be made between
normalized coordinates fall exactly on the unit sphere, and
unnormalized coordinates which do not.  A normalized coordinate is
defined as a vector whose length is 1, i.e.:

.. math::

    \sqrt{x^2 + y^2 + z^2} = 1

To prevent unnecessary recomputation, many methods in this library
assume that the vectors passed in are already normalized.  If this is
not the case, `sphere.vector.normalize_vector` can be used to
normalize an array of vectors.

The library allows the user to work in either degrees or radians.  All
methods that require or return an angular value have a `degrees`
keyword argument.  When `degrees` is `True`, these measurements are in
degrees, otherwise they are in radians.

Spherical polygons
------------------

Spherical polygons are arbitrary areas on the sky sphere enclosed by
great circle arcs.  They are represented by the
`~sphere.polygon.SphericalPolygon` class.

Representation
``````````````

The points defining the polygon are available from the
`~polygon.SphericalPolygon.points` property.  It is a Nx3 array where
each row is an (*x*, *y*, *z*) vector, normalized.  The polygon points
are explicitly closed, i.e., the first and last points are the same.

Where is the inside?
^^^^^^^^^^^^^^^^^^^^

The edges of a polygon serve to separate the “inside” from the
“outside” area.  On a traditional 2D planar surface, the “inside” is
defined as the finite area and the “outside” is the infinite area.
However, since the surface of a sphere is cyclical, i.e., it wraps
around on itself, the a spherical polygon actually defines two finite
areas.  To specify which should be considered the “inside” vs. the
“outside”, the definition of the polygon also has an “inside point”
which is just any point that should be considered inside of the
polygon.

In the following image, the inside point (marked with the red dot)
declares that the area of the polygon is the green region, and not the
white region.

.. image:: inside.png

The inside point of the the polygon can be obtained from the
`~polygon.SphericalPolygon.inside` property.

Cut lines
^^^^^^^^^

If the polygon represents two disjoint areas or the polygon has holes,
those areas will be connected by cut lines.  The following image shows
a polygon made from the union of a number of cone areas which has both
a hole and a disjoint region connected by cut lines.

.. image:: cutlines.png

Creating spherical polygons
```````````````````````````

.. currentmodule:: sphere.polygon

`SphericalPolygon` objects have 4 different constructors:

  - `SphericalPolygon`: Takes an array of (*x*, *y*, *z*)
    points and an inside point.

  - `SphericalPolygon.from_radec`: Takes an array of (*ra*, *dec*)
    points and an inside point.

  - `SphericalPolygon.from_cone`: Creates a polygon from a cone on the
    sky shere.  Takes (*ra*, *dec*, *radius*).

  - `SphericalPolygon.from_wcs`: Creates a polygon from the footprint
    of a FITS image using its WCS header keywords.  Takes a FITS
    filename or a `pyfits.Header` object.

Operations on Spherical Polygons
````````````````````````````````

Once one has a `SphericalPolygon` object, there are a number of
operations available:

  - `~SphericalPolygon.contains_point`: Determines if the given point is inside the polygon.

  - `~SphericalPolygon.intersects_poly`: Determines if one polygon intersects with another.

  - `~SphericalPolygon.area`: Determine the area of a polygon.

  - `~SphericalPolygon.union` and `~SphericalPolygon.multi_union`:
    Return a new polygon that is the union of two or more polygons.

  - `~SphericalPolygon.intersection` and
    `~SphericalPolygon.multi_intersection`: Return a new polygon that
    is the intersection of two or more polygons.

  - `~SphericalPolygon.overlap`: Determine how much a given polygon
    overlaps another.

  - `~SphericalPolygon.to_radec`: Convert (*x*, *y*, *z*) points in the
    polygon to (*ra*, *dec*) points.

  - `~SphericalPolygon.same_points_as`: Determines if one polygon has the
    same points as another. When only sorted unique points are considered
    (default behavior), polygons with same points might not be the same
    polygons because the order of the points matter.

  - `~SphericalPolygon.draw`: Plots the polygon using matplotlib’s
    Basemap toolkit.  This feature is rather bare and intended
    primarily for debugging purposes.

Great circle arcs
-----------------

.. currentmodule:: sphere.great_circle_arc

As seen above, great circle arcs are used to define the edges of the
polygon.  The `sphere.great_circle_arc` module contains a number of
functions that are useful for dealing with them.

- `length`: Returns the angular distance between two points on the sphere.

- `intersection`: Returns the intersection point between two great
  circle arcs.

- `intersects`: Determines if two great circle arcs intersect.

- `angle`: Calculate the angle between two great circle arcs.

- `midpoint`: Calculate the midpoint along a great circle arc.

Skylines
--------

Skylines are designed to capture and manipulate HST WCS image information as
spherical polygons. They are represented by the `~sphere.skyline.SkyLine` class,
which is an extension of `~sphere.polygon.SphericalPolygon` class.

Representation
``````````````
Each skyline has a list of members, `~sphere.skyline.SkyLine.members`, and a
composite spherical polygon, `~sphere.skyline.SkyLine.polygon`, defined by those
members. The polygon has all the functionalities of
`~sphere.polygon.SphericalPolygon`.

What is a skyline member?
^^^^^^^^^^^^^^^^^^^^^^^^^

Each member in `~sphere.skyline.SkyLine.members` belongs to the
`~sphere.skyline.SkyLineMember` class, which contains image name (with path if
given), science extension, and WCS object and polygon of that extension.

For example, an ACS/WFC full-frame image would give 2 members, one from EXT 1
and another from EXT 4.

Creating skylines
`````````````````

`~sphere.skyline.SkyLine` constructor takes an image name and an optional
`extname` keyword, which defaults to "SCI". To create skyline from
single-extension FITS, change `extname` to "PRIMARY".

If `None` is given instead of image name, an empty skyline is created with no
member and an empty spherical polygon.

Operations on skylines
``````````````````````

`~sphere.skyline.SkyLine` has direct access to most of the
`~sphere.polygon.SphericalPolygon` properties and methods *except* for the
following (which are still accessible indirectly via
`~sphere.skyline.SkyLine.polygon`):

  - `~sphere.polygon.SphericalPolygon.from_radec`
  - `~sphere.polygon.SphericalPolygon.from_cone`
  - `~sphere.polygon.SphericalPolygon.from_wcs`
  - `~sphere.polygon.SphericalPolygon.multi_union`
  - `~sphere.polygon.SphericalPolygon.multi_intersection`

In addition, `~sphere.skyline.SkyLine` also has these operations available:

  - `~sphere.skyline.SkyLine.to_wcs`: Return a composite HST WCS object defined
    by all the members.

  - `~sphere.skyline.SkyLine.add_image`: Return a new skyline that is the union
    of two skylines. This should be used, *not* `SkyLine.union` (which is
    actually `~sphere.polygon.SphericalPolygon.union`) that will not include
    members.

  - `~sphere.skyline.SkyLine.find_intersection`: Return a new skyline that is
    the intersection of two skylines. This should be used, *not*
    `SkyLine.intersection` (which is actually
    `~sphere.polygon.SphericalPolygon.intersection`) that will not include
    members.

  - `~sphere.skyline.SkyLine.find_max_overlap` and
    `~sphere.skyline.SkyLine.max_overlap_pair`: Return a pair of skylines that
    overlap the most from a given list of skylines.

  - `~sphere.skyline.SkyLine.mosaic`: Return a new skyline that is a mosaic of
    given skylines that overlap, a list of image names of the skylines used, and
    a list of image names of the excluded skylines. A pair of skylines with the
    most overlap is used as a starting point. Then a skyline that overlaps the
    most with the mosaic is used, and so forth until no overlapping skyline is
    found.