# -*- coding: utf-8 -*-

# Copyright (C) 2011 Association of Universities for Research in
# Astronomy (AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#     1. Redistributions of source code must retain the above
#       copyright notice, this list of conditions and the following
#       disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials
#       provided with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.

"""
This module provides support for working with footprints
on the sky. Primary use case would use the following
generalized steps:

    #. Initialize `SkyLine` objects for each input image.
       This object would be the union of all the input
       image's individual chips WCS footprints.

    #. Determine overlap between all images. The
       determination would employ a recursive operation
       to return the extended list of all overlap values
       computed as [img1 vs [img2,img3,...,imgN],img2 vs
       [img3,...,imgN],...]

    #. Select the pair with the largest overlap, or the
       pair which produces the largest overlap with the
       first input image. This defines the initial
       reference `SkyLine` object.

    #. Perform some operation on the 2 images: for example,
       match sky in intersecting regions, or aligning
       second image with the first (reference) image.

    #. Update the second image, either apply the sky value
       or correct the WCS, then generate a new `SkyLine`
       object for that image.

    #. Create a new reference `SkyLine` object as the union
       of the initial reference object and the newly
       updated `SkyLine` object.

    #. Repeat Steps 2-6 for all remaining input images.

This process will work reasonably fast as most operations
are performed using the `SkyLine` objects and WCS information
solely, not image data itself.

"""
from __future__ import division, print_function, absolute_import

# STDLIB
from copy import copy, deepcopy

# THIRD-PARTY
import pyfits
from stwcs import wcsutil
from stwcs.distortion.utils import output_wcs

# LOCAL
from .polygon import SphericalPolygon

# DEBUG
SKYLINE_DEBUG = True

__all__ = ['SkyLineMember', 'SkyLine']
__version__ = '0.4a'
__vdate__ = '11-Jul-2012'

class SkyLineMember(object):
    """
    Container for `SkyLine` members with these attributes:
    
        * `fname`: Image name (with path if given)
        * `ext`: Tuple of extensions read
        * `wcs`: `HSTWCS` object the composite data
        * `polygon`: `~sphere.polygon.SphericalPolygon` object of the composite data

    """
    def __init__(self, fname, extname):
        """
        Parameters
        ----------
        fname : str
            FITS image.

        extname : str
            EXTNAME to use. SCI is recommended for normal
            HST images. PRIMARY if image is single ext.

        """
        extname = extname.upper()
        ext_list = []
        wcs_list = []
        
        with pyfits.open(fname) as pf:
            for i,ext in enumerate(pf):
                if ext.name.upper() == extname:
                    ext_list.append(i)
                    wcs_list.append(wcsutil.HSTWCS(fname, ext=i))

        # By combining WCS first before polygon, will remove chip gaps

        n_wcs = len(wcs_list)
        if n_wcs > 1:
            self._wcs = output_wcs(wcs_list)
        elif n_wcs == 1:
            self._wcs = wcs_list[0]
        else:
            raise ValueError('%s has no WCS' % fname)

        self._fname = fname
        self._ext = tuple(ext_list)
        self._polygon = SphericalPolygon.from_wcs(self.wcs)

    def __repr__(self):
        return '%s(%r, %r, %r, %r)' % (self.__class__.__name__, self.fname,
                                       self.ext, self.wcs, self.polygon)

    @property
    def fname(self):
        return self._fname

    @property
    def ext(self):
        return self._ext

    @property
    def wcs(self):
        return self._wcs

    @property
    def polygon(self):
        return self._polygon

class SkyLine(object):
    """
    Manage outlines on the sky.

    Each `SkyLine` has a list of `~SkyLine.members` and
    a composite `~SkyLine.polygon` with all the
    functionalities of `~sphere.polygon.SphericalPolygon`.

    """
    def __init__(self, fname, extname='SCI'):
        """
        Parameters
        ----------
        fname : str
            FITS image. `None` to create empty `SkyLine`.

        extname : str
            EXTNAME to use. SCI is recommended for normal
            HST images. PRIMARY if image is single ext.

        """      
        # Convert SCI data to SkyLineMember
        if fname is not None:
            self.members = [SkyLineMember(fname, extname)]
        else:
            self.members = []

        # Put mosaic of all the chips in SkyLine
        n = len(self.members)
        if n == 0:
            self.polygon = SphericalPolygon([])
        elif n == 1:
            self.polygon = copy(self.members[0].polygon)
        else:
            raise ValueError('%s cannot initialize polygon with '
                             'multiple members' % self.__class__.__name__)

    def __getattr__(self, what):
        """Control attribute access to `~sphere.polygon.SphericalPolygon`."""
        if what in ('from_radec', 'from_cone', 'from_wcs',
                    'multi_union', 'multi_intersection',
                    '_find_new_inside',):
            raise AttributeError('\'%s\' object has no attribute \'%s\'' %
                                 (self.__class__.__name__, what))
        else:
            return getattr(self.polygon, what)

    def __copy__(self):
        return deepcopy(self)
    
    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__,
                               self.polygon, self.members)

    @property
    def polygon(self):
        """
        `~sphere.polygon.SphericalPolygon` portion of `SkyLine`
        that contains the composite skyline from `members`
        belonging to *self*.

        """
        return self._polygon

    @polygon.setter
    def polygon(self, value):
        assert isinstance(value, SphericalPolygon)
        self._polygon = copy(value)  # Deep copy

    @property
    def members(self):
        """
        List of `SkyLineMember` objects that belong to *self*.
        Duplicate members are discarded. Members are kept in
        the order of their additions to *self*.

        """
        return self._members

    @members.setter
    def members(self, values):
        self._members = []

        # Not using set to preserve order
        for v in values:
            # Report corrupted members list instead of skipping
            assert isinstance(v, SkyLineMember)

            if v not in self._members:
                self._members.append(v)

    def to_wcs(self):
        """
        Combine `HSTWCS` objects from all `members` and return
        a new `HSTWCS` object. If no `members`, return `None`.

        .. warning:: This cannot return WCS of intersection.

        """
        wcs_list = []
        
        for m in self.members:
            for i in m.ext:
                wcs_list.append(wcsutil.HSTWCS(m.fname, ext=i))

        if len(wcs_list) > 0:
            wcs = output_wcs(wcs_list)
        else:
            wcs = None

        return wcs

    def _rough_id(self):
        """Filename of first member."""
        if len(self.members) > 0:
            return self.members[0].fname
        else:
            return None

    def _draw_members(self, map, **kwargs):
        """
        Draw individual extensions in members.
        Useful for debugging.

        Parameters
        ----------
        map : Basemap axes object

        **kwargs : Any plot arguments to pass to basemap

        """
        for m in self.members:
            for i in m.ext:
                poly = SphericalPolygon.from_wcs(wcsutil.HSTWCS(m.fname, ext=i))
                poly.draw(map, **kwargs)

    def _find_members(self, given_members):
        """
        Find `SkyLineMember` in *given_members* that is in
        *self*. This is used for intersection.

        Parameters
        ----------
        self : obj
            `SkyLine` instance.

        given_members : list
            List of `SkyLineMember` to consider.

        Returns
        -------
        new_members : list
            List of `SkyLineMember` belonging to *self*.

        """
        if len(self.points) > 0:
            out_mem = [m for m in given_members if
                       self.intersects_poly(m.polygon)]
        else:
            out_mem = []
        return out_mem

    def add_image(self, other):
        """
        Return a new `SkyLine` that is the union of *self*
        and *other*.

        .. warning:: `SkyLine.union` only returns `polygon`
            without `members`.

        Parameters
        ----------
        other : `SkyLine` object

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s3 = s1.add_image(s2)

        """
        newcls = self.__class__(None)
        newcls.polygon = self.union(other)
        newcls.members = self.members + other.members
        return newcls

    def find_intersection(self, other):
        """
        Return a new `SkyLine` that is the intersection of
        *self* and *other*.

        .. warning:: `SkyLine.intersection` only returns
            `polygon` without `members`.

        Parameters
        ----------
        other : `SkyLine` object

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s3 = s1.find_intersection(s2)

        """
        newcls = self.__class__(None)
        newcls.polygon = self.intersection(other)
        newcls.members = newcls._find_members(self.members + other.members)
        return newcls

    def find_max_overlap(self, skylines):
        """
        Find `SkyLine` from a list of *skylines* that overlaps
        the most with *self*.

        Parameters
        ----------
        skylines : list
            A list of `SkyLine` instances.

        Returns
        -------
        max_skyline : `SkyLine` instance or `None`
            `SkyLine` that overlaps the most or `None` if no
            overlap found. This is *not* a copy.

        max_overlap_area : float
            Area of intersection.
        
        """
        max_skyline = None
        max_overlap_area = 0.0

        for next_s in skylines:
            try:
                overlap_area = self.intersection(next_s).area()
            except (ValueError, AssertionError):
                if SKYLINE_DEBUG:
                    print('WARNING: Intersection failed for %s and %s. '
                          'Ignoring %s...' % (self._rough_id(),
                                              next_s._rough_id(),
                                              next_s._rough_id()))
                    overlap_area = 0.0
                else:
                    raise

            if overlap_area > max_overlap_area:
                max_overlap_area = overlap_area
                max_skyline = next_s

        return max_skyline, max_overlap_area

    @staticmethod
    def max_overlap_pair(skylines):
        """
        Find a pair of skylines with maximum overlap.

        Parameters
        ----------
        skylines : list
            A list of `SkyLine` instances.

        Returns
        -------
        max_pair : tuple
            Pair of `SkyLine` objects with max overlap
            among given *skylines*. If no overlap found,
            return `None`. These are *not* copies.

        """       
        max_pair = None
        max_overlap_area = 0.0
    
        for i in xrange(len(skylines) - 1):
            curr_s = skylines[i]
            next_s, i_area = curr_s.find_max_overlap(skylines[i+1:])

            if i_area > max_overlap_area:
                max_overlap_area = i_area
                max_pair = (curr_s, next_s)

        return max_pair

    @classmethod
    def mosaic(cls, skylines, verbose=True):
        """
        Mosaic all overlapping *skylines*.

        A pair of skylines with the most overlap is used as
        a starting point. Then a skyline that overlaps the
        most with the mosaic is used, and so forth until no
        overlapping skyline is found.

        Parameters
        ----------
        skylines : list
            A list of `SkyLine` objects.

        verbose : bool
            Print info to screen.

        Returns
        -------
        mosaic : `SkyLine` instance or `None`
            Union of all overlapping *skylines*, or `None` if
            no overlap found.

        included : list
            List of image names added to mosaic in the order
            of addition.

        excluded : list
            List of image names excluded because they do not
            overlap with mosaic.

        """
        out_order = []
        excluded  = []

        if verbose:
            print('***** SKYLINE MOSAIC *****')
        
        starting_pair = cls.max_overlap_pair(skylines)
        if starting_pair is None:
            if verbose:
                print('    Cannot find any overlapping skylines. Aborting...')
            return starting_pair, out_order, excluded

        remaining = list(skylines)

        s1, s2 = starting_pair
        if verbose:
            print('    Starting pair: %s, %s' %
                  (s1._rough_id(), s2._rough_id()))

        mosaic = s1.add_image(s2)
        out_order = [s1._rough_id(), s2._rough_id()]
        remaining.remove(s1)
        remaining.remove(s2)

        while len(remaining) > 0:
            next_skyline, i_area = mosaic.find_max_overlap(remaining)

            if next_skyline is None:
                for r in remaining:
                    if verbose:
                        print('    No overlap: Excluding %s...' % r._rough_id())
                    excluded.append(r._rough_id())
                break

            try:
                new_mos = mosaic.add_image(next_skyline)
            except (ValueError, AssertionError):
                if SKYLINE_DEBUG:
                    print('WARNING: Cannot add %s to mosaic. Skipping it...' %
                          next_skyline._rough_id())
                    excluded.append(next_skyline._rough_id())
                else:
                    raise
            else:
                print('    Adding %s to mosaic...' % next_skyline._rough_id())
                mosaic = new_mos
                out_order.append(next_skyline._rough_id())
            finally:
                remaining.remove(next_skyline)

        return mosaic, out_order, excluded

    @classmethod
    def _find_frosty(cls, show=False):
        s1 = SphericalPolygon.from_cone(0, 35, 8)
        s2 = SphericalPolygon.from_cone(0, 20, 13)
        s3 = SphericalPolygon.from_cone(0, 0, 20)
        ss = SphericalPolygon.multi_union([s1,s2,s3])
        frosty = cls(None)
        frosty.polygon = ss
        if show:
            from mpl_toolkits.basemap import Basemap
            map = Basemap()
            frosty.draw(map)
            print('Frosty the Snowman says Brrr so cold...')
        return frosty
