# -*- coding: utf-8 -*-
"""
ObsPy utilities for relationships between points on a sphere or ellipsoid.

the `obspy.geodetics.geodesic` module includes functions to find the distance
and direction between pairs of points on an ellipsoid (the inverse geodetic
ptoblem) and the location of a point given a starting point, distance and
direction (the direct geodetic problem).

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from collections import namedtuple

from geographiclib.geodesic import Geodesic


PointDistance = namedtuple('PointDistance', ['latitude', 'longitude', 'degree',
                                             'meter', 'azimuth',
                                             'backazimuth'])


def direct_geodetic_degree(latitude, longitude, azimuth, distance_in_degrees,
                           semi_major_axis, flattening):
    """
    Computes the latitude and longitude of a geographic point on an ellipsoid
    given and a starting location, an azimuth and the angular distance.

    :param latitude: Latitude of point A in degrees (positive for northern,
        negative for southern hemisphere)
    :param longitude: Longitude of point A in degrees (positive for eastern,
        negative for western hemisphere)
    :param azimuth: Azimuth from point A to point B, clockwise from north in
        degrees.
    :param distance_in_degrees: Angular distance from point A to point B
        measured at the center of the ellipsoid, in degrees.
    :param semi_major_axis: Radius of Earth in m.
    :param flattening: Flattening of Earth.
    :return: PointDistance named tuple including the location of point B.
    """
    ellipsoid = Geodesic(a=semi_major_axis, f=flattening)
    g = ellipsoid.ArcDirect(latitude, longitude, azimuth,
                            distance_in_degrees)
    result = PointDistance(latitude=g['lat2'], longitude=g['lon2'],
                           degree=g['a12'], meter=g['s12'],
                           azimuth=_angle_wrap(g['azi1']),
                           backazimuth=_angle_wrap(g['azi2']+180.0))
    return result


def direct_geodetic_meter(latitude, longitude, azimuth, distance_in_meters,
                          semi_major_axis, flattening):
    """
    Computes the latitude and longitude of a geographic point on an ellipsoid
    given and a starting location, an azimuth and the distance on the surface
    of the ellipsoid.

    :param latitude: Latitude of point A in degrees (positive for northern,
        negative for southern hemisphere)
    :param longitude: Longitude of point A in degrees (positive for eastern,
        negative for western hemisphere)
    :param azimuth: Azimuth from point A to point B, clockwise from north in
        degrees.
    :param distance_in_meters: distance from point A to point B measured on the
        surface of the ellipsoid, in m.
    :param semi_major_axis: Radius of Earth in m.
    :param flattening: Flattening of Earth.
    :return: PointDistance named tuple including the location of point B.
    """
    ellipsoid = Geodesic(a=semi_major_axis, f=flattening)
    g = ellipsoid.Direct(latitude, longitude, azimuth, distance_in_meters)
    result = PointDistance(latitude=g['lat2'], longitude=g['lon2'],
                           degree=g['a12'], meter=g['s12'],
                           azimuth=_angle_wrap(g['azi1']),
                           backazimuth=_angle_wrap(g['azi2']+180.0))
    return result


def inverse_geodetic(latitude_1, longitude_1, latitude_2, longitude_2,
                     semi_major_axis, flattening):
    """
    Computes the distance between two geographic points on an ellipsoid and
    the forward and backward azimuths between these points.

    :param latitude_1: Latitude of point A in degrees (positive for northern,
        negative for southern hemisphere)
    :param longitude_1: Longitude of point A in degrees (positive for eastern,
        negative for western hemisphere)
    :param latitude_2: Latitude of point B in degrees (positive for northern,
        negative for southern hemisphere)
    :param longitude_2: Longitude of point B in degrees (positive for eastern,
        negative for western hemisphere)
    :param semi_major_axis: Radius of Earth in m.
    :param flattening: Flattening of Earth.
    :return: PointDistance named tuple.
    """
    ellipsoid = Geodesic(a=semi_major_axis, f=flattening)
    g = ellipsoid.Inverse(latitude_1, longitude_1, latitude_2, longitude_2)
    result = PointDistance(latitude=g['lat2'], longitude=g['lon2'],
                           degree=g['a12'], meter=g['s12'],
                           azimuth=_angle_wrap(g['azi1']),
                           backazimuth=_angle_wrap(g['azi2']+180.0))
    return result


def _angle_wrap(angle):
    if angle < 0.0:
        angle = angle + 360.0
    if angle > 360.0:
        angle = angle - 360.0
    return angle
