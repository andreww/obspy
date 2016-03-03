# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import unittest

from obspy.geodetics import (direct_geodetic_meter, inverse_geodetic)


class GeodesicTestCase(unittest.TestCase):
    """
    Test suite for obspy.geodetics.geodesic
    """
    def test_inverse_australia(self):
        """
        Test inverse_geodetic() method with test data from Geocentric Datum of
        Australia. (see http://www.icsm.gov.au/gda/gdatm/gdav2.3.pdf)
        """
        # test data:
        # Point 1: Flinders Peak, Point 2: Buninyong
        lat1 = -(37 + (57 / 60.) + (3.72030 / 3600.))
        lon1 = 144 + (25 / 60.) + (29.52440 / 3600.)
        lat2 = -(37 + (39 / 60.) + (10.15610 / 3600.))
        lon2 = 143 + (55 / 60.) + (35.38390 / 3600.)
        dist = 54972.271
        alpha12 = 306 + (52 / 60.) + (5.37 / 3600.)
        alpha21 = 127 + (10 / 60.) + (25.07 / 3600.)

        WGS84_A = 6378137.0
        WGS84_F = 1 / 298.257223563

        # calculate result
        res = inverse_geodetic(lat1, lon1, lat2, lon2, WGS84_A, WGS84_F)

        # calculate deviations from test data
        dist_err_rel = abs(dist - res.meter) / dist
        alpha12_err = abs(alpha12 - res.azimuth)
        alpha21_err = abs(alpha21 - res.backazimuth)

        self.assertEqual(dist_err_rel < 1.0e-5, True)
        self.assertEqual(alpha12_err < 1.0e-5, True)
        self.assertEqual(alpha21_err < 1.0e-5, True)

    def test_direct_meter_australia(self):
        """
        Test direct_geodetici_meter() with test data from Geocentric Datum of
        Australia. (see http://www.icsm.gov.au/gda/gdatm/gdav2.3.pdf)
        """
        # test data:
        # Point 1: Flinders Peak, Point 2: Buninyong
        lat1 = -(37 + (57 / 60.) + (3.72030 / 3600.))
        lon1 = 144 + (25 / 60.) + (29.52440 / 3600.)
        lat2 = -(37 + (39 / 60.) + (10.15610 / 3600.))
        lon2 = 143 + (55 / 60.) + (35.38390 / 3600.)
        dist = 54972.271
        alpha12 = 306 + (52 / 60.) + (5.37 / 3600.)
        alpha21 = 127 + (10 / 60.) + (25.07 / 3600.)

        WGS84_A = 6378137.0
        WGS84_F = 1 / 298.257223563

        # calculate result
        res = direct_geodetic_meter(lat1, lon1, alpha12, dist, WGS84_A,
                                    WGS84_F)

        # calculate deviations from test data
        lat2_err = abs(lat2 - res.latitude)
        lon2_err = abs(lon2 - res.longitude)
        alpha21_err = abs(alpha21 - res.backazimuth)

        self.assertEqual(lat2_err < 1.0e-5, True)
        self.assertEqual(lon2_err < 1.0e-5, True)
        self.assertEqual(alpha21_err < 1.0e-5, True)

    def test_direct_meter_2_australia(self):
        """
        Test direct_geodetici_meter() with test data from Geocentric Datum of
        Australia. (see http://www.icsm.gov.au/gda/gdatm/gdav2.3.pdf)
        """
        # test data:
        # Point 1: Flinders Peak, Point 2: Buninyong
        lat1 = -(37 + (57 / 60.) + (3.72030 / 3600.))
        lon1 = 144 + (25 / 60.) + (29.52440 / 3600.)
        lat2 = -(37 + (39 / 60.) + (10.15610 / 3600.))
        lon2 = 143 + (55 / 60.) + (35.38390 / 3600.)
        dist = 54972.271
        alpha12 = 306 + (52 / 60.) + (5.37 / 3600.)
        alpha21 = 127 + (10 / 60.) + (25.07 / 3600.)

        WGS84_A = 6378137.0
        WGS84_F = 1 / 298.257223563

        # calculate result
        res = direct_geodetic_meter(lat2, lon2, alpha21, dist, WGS84_A,
                                    WGS84_F)

        # calculate deviations from test data
        lat1_err = abs(lat1 - res.latitude)
        lon1_err = abs(lon1 - res.longitude)
        alpha12_err = abs(alpha12 - res.backazimuth)

        self.assertEqual(lat1_err < 1.0e-5, True)
        self.assertEqual(lon1_err < 1.0e-5, True)
        self.assertEqual(alpha12_err < 1.0e-5, True)


def suite():
    return unittest.makeSuite(GeodesicTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
