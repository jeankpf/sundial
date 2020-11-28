from math import radians, degrees, tan, pi, acos
import orekit

vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()

from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.utils import IERSConventions, Constants, PVCoordinatesProvider
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.time import AbsoluteDate, TimeScalesFactory

def utc2local_time(utc_date, time_zone, summer_time):
    """
    Converts the date from UTC to the local time.

    Parameters
    ----------
    utc_date : Absolute Date
        Date in UTC to be converted.
    time_zone : int
        Time zone of the local area.
    summer_time : bool
        True if it is summer time, False otherwise.

    Returns
    -------
    Absolute Date
        Converted date in local time.

    """
    if summer_time:
        time_zone += 1
    return utc_date.shiftedBy(3600.0*time_zone)

def local_time2utc(local_time, time_zone, summer_time):
    """
    Converts the date from the local time to UTC.

    Parameters
    ----------
    local_date : Absolute Date
        Date in local time to be converted.
    time_zone : int
        Time zone of the local area.
    summer_time : bool
        True if it is summer time, False otherwise.

    Returns
    -------
    Absolute Date
        Converted date in UTC.

    """
    if summer_time:
        time_zone += 1
    return local_time.shiftedBy(-3600.0*time_zone)

def get_elevation(sun_pos_sundial_frame):
    """
    Returns the elevation of the sun in the sundial frame.

    Parameters
    ----------
    sun_pos_sundial_frame : Vector3D
        Position vector of the Sun in the sundial frame.

    Returns
    -------
    float
        Elevation (in radians) of the Sun in the sundial frame.

    """
    proj_sun_pos_sundial_frame = Vector3D(sun_pos_sundial_frame.getX(),
                                          sun_pos_sundial_frame.getY(), 0.0)
    cos_elevation = Vector3D.dotProduct(sun_pos_sundial_frame, proj_sun_pos_sundial_frame)/ \
    (sun_pos_sundial_frame.getNorm()*proj_sun_pos_sundial_frame.getNorm())
    return acos(cos_elevation)

def get_azimuth(sun_pos_sundial_frame):
    """
    Returns the azimuth of the sun in the sundial frame.

    Parameters
    ----------
    sun_pos_sundial_frame : Vector3D
        Position vector of the Sun in the sundial frame.

    Returns
    -------
    float
        Azimuth (in radians) of the Sun in the sundial frame.

    """
    proj_sun_pos_sundial_frame = Vector3D(sun_pos_sundial_frame.getX(),
                                          sun_pos_sundial_frame.getY(), 0.0)
    cos_azimuth = proj_sun_pos_sundial_frame.getY()/proj_sun_pos_sundial_frame.getNorm()
    return acos(cos_azimuth)


class Sundial:

    def __init__(self, lat=48.58, longi=7.75, alt=142.0, loc="Strasbourg Frame",
                 time_zone=1.0, gnomon_length=1.0):
        """
        Initiate the Sundial instance. The default is a sundial located in Strasbourg.

        Parameters
        ----------
        lat : float, optional
            Latitude of the murial sundial. The default is 48.58 (Strasbourg).
        longi : float, optional
            Longitude of the mural sundial. The default is 7.75 (Strasbourg).
        alt : float, optional
            Altitude of the mural sundial. The default is 142.0 (Strasbourg).
        loc : str, optional
            Name of the place where the murial sundial will stand. The default is "Strasbourg".
        time_zone : int, optional
            Time zone of the place where the murial sundial will stand . The default is 1 (Strasbourg).
        gnomon_length : float, optional
            Length of the gnomon of the murial Sundial. The default is 1.0.

        Returns
        -------
        None.

        """
        self.latitude = radians(lat)
        self.longitude = radians(longi)
        self.altitude = alt
        # utc = TimeScalesFactory.getUTC()
        # date = AbsoluteDate(year, month, 1, 0, 0, 0.0, utc)
        # self.date = date
        self.location = loc
        self.time_zone = time_zone
        # self.summer_time = summer_time
        self.gnomon_length = gnomon_length
        self.station_geo_point = GeodeticPoint(self.latitude, self.longitude, self.altitude)
        itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                 Constants.WGS84_EARTH_FLATTENING, itrf)
        self.station_frame = TopocentricFrame(earth, self.station_geo_point, self.location)
        self.pv_sun = CelestialBodyFactory.getSun()
        self.pv_sun = PVCoordinatesProvider.cast_(self.pv_sun)

    def get_lat(self):
        """
        Getter for the latitude of the sundial.

        Returns
        -------
        flaot
            Latitude of the sundial.

        """
        return self.latitude

    def get_long(self):
        """
        Getter for the longitude of the sundial.

        Returns
        -------
        float
            Longitude of the sundial.

        """
        return self.longitude

    def get_alt(self):
        """
        Getter for the altitude of the sundial.

        Returns
        -------
        float
            Altitude of the sundial.

        """
        return self.altitude

    def get_station_geo_point(self):
        """
        Getter for the Geodetic Point of the sundial station on Earth.

        Returns
        -------
        Geodetic Point
             Geodetic Point of the sundial station on Earth.

        """
        return self.station_geo_point

    def get_station_frame(self):
        """
        Getter of the sundial station frame.

        Returns
        -------
        Topocentric Frame
            Sundial station frame.

        """
        return self.station_frame

    def get_sun_pos_station_frame(self, utc_date):
        """
        Returns the position of the Sun in the sundial station frame.

        Parameters
        ----------
        utc_date : Absolute Date
            Date (in UTC) when the position of the sun should be computed.

        Returns
        -------
        Vector3D
            Position vector of the Sun in the sundial station frame.

        """
        j2000 = FramesFactory.getEME2000()
        sun_pos_j2000 = self.pv_sun.getPVCoordinates(utc_date, j2000).getPosition()
        j2000_to_topocentric = j2000.getTransformTo(self.station_frame, utc_date)
        sun_pos_topo = j2000_to_topocentric.transformPosition(sun_pos_j2000)
        return sun_pos_topo

    def get_sun_elevation(self, utc_date):
        """
        Returns the elevation (in radians) of the Sun at a given date in the sundial reference frame.

        Parameters
        ----------
        utc_date : Absolute Date
            Date in UTC.

        Returns
        -------
        float
            Elevation (in radians) of the Sun.

        """
        return self.station_frame.getElevation(self.get_sun_pos_station_frame(utc_date),
                                               self.station_frame, utc_date)

    def get_sun_elevation_degrees(self, utc_date):
        """
        Returns the elevation (in degrees) of the Sun at a given date in the sundial reference frame.

        Parameters
        ----------
        utc_date : Absolute Date
            Date in UTC.

        Returns
        -------
        float
            Elevation (in radians) of the Sun.

        """
        return degrees(self.station_frame.getElevation(self.get_sun_pos_station_frame(utc_date),
                                                       self.station_frame, utc_date))

    def get_gnomon_shadow_length(self, utc_date):
        """
        Returns the length of the gnomon shadow.

        Parameters
        ----------
        utc_date : Absolute Date
            Date (in UTC) when the length should be computed.

        Returns
        -------
        float
            Length of the gnomon shadow.

        """
        return self.gnomon_length/tan(self.get_sun_elevation(utc_date))

    def get_gnomon_shadow_azimtuh(self, utc_date):
        """
        Returns the azimuth angle (in radians) of the gnomon shadow.

        Parameters
        ----------
        utc_date : Absolute Date
            Date (in UTC) when the azimuth angle should be computed.

        Returns
        -------
        float
            Azimuth (in radians) of the gnonom shadow.

        """
        return self.station_frame.getAzimuth(self.get_sun_pos_station_frame(utc_date),
                                             self.station_frame, utc_date)

    def get_gnomon_shadow_azimtuh_degrees(self, utc_date):
        """
        Returns the azimuth angle (in radians) of the gnomon shadow.

        Parameters
        ----------
        utc_date : Absolute Date
            Date (in UTC) when the azimuth angle should be computed.

        Returns
        -------
        float
            Azimuth (in degrees) of the gnonom shadow.

        """
        return degrees(self.station_frame.getAzimuth(self.get_sun_pos_station_frame(utc_date),
                                                     self.station_frame, utc_date))

    def get_gnomon_shadow_angle_wrtx(self, utc_date):
        """
        Returns the angle (in radians) between the gnomon shadow and the X axis (the East).

        Parameters
        ----------
        utc_date : Absolute Date
            Date when the angle should be computed.

        Returns
        -------
        float
            Angle (in radians) between the gnomon shadow and the X axis (the East)..

        """
        return 3*pi/2 - self.get_gnomon_shadow_azimtuh(utc_date)

    def get_gnomon_shadow_angle_wrtx_degrees(self, utc_date):
        """
        Returns the angle (in degrees) between the gnomon shadow and the X axis (the East).

        Parameters
        ----------
        utc_date : Absolute Date
            Date when the angle should be computed.

        Returns
        -------
        float
            Angle (in degrees) between the gnomon shadow and the X axis (the East).

        """
        return degrees(3*pi/2 - self.get_gnomon_shadow_azimtuh(utc_date))

    def get_sun_max_elevation(self, year, month, day):
        """
        Returns the maximum elevation angle (in radians) of the Sun (that is the zenith) a given day.

        Parameters
        ----------
        year : int
            Year.
        month : int
            Month.
        day : int
            Day.

        Returns
        -------
        current_ele0 : TYPE
            DESCRIPTION.
        currentutc_date : TYPE
            DESCRIPTION.

        """
        currentutc_date = AbsoluteDate(year, month, day, 0, 0, 0.0, TimeScalesFactory.getUTC())
        is_max_elevation = False
        time_step = 10.0
        current_ele0 = self.get_sun_elevation(currentutc_date)
        current_ele1 = self.get_sun_elevation(currentutc_date.shiftedBy(time_step))
        if current_ele0 > current_ele1:
            is_max_elevation = True
        while not is_max_elevation:
            current_ele0 = current_ele1
            current_ele1 = self.get_sun_elevation(currentutc_date.shiftedBy(time_step))
            if current_ele0 > current_ele1:
                is_max_elevation = True
                currentutc_date.shiftedBy(-time_step)
            currentutc_date = currentutc_date.shiftedBy(time_step)
        return current_ele0, currentutc_date


    def get_zenith_local_time(self, year, month, day, summer_time):
        _, zenith_utc = self.get_sun_max_elevation(year, month, day)
        return utc2local_time(zenith_utc, self.time_zone, summer_time)

    def get_zenith_utc(self, year, month, day):
        _, zenith_utc = self.get_sun_max_elevation(year, month, day)
        return zenith_utc

    def get_gnomon_shadow_angle_zenith(self, year, month, day):
        """
        Returns the angle (in radians) between the gnomon shadow and the X axis (the East) at the zenith.

        Parameters
        ----------
        year : int
            Year.
        month : int
            Month.
        day : TYPE
            Day.

        Returns
        -------
        TYPE
            Angle (in radians) between the gnomon shadow and the X axis (the East) at the zenith.

        """
        utc_date = self.get_zenith_utc(year, month, day)
        return self.get_gnomon_shadow_angle_wrtx(utc_date)

    def get_all_gnomon_shadow_angle_degrees(self, year, month, day):
        """
        Returns all (each hour) angles (in degrees) between the gnomon shadow and the X axis (the East)
        beginning with the zenith angle.

        Parameters
        ----------
        year : int
            Year.
        month : int
            Month.
        day : int
            Day.

        Returns
        -------
        gnomon_shadow_angles : list of float
            All (each hour) angles (in degrees) between the gnomon shadow and the X axis (the East) beginning
            beginning with the zenith angle.

        """
        zenith_utc = self.get_zenith_utc(year, month, day)
        gnomon_shadow_angles = [degrees(self.get_gnomon_shadow_angle_zenith(year, month, day))]
        assert self.get_gnomon_shadow_angle_zenith(year, month, day) == \
            self.get_gnomon_shadow_angle_wrtx(zenith_utc)
        for i in range(1, 12):
            gnomon_shadow_angles.append(self.get_gnomon_shadow_angle_wrtx_degrees(zenith_utc.shiftedBy(i*3600.0)))
        return gnomon_shadow_angles

    def get_all_gnomon_shadow_angle(self, year, month, day):
        """
        Returns all (each hour) angles (in radians) between the gnomon shadow and the X axis (the East)
        beginning with the zenith angle.

        Parameters
        ----------
        year : int
            Year.
        month : int
            Month.
        day : int
            Day.

        Returns
        -------
        gnomon_shadow_angles : list of float
            All (each hour) angles (in radians) between the gnomon shadow and the X axis (the East) beginning
            beginning with the zenith angle.

        """
        zenith_utc = self.get_zenith_utc(year, month, day)
        gnomon_shadow_angles = [self.get_gnomon_shadow_angle_zenith(year, month, day)]
        assert self.get_gnomon_shadow_angle_zenith(year, month, day) == \
            self.get_gnomon_shadow_angle_wrtx(zenith_utc)
        for i in range(1, 12):
            gnomon_shadow_angles.append(self.get_gnomon_shadow_angle_wrtx(zenith_utc.shiftedBy(i*3600.0)))
        return gnomon_shadow_angles
