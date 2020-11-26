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
    if summer_time:
        time_zone += 1
    return utc_date.shiftedBy(3600.0*time_zone)

def local_time2utc(local_time, time_zone, summer_time):
    if summer_time:
        time_zone += 1
    return local_time.shiftedBy(-3600.0*time_zone)

def get_elevation(sun_pos_sundial_frame):
    proj_sun_pos_sundial_frame = Vector3D(sun_pos_sundial_frame.getX(),
                                          sun_pos_sundial_frame.getY(), 0.0)
    cos_elevation = Vector3D.dotProduct(sun_pos_sundial_frame, proj_sun_pos_sundial_frame)/ \
    (sun_pos_sundial_frame.getNorm()*proj_sun_pos_sundial_frame.getNorm())
    return acos(cos_elevation)

def get_azimuth(sun_pos_sundial_frame):
    proj_sun_pos_sundial_frame = Vector3D(sun_pos_sundial_frame.getX(),
                                          sun_pos_sundial_frame.getY(), 0.0)
    cos_azimuth = proj_sun_pos_sundial_frame.getY()/proj_sun_pos_sundial_frame.getNorm()
    return acos(cos_azimuth)


class Sundial:

    def __init__(self, lat=48.58, longi=7.75, alt=142.0, loc="Strasbourg Frame",
                 time_zone=1.0, gnomon_length=1.0):
        self.latitude = radians(lat)
        self.longitude = radians(longi)
        self.altitude = alt
        # utc = TimeScalesFactory.getUTC()
        # date = AbsoluteDate(year, month, 1, 0, 0, 0.0, utc)
        # self.date = date
        self.location= loc
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
        return self.latitude

    def get_long(self):
        return self.longitude

    def get_alt(self):
        return self.altitude

    def get_station_geo_point(self):
        return self.station_geo_point

    def get_station_frame(self):
        return self.station_frame

    def get_sun_pos_station_frame(self, utc_date):
        j2000 = FramesFactory.getEME2000()
        sun_pos_j2000 = self.pv_sun.getPVCoordinates(utc_date, j2000).getPosition()
        j2000_to_topocentric = j2000.getTransformTo(self.station_frame, utc_date)
        sun_pos_topo = j2000_to_topocentric.transformPosition(sun_pos_j2000)
        return sun_pos_topo

    def get_sun_elevation(self, utc_date):
        #Radians
        return self.station_frame.getElevation(self.get_sun_pos_station_frame(utc_date),
                                               self.station_frame, utc_date)

    def get_sun_elevation_degrees(self, utc_date):
        #Degrees
        return degrees(self.station_frame.getElevation(self.get_sun_pos_station_frame(utc_date),
                                                       self.station_frame, utc_date))

    def get_gnomon_shadow_length(self, utc_date):
        return self.gnomon_length/tan(self.get_sun_elevation(utc_date))

    def get_gnomon_shadow_azimtuh(self, utc_date):
        return self.station_frame.getAzimuth(self.get_sun_pos_station_frame(utc_date),
                                             self.station_frame, utc_date)

    def get_gnomon_shadow_azimtuh_degrees(self, utc_date):
        return degrees(self.station_frame.getAzimuth(self.get_sun_pos_station_frame(utc_date),
                                                     self.station_frame, utc_date))

    def get_gnomon_shadow_angle_wrtx(self, utc_date):
        return 3*pi/2 - self.get_gnomon_shadow_azimtuh(utc_date)

    def get_gnomon_shadow_angle_wrtx_degrees(self, utc_date):
        return degrees(3*pi/2 - self.get_gnomon_shadow_azimtuh(utc_date))

    def get_sun_max_elevation(self, year, month, day):
        currentutc_date = AbsoluteDate(year, month, day, 0, 0 , 0.0, TimeScalesFactory.getUTC())
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
        utc_date = self.get_zenith_utc(year, month, day)
        return self.get_gnomon_shadow_angle_wrtx(utc_date)

    def get_all_gnomon_shadow_angle_degrees(self, year, month, day):
        zenith_utc = self.get_zenith_utc(year, month, day)
        gnomon_shadow_angles = [degrees(self.get_gnomon_shadow_angle_zenith(year, month, day))]
        assert self.get_gnomon_shadow_angle_zenith(year, month, day) == \
            self.get_gnomon_shadow_angle_wrtx(zenith_utc)
        for i in range(1, 12):
            gnomon_shadow_angles.append( \
                                        self.get_gnomon_shadow_angle_wrtx_degrees(zenith_utc.shiftedBy(i*3600.0)))
        return gnomon_shadow_angles

    def get_all_gnomon_shadow_angle(self, year, month, day):
        zenith_utc = self.get_zenith_utc(year, month, day)
        gnomon_shadow_angles = [self.get_gnomon_shadow_angle_zenith(year, month, day)]
        assert self.get_gnomon_shadow_angle_zenith(year, month, day) == \
            self.get_gnomon_shadow_angle_wrtx(zenith_utc)
        for i in range(1, 12):
            gnomon_shadow_angles.append(self.get_gnomon_shadow_angle_wrtx(zenith_utc.shiftedBy(i*3600.0)))
        return gnomon_shadow_angles
