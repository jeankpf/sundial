#### Programme ne marche pas pour le mural, impression que le resultat doit etre tourner de 180 deg, a verifier,



import orekit
import numpy

vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()

from org.orekit.frames import Frame, FramesFactory, TopocentricFrame, Transform
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.utils import IERSConventions, Constants, PVCoordinatesProvider, AngularCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D, Rotation, RotationConvention
from org.orekit.time import AbsoluteDate, TimeScalesFactory

from math import radians, degrees, tan, cos, sin, pi

from sundial import Sundial, UTC2LocalTime, localTime2UTC, get_elevation, get_azimuth

class MuralSundial(Sundial):

    def __init__(self, lat=48.58, longi=7.75, alt=142.0, loc="Strasbourg", timeZone=1, gnomonLength=1.0, mural_angle=45.0, height= 1.0, orient="Nord"):
        super().__init__(lat, longi, alt, loc, timeZone, gnomonLength)
        self.mural_angle = radians(mural_angle)
        self.height = height
        if orient == "Nord":
            orientation = -1
        elif orient == "Sud":
            orientation = 1
        else:
            orientation = 0
            print("Cannot read the orientation, put to 0")
        z_translation = Vector3D(0.0, 0.0, height)
        whenever = AbsoluteDate()
        z_transformation = Transform(whenever, z_translation)
        z_vector = Vector3D(0.0, 0.0, 1.0)
        horizontal_rotation = Rotation(z_vector, mural_angle, RotationConvention.VECTOR_OPERATOR)
        horizontal_ac = AngularCoordinates(horizontal_rotation, Vector3D(0.0, 0.0, 0.0))
        z_rotation_tf = Transform(whenever, horizontal_ac)
        mural_transform = Transform(whenever, z_transformation, z_rotation_tf)
        x_vector = Vector3D(1.0, 0.0, 0.0)
        xm_rotation = Rotation(x_vector, orientation*pi/2, RotationConvention.VECTOR_OPERATOR)
        xm_rot_ac = AngularCoordinates(xm_rotation, Vector3D(0.0, 0.0, 0.0))
        xm_rot_tf = Transform(whenever, xm_rot_ac)
        sundial_transform = Transform(whenever, mural_transform, xm_rot_tf)
        self.sundial_frame = Frame(self.stationFrame, sundial_transform, "Sundial Frame")

    def get_mural_angle(self):
        return self.mural_angle

    def get_height(self):
        return self.height

    def get_sundial_frame(self):
        return self.sundial_frame

    def get_sun_pos_sundial_frame(self, UTCdate):
        j2000 = FramesFactory.getEME2000()
        sunPosInJ2000 = self.pvSun.getPVCoordinates(UTCdate, j2000).getPosition()
        j2000_to_sundial = j2000.getTransformTo(self.sundial_frame, UTCdate)
        sun_pos_sundial_frame = j2000_to_sundial.transformPosition(sunPosInJ2000)
        return sun_pos_sundial_frame

    def get_sun_elevation_sundial_frame(self, UTCdate):
        #Radians
        return get_elevation(self.get_sun_pos_sundial_frame(UTCdate))

    def get_sun_elevation_sundial_frame_degrees(self, UTCdate):
        #Degrees
        return degrees(get_elevation(self.get_sun_pos_sundial_frame(UTCdate)))

    def get_gnomon_shadow_length(self, UTCdate):
        return self.gnomonLength/tan(self.get_sun_elevation_sundial_frame_degrees(UTCdate))

    def get_gnomon_shadow_azimuth(self, UTCdate):
        return get_azimuth(self.get_sun_pos_sundial_frame(UTCdate))

    def get_gnomon_shadow_azimuth_degrees(self, UTCdate):
        return degrees(get_azimuth(self.get_sun_pos_sundial_frame(UTCdate)))

    def get_gnomon_shadow_angleWRTX(self, UTCdate):
        return 3*pi/2 - self.get_gnomon_shadow_azimuth(UTCdate)

    def get_gnomon_shadow_angleWRTX_degrees(self, UTCdate):
        return degrees(3*pi/2 - self.get_gnomon_shadow_azimuth(UTCdate))

    def get_gnomon_shadow_angle_zenith(self, year, month, day):
        UTCdate = self.getZenithUTC(year, month, day)
        return self.get_gnomon_shadow_angleWRTX(UTCdate)

    def get_gnomon_shadow_angle_zenith_degrees(self, year, month, day):
        UTCdate = self.getZenithUTC(year, month, day)
        return self.get_gnomon_shadow_angleWRTX_degrees(UTCdate)

    def get_all_gnomon_shadow_angle_degrees(self, year, month, day):
        zenithUTC = self.getZenithUTC(year, month, day)
        gnomonShadowAngles = [degrees(self.get_gnomon_shadow_angle_zenith(year, month, day))]
        # assert self.getGnomonShadowAngleZenith(year, month, day) == self.getGnomonShadowAngleWRTX(zenithUTC)
        for i in range(1, 12):
            gnomonShadowAngles.append(self.get_gnomon_shadow_angleWRTX_degrees(zenithUTC.shiftedBy(i*3600.0)))
        return gnomonShadowAngles