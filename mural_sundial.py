from math import radians, degrees, tan, pi
from sundial import Sundial, get_elevation, get_azimuth
import orekit

vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()

from org.orekit.frames import Frame, FramesFactory, Transform
# from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.utils import AngularCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D, Rotation, RotationConvention
from org.orekit.time import AbsoluteDate


class MuralSundial(Sundial):
    """
    """

    def __init__(self, lat=48.58, longi=7.75, alt=142.0, loc="Strasbourg", time_zone=1,
                 gnomon_length=1.0, mural_angle=45.0, height=1.0, orient="Nord"):
        """
        

        Parameters
        ----------
        lat : float, optional
            DESCRIPTION. The default is 48.58.
        longi : float, optional
            DESCRIPTION. The default is 7.75.
        alt : float, optional
            DESCRIPTION. The default is 142.0.
        loc : str, optional
            DESCRIPTION. The default is "Strasbourg".
        time_zone : int, optional
            DESCRIPTION. The default is 1.
        gnomon_length : float, optional
            DESCRIPTION. The default is 1.0.
        mural_angle : float, optional
            DESCRIPTION. The default is 45.0.
        height : float, optional
            DESCRIPTION. The default is 1.0.
        orient : str, optional
            DESCRIPTION. The default is "Nord".

        Returns
        -------
        None.

        """
        super().__init__(lat, longi, alt, loc, time_zone, gnomon_length)
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
        self.sundial_frame = Frame(self.station_frame, sundial_transform, "Sundial Frame")

    def get_mural_angle(self):
        """
        getter for the mural angle

        Returns
        -------
        float
            Mural angle.

        """
        return self.mural_angle

    def get_height(self):
        return self.height

    def get_sundial_frame(self):
        return self.sundial_frame

    def get_sun_pos_sundial_frame(self, utc_date):
        j2000 = FramesFactory.getEME2000()
        sun_pos_j2000 = self.pv_sun.getPVCoordinates(utc_date, j2000).getPosition()
        j2000_to_sundial = j2000.getTransformTo(self.sundial_frame, utc_date)
        sun_pos_sundial_frame = j2000_to_sundial.transformPosition(sun_pos_j2000)
        return sun_pos_sundial_frame

    def get_sun_elevation_sundial_frame(self, utc_date):
        #Radians
        return get_elevation(self.get_sun_pos_sundial_frame(utc_date))

    def get_sun_elevation_sundial_frame_degrees(self, utc_date):
        #Degrees
        return degrees(get_elevation(self.get_sun_pos_sundial_frame(utc_date)))

    def get_gnomon_shadow_length(self, utc_date):
        return self.gnomon_length/tan(self.get_sun_elevation_sundial_frame_degrees(utc_date))

    def get_gnomon_shadow_azimuth(self, utc_date):
        return get_azimuth(self.get_sun_pos_sundial_frame(utc_date))

    def get_gnomon_shadow_azimuth_degrees(self, utc_date):
        return degrees(get_azimuth(self.get_sun_pos_sundial_frame(utc_date)))

    def get_gnomon_shadow_angle_wrtx(self, utc_date):
        return 3*pi/2 - self.get_gnomon_shadow_azimuth(utc_date)

    def get_gnomon_shadow_angle_wrtx_degrees(self, utc_date):
        return degrees(3*pi/2 - self.get_gnomon_shadow_azimuth(utc_date))

    def get_gnomon_shadow_angle_zenith(self, year, month, day):
        utc_date = self.get_zenith_utc(year, month, day)
        return self.get_gnomon_shadow_angle_wrtx(utc_date)

    def get_gnomon_shadow_angle_zenith_degrees(self, year, month, day):
        utc_date = self.get_zenith_utc(year, month, day)
        return self.get_gnomon_shadow_angle_wrtx_degrees(utc_date)

    def get_all_gnomon_shadow_angle_degrees(self, year, month, day):
        zenith_utc = self.get_zenith_utc(year, month, day)
        gnomon_shadow_angles = [degrees(self.get_gnomon_shadow_angle_zenith(year, month, day))]
        for i in range(1, 12):
            gnomon_shadow_angles.append(self.get_gnomon_shadow_angle_wrtx_degrees(zenith_utc.shiftedBy(i*3600.0)))
        return gnomon_shadow_angles
