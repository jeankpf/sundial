# matplotlib inline
# pylab inline

import orekit
import numpy

vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()

from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.utils import IERSConventions, Constants, PVCoordinatesProvider
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.time import AbsoluteDate, TimeScalesFactory

from math import radians, degrees, tan, cos, sin, pi

def UTC2LocalTime(UTCdate, timeZone, summerTime):
    if summerTime:
        timeZone += 1
    return UTCdate.shiftedBy(3600.0*timeZone)

def localTime2UTC(localTime, timeZone, summerTime):
    if summerTime:
        timeZone += 1
    return localTime.shiftedBy(-3600.0*timeZone)

class Sundial:

    def __init__(self, lat=48.58, longi=7.75, alt=142.0, loc="Strasbourg", timeZone=1, gnomonLength=1):
        self.latitude = radians(lat)
        self.longitude = radians(longi)
        self.altitude = alt
        # utc = TimeScalesFactory.getUTC()
        # date = AbsoluteDate(year, month, 1, 0, 0, 0.0, utc)
        # self.date = date
        self.location= loc
        self.timeZone = timeZone
        # self.summerTime = summerTime
        self.gnomonLength = gnomonLength
        self.stationGeoPoint = GeodeticPoint(self.latitude, self.longitude, self.altitude)
        ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, ITRF)
        self.stationFrame = TopocentricFrame(earth, self.stationGeoPoint, self.location)
        self.pvSun = CelestialBodyFactory.getSun()
        self.pvSun = PVCoordinatesProvider.cast_(self.pvSun)

    def getLat(self):
        return self.latitude

    def getLong(self):
        return self.longitude

    def getAlt(self):
        return self.altitude

    def getStationGeoPoint(self):
        return self.stationGeoPoint

    def getStationFrame(self):
        return self.stationFrame

    def getSunPosInStationFrame(self, UTCdate):
        j2000 = FramesFactory.getEME2000()
        sunPosInJ2000 = self.pvSun.getPVCoordinates(UTCdate, j2000).getPosition()
        j2000ToTopocentric = j2000.getTransformTo(self.stationFrame, UTCdate)
        sunPosInTopo = j2000ToTopocentric.transformPosition(sunPosInJ2000)
        return sunPosInTopo

    def getSunElevation(self, UTCdate):
        #Radians
        return self.stationFrame.getElevation(self.getSunPosInStationFrame(UTCdate), self.stationFrame, UTCdate)

    def getSunElevationDegrees(self, UTCdate):
        #Degrees
        return degrees(self.stationFrame.getElevation(self.getSunPosInStationFrame(UTCdate), self.stationFrame, UTCdate))

    def getGnomonShadowLength(self, UTCdate):
        return self.gnomonLength/tan(self.getSunElevation(UTCdate))

    def getGnomonShadowAzimtuh(self, UTCdate):
        return self.stationFrame.getAzimuth(self.getSunPosInStationFrame(UTCdate), self.stationFrame, UTCdate)

    def getGnomonShadowAzimtuhDegrees(self, UTCdate):
        return degrees(self.stationFrame.getAzimuth(self.getSunPosInStationFrame(UTCdate), self.stationFrame, UTCdate))

    def getGnomonShadowAngleWRTX(self, UTCdate):
        return 3*pi/2 - self.getGnomonShadowAzimtuh(UTCdate)

    def getGnomonShadowAngleWRTXDegrees(self, UTCdate):
        return degrees(3*pi/2 - self.getGnomonShadowAzimtuh(UTCdate))

    def getSunMaxElevation(self, year, month, day):
        currentUTCDate = AbsoluteDate(year, month, day, 0, 0 , 0.0, TimeScalesFactory.getUTC())
        isMaxElevation = False
        timeStep = 10.0
        currentEle0 = self.getSunElevation(currentUTCDate)
        currentEle1 = self.getSunElevation(currentUTCDate.shiftedBy(timeStep))
        if currentEle0 > currentEle1:
            isMaxElevation = True
        while not isMaxElevation:
            currentEle0 = currentEle1
            currentEle1 = self.getSunElevation(currentUTCDate.shiftedBy(timeStep))
            if currentEle0 > currentEle1:
                isMaxElevation = True
                currentUTCDate.shiftedBy(-timeStep)
            currentUTCDate = currentUTCDate.shiftedBy(timeStep)
        return currentEle0, currentUTCDate


    def getZenithLocalTime(self, year, month, day, summertime):
        _, zenithUTC = self.getSunMaxElevation(year, month, day)
        return UTC2LocalTime(zenithUTC, self.timeZone, summertime)

    def getZenithUTC(self, year, month, day):
        _, zenithUTC = self.getSunMaxElevation(year, month, day)
        return zenithUTC

    def getGnomonShadowAngleZenith(self, year, month, day):
        UTCdate = self.getZenithUTC(year, month, day)
        return self.getGnomonShadowAngleWRTX(UTCdate)

    def getAllGnomonShadowAngleDegrees(self, year, month, day):
        zenithUTC = self.getZenithUTC(year, month, day)
        gnomonShadowAngles = [degrees(self.getGnomonShadowAngleZenith(year, month, day))]
        assert self.getGnomonShadowAngleZenith(year, month, day) == self.getGnomonShadowAngleWRTX(zenithUTC)
        for i in range(1, 12):
            gnomonShadowAngles.append(self.getGnomonShadowAngleWRTXDegrees(zenithUTC.shiftedBy(i*3600.0)))
        return gnomonShadowAngles

    def getAllGnomonShadowAngle(self, year, month, day):
        zenithUTC = self.getZenithUTC(year, month, day)
        gnomonShadowAngles = [self.getGnomonShadowAngleZenith(year, month, day)]
        assert self.getGnomonShadowAngleZenith(year, month, day) == self.getGnomonShadowAngleWRTX(zenithUTC)
        for i in range(1, 12):
            gnomonShadowAngles.append(self.getGnomonShadowAngleWRTX(zenithUTC.shiftedBy(i*3600.0)))
        return gnomonShadowAngles