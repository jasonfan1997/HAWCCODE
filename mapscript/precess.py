# Functions to allow the plotters to make a precession.

import astropy.units as apUnits
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.time import Time

from healpy.rotator import Rotator
from healpy.rotator import vec2dir
from healpy.rotator import dir2vec

from numpy import array
from numpy import deg2rad
from numpy import rad2deg
from numpy import arctan2
from numpy import arcsin
from numpy import fabs

def mjd2J2000ang(mjd, coord=None, rot=None):

    return EulerAngles(equinox=Time(mjd,format='mjd'), targetEquinox='J2000', coord=coord, rot=rot)

def EulerAngles(equinox, targetEquinox, coord=None, rot=None):
    #Return the Euler angles (passive Z(-Y)X) to go from equinox to targetEquinox in degrees
    #Rot is an extra rotation, after precessing
    #If coord is provided, the inverse coordinate transformation is applied first, then the precession, then the direct coordinate transformation, and then the rotation is calculated

    deg=apUnits.degree;

    #Unit vector in current epoch
    x0=array([1,0,0])
    y0=array([0,1,0])
    z0=array([0,0,1])

    #Go to original coordinates
    rCoord=Rotator(coord=coord)
    x=rCoord.I(x0)
    y=rCoord.I(y0)
    z=rCoord.I(z0)
    
    #Change to longitude and latitude
    xLonLat=vec2dir(x,lonlat=True)
    yLonLat=vec2dir(y,lonlat=True)
    zLonLat=vec2dir(z,lonlat=True)

    #Put them in astropy
    xAP=SkyCoord(ra=xLonLat[0]*deg,dec=xLonLat[1]*deg,frame=FK5(equinox=equinox))
    yAP=SkyCoord(ra=yLonLat[0]*deg,dec=yLonLat[1]*deg,frame=FK5(equinox=equinox))
    zAP=SkyCoord(ra=zLonLat[0]*deg,dec=zLonLat[1]*deg,frame=FK5(equinox=equinox))

    #Precess
    xpAP=xAP.transform_to(FK5(equinox=targetEquinox))
    ypAP=yAP.transform_to(FK5(equinox=targetEquinox))
    zpAP=zAP.transform_to(FK5(equinox=targetEquinox))

    #Go back to vectors
    xp=dir2vec(theta=xpAP.ra.degree, phi=xpAP.dec.degree, lonlat=True)
    yp=dir2vec(theta=ypAP.ra.degree, phi=ypAP.dec.degree, lonlat=True)
    zp=dir2vec(theta=zpAP.ra.degree, phi=zpAP.dec.degree, lonlat=True)

    #Go back to target coordinates
    xp=rCoord(xp)
    yp=rCoord(yp)
    zp=rCoord(zp)

    #Perform final rotation from argument rot
    rf = Rotator(rot=rot)
    xp=array(rf(xp))
    yp=array(rf(yp))
    zp=array(rf(zp))
    
    #Get Euler angles
    if fabs(zp[0])==1:
        alpha=0.0
        beta=zp[0]*pi/2
        gamma=arctan2(xp[1],xp[2])
    else:
        alpha=arctan2(yp[0],xp[0])
        beta=arcsin(zp[0])
        gamma=arctan2(zp[1],zp[2])
        
    return (rad2deg(alpha),rad2deg(beta),rad2deg(gamma))

def precess(ra,dec,equinox,targetEquinox):
    #Function to check that the above one actually works

    vec=dir2vec(ra,dec,lonlat=True)

    R=Rotator(rot=EulerAngles(equinox,targetEquinox))

    vec=R(vec)

    [oRA,oDec]=vec2dir(vec,lonlat=True)

    while oRA<0:
        oRA = oRA+360
    
    
    return (oRA,oDec)
