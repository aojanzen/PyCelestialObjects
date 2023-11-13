#!/usr/bin/env python3

"""
Dobby_celestial_mechanics.py

Part of Dobby, a tool to turn a Dobsonian telescope into a go-to platform

Provides data structures and methods for calculations of celestial mechanics,
e.g. conversion from right ascension and declination to altitude and azimuth

Author: Dr. Andreas Janzen, janzen@gmx.net
Date: 2023-11-12
"""


import datetime
import math


TimePlace = collections.namedtuple("TimePlace",
        "Year, Month, Day, Hours, Minutes, Seconds, Latitude, Longitude")


# Use Bad Camberg, Germany, as preset for location site
LAT = 50.314 # 50.314°N
LONG = 8.255 # 8.255°E

# Use current date and system time as presets
current = datetime.datetime.now()
YEAR, MONTH, DAY = current.year, current.month, current.day
HOURS, MINUTES, SECONDS = current.hour, current.minute, current.second

TIME_ZONE = 1 # local time is UT + 1h
DST = False # daylight saving time yes/no?

if DST:
    UT_COR = -(TIME_ZONE + 1)
else:
    UT_COR = -TIME_ZONE

DEFAULT = TimePlace(
    YEAR,
    MONTH,
    DAY,
    HOURS,
    MINUTES,
    SECONDS,
    LAT,
    LONG,
)


###############################################################################
# Helper functions
###############################################################################

def HMS_to_decimal_time(hours, minutes, seconds):
    """ Converts hours, minutes, seconds to a time measured in fractional hours
    """
    return hours + minutes/60 + seconds/3600


def rad2deg(angle):
    """ Converts an angle from radians to degrees
    """
    return 180 * angle / math.pi


def deg2rad(angle):
    """ Converts an angle from degrees to radians
    """
    return math.pi * angle / 180


def deg2time(deg):
    """ Converts a rotation angle to hours, minutes, seconds
    """
    rel_day = deg / 360.0
    hours, frac_hours = divmod(rel_day * 24.0, 1)
    minutes, frac_minutes = divmod(frac_hours * 60, 1)
    seconds = frac_minutes * 60
    return int(hours), int(minutes), int(seconds)


def RAStr2RA(RAStr):
    """ Converts an angle given in the format hours:minutes:seconds to a
    decimal angle, using the conversion 1 hr --> 15 deg
    """
    tmp = RAStr.split(":")
    h, m, s = float(tmp[0][:-1]), float(tmp[1][:-1]), float(tmp[2][:-1])
    return (h + m/60.0 + s/3600.0)*15.0


def DECStr2DEC(DECStr):
    """ Converts an angle given as a string of degrees, arcminutes and arcseconds
    to a decimal angle. The sign of the angle is taken into account.
    """
    tmp = DECStr.split(":")
    deg = tmp[0][:-3]
    arcmin = tmp[1][:-2]
    arcsec = tmp[2][:-2]
    if deg[0] == "-":
        multiplier = -1.0
        deg = deg[1:]
    else:
        multiplier = 1.0
    return multiplier * (float(deg) + float(arcmin)/60.0 + float(arcsec)/3600.0)


def cartesian(star):
    """ Returns the position of a star in cartesian coordinates
    """
    Alt = math.radians(star.Alt)
    Az = math.radians(star.Az)

    x = math.cos(Alt) * math.cos(Az)
    y = math.cos(Alt) * math.sin(Az)
    z = math.sin(Alt)

    return x, y, z


def vector_product(v1, v2):
    """ Returns the vector product of two vectors in cartesian coordinates
    """
    r0 = v1[1] * v2[2] - v1[2] * v2[1]
    r1 = v1[2] * v2[0] - v1[0] * v2[2]
    r2 = v1[0] * v2[1] - v1[1] * v2[0]

    return r0, r1, r2


def vector_norm(v):
    """ Returns the cartesian norm of a vector
    """
    sum_sq = 0
    for component in v:
        sum_sq += component**2
    return math.sqrt(sum_sq)


def area(s1, s2, s3):
    """ Returns the area of the planar triangle created by the position vectors
    of three stars s1, s2 and s3 on a sphere with unity radius
    """
    v1 = cartesian(s1)
    v2 = cartesian(s2)
    v3 = cartesian(s3)

    edge1 = (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])
    edge2 = (v1[0]-v3[0], v1[1]-v3[1], v1[2]-v2[2])

    area = 0.5 * vector_norm( vector_product(edge1, edge2) )

    return area


def angle(s1, s2):
    """ Returns the angle between two stars s1, s2
    """
    alt1 = math.radians(s1.Alt)
    az1 = math.radians(s1.Az)

    alt2 = math.radians(s2.Alt)
    az2 = math.radians(s2.Az)

    cos_angle = math.cos(alt1)*math.cos(alt2)*math.cos(az1-az2) + \
            math.sin(alt1)*math.sin(alt2)

    return math.acos(cos_angle)


###############################################################################
# Time calculations
###############################################################################

def days_from_J2000_until_year(year):
    """ Returns the number of days from J2000 (January, 1st, 2000, 12 AM
    on 0th meridian (i.e. 12AM UT)) until the beginning of the present year.
    Takes into account whether the current year is a leap year or not.
    """
    return (year-1998)*365 + (year-1998+1)//4 - 731.5


def days_until_month_begins(year, month):
    """ Prints the number of days from the beginning of the year until
    the beginning of the current month, e.g. 31 for February. Takes
    into account whether the current year is a leap year or not.
    """
    DAYS_TIL_MONTH = {
        1   : 0,
        2   : 31,
        3   : 59,
        4   : 90,
        5   : 120,
        6   : 151,
        7   : 181,
        8   : 212,
        9   : 243,
        10  : 273,
        11  : 304,
        12  : 334
    }
    if year % 4 == 0 and month > 2:
        return DAYS_TIL_MONTH[month] + 1
    else:
        return DAYS_TIL_MONTH[month]


def days_from_J2000(year, month, day, hours, minutes, seconds):
    """ Returns the number of days as a decimal fraction from J2000
    (January, 1, 2000, 12 AM UT) until the specified time.
    days1, days2, days3:
    - #days from J2000 until the beginning of the specified year
    - #days from the beginning of the year to the beginning of the specified day
    - fraction that has passed of the specified day until the specified time
    """
    days1 = days_from_J2000_until_year(year)
    days2 = days_until_month_begins(year, month) + day
    days3 = HMS_to_decimal_time(hours, minutes, seconds) / 24.0
    return days1 + days2 + days3


def local_siderial_time(days_from_J2000, longitude, time):
    """ Returns the local siderial time in degrees at the specified
    longitude and universal time for a number of days since J2000
    given as a decimal fraction including the fraction of the specified
    day. The approximation is within 0.3 seconds of time for dates within
    100 years of J2000.
    """
    LST = 100.46 + 0.985647*days_from_J2000 + longitude + 15*time
    return LST % 360.0


def hour_angle(RA, LST):
    """ Returns the hour angle of an object with right ascension RA at
    local siderial time LST. HA and RA are measured in degrees.
    """
    return (LST - RA) % 360


def RA_DEC_to_ALT_AZ(raD, decD, haD, latD):
    """ Converts (RA, DEC) coordinates to (ALT, AZ) for a given hour angle
    and latitude. All angles are measured in degrees.
    """
    dec = deg2rad(decD)
    ha = deg2rad(haD)
    lat = deg2rad(latD)

    sin_alt = (math.sin(dec) * math.sin(lat)) + (math.cos(dec) * math.cos(lat) * math.cos(ha))
    alt = math.asin(sin_alt)

    cos_a = (math.sin(dec) - math.sin(alt) * math.sin(lat)) / (math.cos(alt) * math.cos(lat))
    a = math.acos(cos_a)

    alt = rad2deg(alt)
    if math.sin(ha) < 0:
        az = rad2deg(a)
    else:
        az = 360 - rad2deg(a)

    return alt, az


def set_time_and_place():
    """ Set time, date, latitude and longitude of observation site
    """
    print("\nDefault time and date:")
    print(f"{DEFAULT.Year}-{DEFAULT.Month}-{DEFAULT.Day},",
          f"{DEFAULT.Hours}:{DEFAULT.Minutes}:{DEFAULT.Seconds},", end = " ")
    if DEFAULT.Latitude < 0.0:
        print(f"{-DEFAULT.Latitude}°S / ", end = "")
    else:
        print(f"{DEFAULT.Latitude}°N / ", end = "")
    if DEFAULT.Longitude < 0.0:
        print(f"{-DEFAULT.Longitude}°W")
    else:
        print(f"{DEFAULT.Longitude}°E")

    change = input("\nDo you want to change the default time and place? (y/n) > ")
    if change.lower() in ["n", "no"]:
        return DEFAULT

    while True:
        try:
            latitude = float(input("Geographic latitude (North positive)? (°) > "))
            longitude = float(input("Geographic longitude (East positive)? (°) > "))
            if not (-90.0 <= latitude <= 90.0) and (-180.0 <= longitude <= 180.0):
                print("\n### Error in geographical position")
                continue

            year, month, day = input("Current date (yyyy-mm-dd) > ").split("-")
            yi = int(year)
            mi = int(month)
            di = int(day)
            if not ((2023 <= yi <= 2073) and (1 <= mi <= 12) and (1 <= di <= 31)):
                print("\n### Error in date\n")
                continue

            time_zone = int(input("What is your local time zone (UT+...)? > "))
            if not (0 <= time_zone <= 23):
                print("\n### Error in time zone\n")
                continue

            dst = input("Is it daylight saving time? (y/n) > ")
            if dst.lower() in ["y", "yes"]:
                ut_cor = -(time_zone + 1)
            elif dst.lower() in ["n", "no"]:
                ut_cor = -time_zone
            else:
                print("\n### Error in daylight saving time\n")
                continue

            hs, ms, ss = input("Local time (hh:mm:ss) > ").split(":")
            h = int(hs)
            m = int(ms)
            s = int(ss)
            if not ((0 <= h <= 23) and (0 <= m <= 59) and (0 <= s <= 59)):
                print("\n### Error in time\n")
                continue

            break

        except:
            print("### Incorrect data. Please try again!")

    return TimePlace(yi, mi, di, h + ut_cor, m, s, latitude, longitude)

