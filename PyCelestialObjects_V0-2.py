#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyCelestialObjects

Version 0.1, 2023-10-21

Dr. Andreas Janzen, janzen@gmx.net, 2023-10-23

Converts positions of celestial objects from right ascension (RA, measured in
hours and minutes) and declination (DEC, measured in degrees and minutes) to
altitude (Alt, measured in degrees) and azimuth (Az, measured in degrees).

Positions of celestial objects are taken from "NGC 2000.0, The Complete New
General Catalogue and Index Catalogue of Nebulae and Star Clusters" by
J.L.E. Dreyer (from https://heasarc.gsfc.nasa.gov/W3Browse/all/ngc2000.html)
"""


import collections
import csv
import datetime
import itertools
import math
import sys


Star = collections.namedtuple("Star", "Rank, Name, Beyer, Mag_v, RA, DEC, Alt, Az")
TimePlace = collections.namedtuple("TimePlace",
        "Year, Month, Day, Hours, Minutes, Seconds, Latitude, Longitude")


# Input file contains 48 brightest stars (according to Hipparcos project) with
# positions manually taken from Dreyer's "NGC 2000.0" catalogue
INPUT_FILE = "48_Brightest_Stars_Hipparcos.csv"
DSO_FILE = "ngc2000_final.csv"

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


# Altitude angle below which stars are considered unsuitable for alignment,
# e.g. due to limited visibility (houses, trees, etc.) or atmospheric seeing
LOW_ANGLE = 20.0


# Menu entries
MENUS = {
        "MAIN": """Align mount
                   Navigate to object
                   Show object information
                   Set time and place
                   Exit program""",

        "STAR ALIGNMENT": """Print bright stars
                             1-star alignment
                             2-star alignment
                             3-star alignment
                             Back""",
}

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


###############################################################################
# Data handling
###############################################################################

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


def load_bright_stars(fname = INPUT_FILE):
    """ Loads input file with bright star positions, converts (RA, DEC) positions
    to (Alt, Az), and returns list with all stars stored as Star namedtuples
    """
    all_stars = list()

    days_J2000 = days_from_J2000(YEAR, MONTH, DAY, HOURS + UT_COR, MINUTES, SECONDS)
    lst = local_siderial_time(days_J2000, LONG,
                              HMS_to_decimal_time(HOURS + UT_COR, MINUTES, SECONDS))

    with open(fname) as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for line in reader:
            if line[0] == "Rank":
                continue
            else:
                rank = int(line[0])
                name = line[1]
                beyer = line[2]
                mag = float(line[3])
                ra = RAStr2RA(line[4])
                dec = DECStr2DEC(line[5])
                ha = hour_angle(ra, lst)
                Alt, Az = RA_DEC_to_ALT_AZ(ra, dec, ha, LAT)
                current = Star(rank, name, beyer, mag, ra, dec, Alt, Az)
                all_stars.append(current)

    return all_stars


def best_stars(stars):
    """ Prints a list of star combinations sorted by the size of the planar
    triangle that the stars span on the sky.
    """
    NUMBER = 15
    triples = dict()

    for i,j,k in itertools.combinations(range(len(stars)), 3):
        triples[(i,j,k)] = area(stars[i], stars[j], stars[k])

    triples_sorted = sorted(triples.items(), key = lambda tr:tr[1], reverse = True)

    heading = f"{NUMBER} best suited triples out of {len(triples_sorted)} combinations:"
    print("\n" + heading)
    print()

    index = 1
    for triple in triples_sorted:
        txt = f"{stars[triple[0][0]].Name:<15} {stars[triple[0][1]].Name:<15}" + \
              f"{stars[triple[0][2]].Name:<15} area: {triple[1]:>5.3f}"
        print(f" {index:>2d}. ", end="")
        print(txt)

        index += 1
        if index > NUMBER:
            break

    return None


###############################################################################
# User Interaction
###############################################################################

def print_welcome_message():
    """ Prints welcome screen
    """
    msg = """
    ###########################################################################
    ###                                                                     ###
    ###              PyCelestialObjects V0.1                                ###
    ###                                                                     ###
    ###              Dr. Andreas Janzen, November 2023                      ###
    ###                                                                     ###
    ###########################################################################
    """
    print(msg)

    return None


def print_menu(menu):
    """ Prints a menu from global variable MENUS
    """
    heading = f"{menu} MENU"
    print("\n\n" + heading)
    print("="* len(heading) + "\n")

    entries = [entry.strip() for entry in MENUS[menu].split("\n")]
    for k, v in enumerate(entries, 1):
        print(f"  ({k}) {v}")
    print()

    while True:
        try:
            choice = int(input("> "))
            if 1 <= choice <= len(entries):
                break
        except:
            pass

    return choice


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


def get_suitable_stars():
    """ Returns a list of currently visible bright stars. Stars that are low in
    the sky, i.e. below 20° altitude, are not included, even if they are
    visible from the observation site.
    """
    all_stars = load_bright_stars(INPUT_FILE)
    suitable = list()

    for star in all_stars:
        if star.Alt >= LOW_ANGLE:
            suitable.append(star)

    sorted_stars = sorted(suitable, key = lambda s:s.Alt, reverse = True)

    return sorted_stars # suitable


def print_star_list(stars, description):
    """ Print a list of stars with their full information set
    """
    heading = f"\n\nList of {description} ({len(stars)}):"
    print(heading)
    print("=" * (len(heading)-2))

    print(f"\nDate: {YEAR:4d}-{MONTH:02d}-{DAY:02d}, local time:",
          f"{HOURS:02d}:{MINUTES:02d}:{SECONDS:02d} (UT+{-UT_COR}h),",
          f"position: {LAT:>5.2f}°", end="")
    if LAT >= 0.0:
        print("N", end="")
    else:
        print("S", end="")
    print(f" / {LONG:>5.2f}°", end="")
    if LONG >= 0.0:
        print("E\n")
    else:
        print("W\n")

    for star in stars:
        if star.Alt <= 0.0:
            continue
        elif star.Alt < LOW_ANGLE:
            print(" (", end="")
            low += 1
        else:
            print("  ", end="")
        print(f"{star.Name:<15} {star.Beyer:>11}     {star.Mag_v:5.2f} mag    ",
              f"Az = {star.Az:>6.2f}° / Alt = {star.Alt:>5.2f}°", end="")
        if star.Alt < 20.0:
            print(")")
        else:
            print()

    return None


###############################################################################
# Main Program
###############################################################################


def main():
    """ PyCelestialObjects main function
    """
    print_welcome_message()

    while True:
        choice = print_menu("MAIN")
        if choice == 1: # Align mount
            suitable = get_suitable_stars()
            print_star_list(suitable, "suitable stars")
            best_stars(suitable)
        elif choice == 2: # Navigate to object -- requires hardware programming
            print("\n### Function is not implemented yet.")
        elif choice == 3: # Show object information
            print("\n### Function is not implemented yet.")
        elif choice == 4: # Basic settings
            set_time_and_place()
        elif choice == 5:
            sys.exit(0)
        else:
            print("\n\nHow could I even get here??")
            sys.exit(0)


if __name__ == "__main__":
    main()

