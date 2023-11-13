#!/usr/bin/env python3

"""
Dobby_bright_stars.py

Part of Dobby, a tool to turn a Dobsonian telescope into a go-to platform

Provides data structures and methods to load and store a list of bright stars
used as alignment objects in the night sky

Author: Dr. Andreas Janzen, janzen@gmx.net
Date: 2023-11-12
"""


import collections
import csv

import Dobby_celestial_mechanics as cmech


Star = collections.namedtuple("Star", "Rank, Name, Beyer, Mag_v, RA, DEC, Alt, Az")


# Input file contains 48 brightest stars (according to Hipparcos project) with
# positions manually taken from Dreyer's "NGC 2000.0" catalogue
INPUT_FILE = "48_Brightest_Stars_Hipparcos.csv"

# Altitude angle below which stars are considered unsuitable for alignment,
# e.g. due to limited visibility (houses, trees, etc.) or atmospheric seeing
LOW_ANGLE = 20.0

# Global variable that contains the list of bright stars
star_list = list()


###############################################################################
# Data handling
###############################################################################

def load_bright_stars(fname = INPUT_FILE):
    """ Loads input file with bright star positions, converts (RA, DEC) positions
    to (Alt, Az), and adds them to the global variable star_list as Star namedtuples
    """
    global star_list # allows this function to change the global variable

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

