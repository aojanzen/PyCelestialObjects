#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dobby

Version 0.1, 2023-10-21

Dr. Andreas Janzen, janzen@gmx.net, 2023-10-23

Converts positions of celestial objects from right ascension (RA, measured in
hours and minutes) and declination (DEC, measured in degrees and minutes) to
altitude (Alt, measured in degrees) and azimuth (Az, measured in degrees).

Positions of celestial objects are taken from "NGC 2000.0, The Complete New
General Catalogue and Index Catalogue of Nebulae and Star Clusters" by
J.L.E. Dreyer (from https://heasarc.gsfc.nasa.gov/W3Browse/all/ngc2000.html)
"""


import itertools
import sys

import Dobby_bright_stars as stars
import Dobby_deep_sky_objects as dso
import Dobby_celestial_mechanics as cmech


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
# User Interaction
###############################################################################

def print_welcome_message():
    """ Prints welcome screen
    """
    msg = """
    ###########################################################################
    ###                                                                     ###
    ###              Dobby V0.2                                             ###
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


###############################################################################
# Main Program
###############################################################################


def main():
    """ Dobby main function
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

