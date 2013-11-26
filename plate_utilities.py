#! /usr/bin/env python

"""
______
 -- Weronika Patena, 2013
USAGE: _____
"""

# standard library
from __future__ import division
import sys
import unittest
# other packages
# my modules
from robotic_plate_transfer import Plate_type

PLATE_PREFIXES_REMOVE = ['plate']

def mutant_number_from_plate_well(plate, well, plate_size=384):
    """ Get sequential mutant number from plate+well: first plate is 1-384, second is 385-768, etc. """
    plate_type = Plate_type(plate_size)
    # get the plate number - we don't know if it's a string or an int, strip prefixes if string, etc
    if hasattr(plate, 'startswith'):
        for prefix in PLATE_PREFIXES_REMOVE:
            if plate.lower().startswith(prefix):
                plate = plate[len(prefix):]
    plate_number = int(plate)
    # get the well number - adding 1 because the original version is 0-based
    well_number = plate_type.get_well_number_from_ID(well) + 1
    # get mutant ID
    return (plate_number-1)*plate_size + well_number
    

# TODO should probably merge this with robotic_plate_transfer.Plate_type - the unit-tests are very similar, etc...

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__mutant_number_from_plate_well(self):
        assert mutant_number_from_plate_well(1, 'A1', 384) == 1
        assert mutant_number_from_plate_well(1, 'A1', 384) == 1
        assert mutant_number_from_plate_well('plate1', 'A1', 384) == 1
        assert mutant_number_from_plate_well('plate1', 'A2', 384) == 2
        assert mutant_number_from_plate_well('plate1', 'A24', 384) == 24
        assert mutant_number_from_plate_well('plate1', 'B1', 384) == 25
        assert mutant_number_from_plate_well('plate1', 'P24', 384) == 384
        assert mutant_number_from_plate_well('plate2', 'A1', 384) == 385
        assert mutant_number_from_plate_well('plate2', 'A2', 384) == 386
        assert mutant_number_from_plate_well('plate3', 'A1', 384) == 769
        assert mutant_number_from_plate_well(1, 'A1', 96) == 1
        assert mutant_number_from_plate_well(2, 'A1', 96) == 97
        assert mutant_number_from_plate_well(1, 'A1', 6) == 1
        assert mutant_number_from_plate_well(2, 'A1', 6) == 7


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
