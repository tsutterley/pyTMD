#!/usr/bin/env python
u"""
constituents.py
Written by Tyler Sutterley (07/2023)
Basic tide model constituent class

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2023: output constituent from get and pop as copy
    Updated 03/2023: add basic variable typing to function inputs
    Written 12/2022
"""
from __future__ import division, annotations

import copy
import numpy as np

class constituents:
    """
    Class for tide model constituents

    Attributes
    ----------
    fields: list
        list of tide model constituents
    """
    def __init__(self, **kwargs):
        # set initial attributes
        self.fields = []
        self.__index__ = 0
        # set optional fields
        for key, val in kwargs.items():
            setattr(self, key, val)

    def append(self, field: str, constituent: np.ndarray):
        """
        Append a tide model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name
        constituent: np.ndarray
            Tide model constituent (complex form)
        """
        # append field
        self.fields.append(field)
        setattr(self, field, constituent)
        return self

    def get(self, field: str):
        """
        Get a tide model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name

        Returns
        -------
        constituent: np.ndarray
            Tide model constituent (complex form)
        """
        constituent = getattr(self, field)
        return copy.copy(constituent)

    def pop(self, field: str):
        """
        Retrieve a tide model constituent and remove from list

        Parameters
        ----------
        field: str
            Tide model constituent name

        Returns
        -------
        constituent: np.ndarray
            Tide model constituent (complex form)
        """
        self.fields.remove(field)
        constituent = getattr(self, field)
        delattr(self, field)
        return copy.copy(constituent)

    def update(self, field: str, constituent: np.ndarray):
        """
        Update a tide model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name
        constituent: np.ndarray
            Tide model constituent (complex form)
        """
        # raise exception if field not in list
        if not hasattr(self, field):
            raise KeyError(f'Constituent {field}')
        # update the constituent
        setattr(self, field, constituent)
        return self

    def amplitude(self, field: str):
        """
        Calculate the amplitude of a tide model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name

        Returns
        -------
        amp: np.ndarray
            Tide model constituent amplitude
        """
        constituent = getattr(self, field)
        # calculate constituent amplitude
        amp = np.sqrt(constituent.real**2 + constituent.imag**2)
        # update mask and fill values
        amp.mask = np.copy(constituent.mask)
        amp.data[amp.mask] = amp.fill_value
        return amp

    def phase(self, field: str):
        """
        Calculate the phase of a tide model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name

        Returns
        -------
        ph: float
            Tide model constituent phase (degrees)
        """
        constituent = getattr(self, field)
        # calculate constituent phase and convert to degrees
        ph = 180.0*np.arctan2(-constituent.imag, constituent.real)/np.pi
        ph.data[ph.data < 0] += 360.0
        # update mask and fill values
        ph.mask = np.copy(constituent.mask)
        ph.data[ph.mask] = ph.fill_value
        return ph

    def __len__(self):
        """Number of constituents
        """
        return len(self.fields)

    def __iter__(self):
        """Iterate over constituents
        """
        self.__index__ = 0
        return self

    def __next__(self):
        """Get the next constituent
        """
        try:
            field = self.fields[self.__index__]
        except IndexError as exc:
            raise StopIteration from exc
        # get the model constituent
        constituent = getattr(self, field)
        self.__index__ += 1
        return (field, constituent)
