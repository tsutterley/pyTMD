#!/usr/bin/env python
u"""
constituents.py
Written by Tyler Sutterley (08/2024)
Basic tide model constituent class

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 08/2024: add GOT prime nomenclature for 3rd degree constituents
    Updated 07/2024: add function to parse tidal constituents from strings
    Updated 05/2024: make subscriptable and allow item assignment
    Updated 01/2024: added properties for Doodson and Cartwright numbers
    Updated 08/2023: added default for printing constituent class
    Updated 07/2023: output constituent from get and pop as copy
    Updated 03/2023: add basic variable typing to function inputs
    Written 12/2022
"""
from __future__ import division, annotations

import re
import copy
import numpy as np
import pyTMD.arguments

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

    @property
    def doodson_number(self):
        """constituent Doodson number
        """
        doodson_numbers = []
        # for each constituent ID
        for f in self.fields:
            try:
                # try to get the Doodson number
                n = pyTMD.arguments.doodson_number(f)
            except (AssertionError, ValueError) as exc:
                n = None
            # add Doodson number to the combined list
            doodson_numbers.append(n)
        # return the list of Doodson numbers
        return doodson_numbers

    @property
    def cartwright_number(self):
        """constituent Cartwright numbers
        """
        cartwright_numbers = []
        # for each constituent ID
        for f in self.fields:
            try:
                # try to get the Cartwright numbers
                n = pyTMD.arguments.doodson_number(f, formalism='Cartwright')
            except (AssertionError, ValueError) as exc:
                n = None
            # add Cartwright numbers to the combined list
            cartwright_numbers.append(n)
        # return the list of Cartwright numbers
        return cartwright_numbers

    @staticmethod
    def parse(constituent: str) -> str:
        """
        Parses for tidal constituents using regular expressions and
        remapping of known cases

        Parameters
        ----------
        constituent: str
            Unparsed tidal constituent name
        """
        # list of tidal constituents (not all are included in tidal program)
        # include negative look-behind and look-ahead for complex cases
        cindex = [r'(?<!s)sa','ssa','mm','msf',r'mt(?!m)(?!ide)','mf','alpha1',
            '2q1','sigma1',r'(?<!2)q1','rho1',r'(?<!rh)(?<!o)o1','tau1',
            'm1','chi1','pi1','p1','s1','k1','psi1','phi1','theta1','j1',
            'oo1','2n2','mu2',r'(?<!2)n2','nu2',r'(?<!2s)m2(?!a)(?!b)',
            'm2a','m2b','lambda2','l2','t2',r'(?<!mn)(?<!mk)(?<!ep)s2(?!0)',
            'r2','k2','eta2','mns2','2sm2','m3','mk3','s3','mn4','m4',
            'ms4','mk4',r'(?<!m)s4','s5','m6','s6','s7','s8','m8','mks2',
            'msqm','mtm',r'(?<!m)n4','eps2','z0']
        # compile regular expression
        # adding GOT prime nomenclature for 3rd degree constituents
        rx = re.compile(r'(' + '|'.join(cindex) + r')(\')?', re.IGNORECASE)
        # check if tide model is a simple regex case
        if rx.search(constituent):
            return "".join(rx.findall(constituent)[0]).lower()
        # known remapped cases
        mapping = [('2n','2n2'), ('e2','eps2'), ('la2','lambda2'),
            ('sig1','sigma1')]
        # iterate over known remapped cases
        for m in mapping:
            # check if tide model is a remapped case
            if m[0] in constituent.lower():
                return m[1]
        # raise a value error if not found
        raise ValueError(f'Constituent not found in {constituent}')

    def __str__(self):
        """String representation of the ``constituents`` object
        """
        properties = ['pyTMD.constituents']
        fields = ', '.join(self.fields)
        properties.append(f"    constituents: {fields}")
        return '\n'.join(properties)

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

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
