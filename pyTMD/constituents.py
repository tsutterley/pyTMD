#!/usr/bin/env python
u"""
constituents.py
Written by Tyler Sutterley (12/2022)
Basic tide model constituent class

UPDATE HISTORY:
    Written 12/2022
"""

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
        # set optional fields
        for key, val in kwargs.items():
            setattr(self, key, val)

    def append(self, field, constituent):
        """
        Append tide model constituents

        Parameters
        ----------
        field: str
            Tide model constituent name
        constituent: float
            Tide model constituent (complex form)
        """
        # append field
        self.fields.append(field)
        setattr(self, field, constituent)
        return self

    def get(self, field):
        """
        Get model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name

        Returns
        -------
        constituent: float
            Tide model constituent (complex form)
        """
        return getattr(self, field)
