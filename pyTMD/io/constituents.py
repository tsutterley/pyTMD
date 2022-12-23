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
        self.__index__ = 0
        # set optional fields
        for key, val in kwargs.items():
            setattr(self, key, val)

    def append(self, field, constituent):
        """
        Append a tide model constituent

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
        Get a tide model constituent

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

    def pop(self, field):
        """
        Retrieve a tide model constituent and remove from list

        Parameters
        ----------
        field: str
            Tide model constituent name

        Returns
        -------
        constituent: float
            Tide model constituent (complex form)
        """
        self.fields.remove(field)
        constituent = getattr(self, field)
        delattr(self, field)
        return constituent

    def update(self, field, constituent):
        """
        Update a tide model constituent

        Parameters
        ----------
        field: str
            Tide model constituent name
        constituent: float
            Tide model constituent (complex form)
        """
        # raise exception if field not in list
        if not hasattr(self, field):
            raise KeyError(f'Constituent {field}')
        # update the constituent
        setattr(self, field, constituent)
        return self

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
        except IndexError as esc:
            raise StopIteration from esc
        # get the model constituent
        constituent = getattr(self, field)
        self.__index__ += 1
        return (field, constituent)
