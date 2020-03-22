convert_calendar_decimal.py
===========================

 - Converts from calendar date into decimal years taking into account leap years  

#### Calling Sequence
```
from gravity_toolkit.convert_calendar_decimal import convert_calendar_decimal
t_date = convert_calendar_decimal(year, month, DAY=day, \
    HOUR=hour, MINUTE=minute, SECOND=second)
```

#### Inputs
 1. `year`: calendar year  
 2. `month`: calendar month (1 = January, ..., 12 = December)  

#### Options
 - `DAY`: Number of day of the month
 - `HOUR`: hour of the day
 - `MINUTE`: minute of the hour
 - `SECOND`: second (and fractions of a second) of the minute
 - `DofY`: day of the year (January 1 = 1)

#### Outputs
 - `t_date`: date in decimal format (years)
