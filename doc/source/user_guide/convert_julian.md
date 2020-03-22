convert_julian.py
=================

 - Return the calendar date and time given Julian date

#### Calling Sequence
```
from gravity_toolkit.convert_julian import convert_julian
YEAR,MONTH,DAY,HOUR,MINUTE,SECOND = convert_julian(JD, FORMAT='tuple')
```

#### Inputs
 1. `JD`: Julian Day of the specified calendar date (days since -4712-01-01T12:00:00)  

#### Options
 - `ASTYPE`: convert output to variable type  
 - `FORMAT`: format of output variables  
    - `'dict'`: dictionary with variable keys as listed above  
    - `'tuple'`: tuple with variable order YEAR,MONTH,DAY,HOUR,MINUTE,SECOND  
    - `'zip'`: aggregated variable sets  

#### Outputs
 - `year`: Number of the desired year
 - `month`: Number of the desired month (1 = January, ..., 12 = December)
 - `day`: Number of day of the month
 - `hour`: hour of the day
 - `minute`: minute of the hour
 - `second`: second (and fractions of a second) of the minute
