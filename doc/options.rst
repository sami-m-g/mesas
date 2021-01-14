Options
=======
A number of optional parameters can be set on a ``Model``. These can be set as keywords when the instance is created::

    my_model = Model(..., option1='value', option2='another')

or by assigning a ``dict`` of valid key-value pairs::

    my_model.options = {option1='value', option2='another'}

In the latter case, only the options specifically referenced will be changed. The others will retain their previous values (i.e. the defaults unless they have been previously changed).

Default options are

'dt': 1
'verbose': False
'debug': False
'warning': True
'jacobian': False
'n_substeps': 1
'max_age': None
'sT_init': None
'influx': 'J'
'ST_smallest_segment': 1./100
'ST_largest_segment': 1000000.

