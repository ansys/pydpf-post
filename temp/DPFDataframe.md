# DPF Dataframe pseudocode

## Naming
Following remarks from James Derrick during the PyDPF AFT meeting of January 25:

There is no real need to identify our own dataframe class as "DPFDataframe".
The namespace should take care of conflicts not arising as it would really be
``ansys.dpf.post.dataframe.DataFrame``.
The same way as ``cudf`` has its ``cudf.DataFrame``, ``cudf.Series``... etc,
or as ``Koalas`` has ``ks.DataFrame``, ``ks.Series``...

The user should be able to do ``df = dpf.DataFrame()``.


## Definition and comparison with existing dataframes
### From Pandas 
"Two-dimensional, size-mutable, potentially heterogeneous tabular data."

### From cuDF 
"cuDF is a Python GPU DataFrame library (built on the Apache Arrow columnar memory
format) for loading, joining, aggregating, filtering, and otherwise manipulating tabular data
using a DataFrame style API."

### In our case
A DPF DataFrame offers a DataFrame style API for manipulating DPF data.
It differs from other dataframes as it differentiates between spatial
and time-like/set (sequence?) dimensions to be able to link data to the underlying simulation.

A value is by default indexed using the mesh entity ID, a single or multi-index 
(in the case of ElementalNodal data).
The set axis for time/set/frequency ID, is by default the main index of the columns,
with results and then components as sub-indexes.
Index types could be used to differentiate them (see 
[here](https://pandas.pydata.org/docs/user_guide/advanced.html#index-types)).
These could be ``SpatialIndex`` and ``TimeIndex``.


A DPF ``Dataframe`` wraps a ``FieldsContainer`` (or several?). Operations on the data, and thus on 
the ``FieldsContainer`` should be handled by DPF operators and workflows server-side whenever 
possible. When working with gRPC however, this would mean going back-and-forth, so calls to the 
server should be minimized.

## Create from scratch

Instantiating a dpf.Dataframe should mimic instantiation of a pandas.DataFrame.
```py
import ansys.dpf.post as dpf
df = dpf.DataFrame(
    data=data,
    index=index,
    columns=columns,
    copy=True,
)
```

with:
- ``data`` either a FieldsContainer, a DPFArray?, a numpy.ndarray, a pandas.DataFrame. dict?
- ``index`` a multi-index with spatial and timefreq IDs, used when 'data' is not a DPF object.
- ``columns`` a list of column labels to use if 'data' does not define them.
- ``copy`` a boolean or None, whether to copy data from inputs.

We most likely want a ``ansys.dpf.post.read_csv`` which instantiates a DataFrame

```py
df = dpf.read_csv(filepath, sep=';', delimiter=None, header='infer', decimal=',', ...)
```

## Create from a Simulation

Querying a result from a ``Simulation`` returns a ``DataFrame``. 
```py
import ansys.dpf.post as dpf
from ansys.dpf.post import examples

simulation = dpf.load_simulation(examples.static_rst)
df = simulation.stress_principal_nodal()
```


## View Data

Let's have a dataframe ``df``:
```pycon
>>> df
    step|     1                     |     2                     |
  result|     S      E      U       |     S      E      U       |
    comp|    S1     E2     UX     UY|    S1     E2     UX     UY|
 node ID
       1|   0.4    0.2    0.3    0.3|   0.5    0.3    0.3    0.4|
       2|   0.3    0.3    0.4    0.3|   0.5    0.3    0.3    0.4|
       3|   0.2    0.4    0.5    0.3|   0.5    0.3    0.3    0.4|
       4|   0.1    0.5    0.6    0.3|   0.5    0.3    0.3    0.4|

[ 4 rows x 8 (2 steps x 4 results) columns]
```

### Large arrays

See [here](https://pandas.pydata.org/docs/user_guide/dsintro.html#console-display)

Options exist in Pandas to configure the console display:
- "pd.set_option("display.width", 40)  # default is 80"
- "pd.set_option("display.max_colwidth", 100)"

Printing ``df`` would print the 2 first and 2 last "rows" by default, 
with the 3 first and 3 last "columns", or as much as can fit in the current ``display.width`` 
using the current ``display.max_colwidth``.
For example, for a big DataFrame with ``display_width=150`` (effective is 121):
```pycon
>>> df
    step>     1                                                 |     2                                                 
  result>     U            |     S                              |     U            |     S                              
    comp>    UX    UY    UZ|   SXX   SYY   SZZ   SXY   SYZ   SXZ|    UX    UY    UZ|   SXX   SYY   SZZ   SXY   SYZ   SXZ
 node ID.
       1|   0.0   0.1   0.2|   0.3   0.4   0.5   0.3   0.4   0.5|   0.0   0.1   0.2|   0.3   0.4   0.5   0.3   0.4   0.5
       2|   0.3   0.1   0.1|   0.3   0.4   0.5   0.3   0.4   0.5|   0.3   0.1   0.1|   0.3   0.4   0.5   0.3   0.4   0.5
     ...|   ...   ...   ...|   ...   ...   ...   ...   ...   ...|   ...   ...   ...|   ...   ...   ...   ...   ...   ...
      99|   0.3   0.3   0.3|   0.3   0.3   0.3   0.3   0.3   0.3|   0.3   0.3   0.3|   0.3   0.3   0.3   0.3   0.3   0.3
     100|   0.2   0.3   0.3|   0.3   0.3   0.3   0.3   0.3   0.3|   0.2   0.3   0.3|   0.3   0.3   0.3   0.3   0.3   0.3

[ 100 rows x 9 (1 step x 9 results) columns ]
```
With ``display_width=80`` (effective is 77):
```pycon
>>> df
    step>     1                            ...     2                         
  result>     U            |     S         ...     S               
    comp>    UX    UY    UZ|   SXX   SYY   ...   SYY   SZZ   SXY   SYZ   SXZ
 node ID.
       1|   0.0   0.1   0.2|   0.3   0.4   ...   0.5   0.5   0.3   0.4   0.5
       2|   0.3   0.1   0.1|   0.3   0.4   ...   0.5   0.5   0.3   0.4   0.5
     ...|   ...   ...   ...|   ...   ...   ...   ...   ...   ...   ...   ...
      99|   0.3   0.3   0.3|   0.3   0.3   ...   0.5   0.3   0.3   0.3   0.3
     100|   0.2   0.3   0.3|   0.3   0.3   ...   0.5   0.3   0.3   0.3   0.3

[ 100 rows x 9 (1 step x 9 results) columns ]
```

#### Memory usage
Pandas offers a useful method for its DataFrame called ``DataFrame.memory_usage``. 
See [here](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.memory_usage.html#pandas-dataframe-memory-usage).
Having the same would be useful, especially to assess and/or explain performance.
It should also be printed-out in the ``DataFrame.info`` method.

### Head and tail 
Viewing the first rows.
A DPF DataFrame is sorted by the time/freq value by default.
The mesh entity IDs are indeed not necessarily sorted:
```pycon
>>> df.head(2)
  step  node|    S1    S2    S3
     1     3|   0.2   0.4   0.5
     1     2|   0.3   0.3   0.4
```

Viewing the last rows:
```pycon
>>> df.tail(2)
  step  node|    S1    S2    S3
     2     3|   0.8   0.4   0.5
     2     2|   0.9   0.3   0.4
```

### Sorting
Sorting by mesh entity ID:
```pycon
>>> df.sort_values(by='node', ascending=True)
  step  node|    S1    S2    S3
     1     1|   0.4   0.2   0.3
     1     2|   0.3   0.3   0.4
     1     3|   0.2   0.4   0.5
     1     4|   0.1   0.5   0.6
```

Sorting values by column value:
```pycon
>>> df.sort_values(by='S1', ascending=True, )
  step  node|    S1    S2    S3
     2     4|   0.0   0.0   0.1
     1     4|   0.1   0.5   0.6
     1     3|   0.2   0.4   0.5
     1     2|   0.3   0.3   0.4
```

### Index and columns
Display the ``DataFrame.index`` or ``DataFrame.columns``:
```pycon
>>> df.index
{"spatial": Index([1, 2, 3, 4, 5], dtype='node'), 
"time": Index([1, 2], dtype='step')}
>>> df.columns
Index(['S1', 'S2', 'S2'], dtype=[<stress>, <stress>, <stress>], units=['Pa', 'Pa', 'Pa'])
```

Name the 'time' index depending on the simulation type:
- static: set ID
- transient: time-step ID
- modal: mode ID
- harmonic: (frequency, phase) set ID?

### Get info on the DataFrame

The length of the whole DataFrame in Pandas returns the "number of rows".

For DPF however, returning the number of rows has no meaning since Fields can have different lengths.
Thus ``len()`` should work on axes, not on a full DataFrame.
Or we can output the number of "sets IDs * mesh entity IDs".
To get info on a full DataFrame:
 - ``df.info``
 - ``df.shape``
 - ``df.columns``
 - ``df.axes``
 - ``len(df.axes[i])``

Since the notion of ``Index`` is probably required anyway (exposing the DPF ``Scoping``),
then the ``__len__`` dunder method can be defined at this level. 
```pycon
>>> len(df)
Error: a DPF DataFrame has no length as it 
>>> df.columns
Index(['S1', 'S2', 'S2'], dtype=[<stress>, <stress>, <stress>], units=['Pa', 'Pa', 'Pa'])
```

## Select Data

Note appearing in pandas documentation:

```
While standard Python / NumPy expressions for selecting and setting are intuitive and come in handy 
for interactive work, for production code, we recommend the optimized pandas data access methods, 
DataFrame.at(), DataFrame.iat(), DataFrame.loc() and DataFrame.iloc().
```
We, too, could provide ``DataFrame.at()``, ``DataFrame.iat()``, ``DataFrame.loc()`` and 
``DataFrame.iloc()`` methods, while overriding regular Python dunder methods (``__get_item__``, ...).

Should priority be given to optimized ``DataFrame.at()``, ``DataFrame.iat()``, ``DataFrame.loc()`` and 
``DataFrame.iloc()`` methods?

### Get Data

Selecting a single column yields a one-dimensional ``DataFrame`` (no ``Series``).
Could be equivalent to ``df.S1``.
```pycon
>>> df["S1"]
  step  node|    S1
     1     1|   0.4
     1     2|   0.3
     1     3|   0.2
     1     4|   0.1
dtype: <stress>, Unit: 'Pa'
```

Selecting via ``[]`` (``__get_item__``) slices the rows index by index.
--> Complicated for multi-index. We would diverge from Pandas as slicing a multi-indexed 
``pandas.DataFrame`` requires explicit use of either ``DataFrame.loc`` or ``DataFrame.select``. 
See [here](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.loc.html#pandas-dataframe-loc)
```pycon
>>> df.iloc[[0, -1], [:, S1]]
  step  node|    S1    S2    S3
     1     1|   0.4   0.2   0.3
     1     2|   0.3   0.3   0.4
```

### Select by label (ID)
Extract ``S2`` and ``S3``:
```pycon
>>> df.loc[[:, :], ["S2", "S3"]]
  step  node|    S2    S3
     1     1|   0.2   0.3
     1     2|   0.3   0.4
     1     3|   0.4   0.5
     1     4|   0.5   0.6
```

Get a scalar value:
```pycon
>>> df.loc[[1, 1], "S2"]
0.2
```
And its equivalent:
```pycon
>>> df.at[[1, 1], "S2"]
0.2
```

### Select by position (index)
Use ``DataFrame.iloc()`` or ``DataFrame.iat()`` ("fast-access version").

Select via the position of the passed integers:
```pycon
>>> df.iloc[[1], [1, [1:2]]]
    S1   0.4   
    S2   0.2
set: 1, node: 1 
```

Get a value explicitly:
```pycon
>>> df.iloc[[1, 1], 1]
0.4
```

### Boolean indexing

One of the strengths of the ``pandas.DataFrame`` is the possibility to select using conditions:
```pycon
>>> df[df['S1']>0.2]
  step  node|    S1    S2    S3
     1     1|   0.4   0.2   0.3
     1     2|   0.3   0.3   0.4
```
We should most likely restrict this to row filtering.

Indeed, selecting sparse values from a ``pandas.DataFrame`` where a boolean condition is met 
requires using ``np.nan`` to get a Pandas-like behavior:
```pycon
>>> df[df>0.3]
  step  node|    S1    S2    S3
     1     1|   0.4   NaN   NaN
     1     2|   NaN   NaN   0.4
     1     3|   NaN   0.4   0.5
     1     4|   NaN   0.5   0.6
```
This would require us to deal with ``np.nan`` everywhere, which we said is a bad idea.

Not having ``np.nan`` thus means we cannot filter by value on a multidimensional dataframe.
Only filtering whole rows is feasible.

## Merge, join, concatenate and compare

We need to propose functionalities to combine or compare different dataframes.

``dpf.DataFrame.concat``

``dpf.DataFrame.merge``

``dpf.DataFrame.compare``


## Operations

### Stats

IMO we should propose basic functions, doable with DPF, but not try to compete with Pandas.

See [here](https://pandas.pydata.org/docs/user_guide/basics.html#descriptive-statistics).

The method ``df.mean`` would allow for arguments such as for ``pandas.DataFrame.mean``.
The first argument is ``axis``, which defines which axis to apply the operation on 
(default ``axis=0``, meaning average along the rows/for each column.

The default behavior would be to apply the mean on all present set and mesh IDs.
```pycon
>>> df.mean(axis=0)
    S1   0.3   
    S2   0.4
    S3   0.5
```
The type of this result if to be defined. Do we need a ``dpf.Series``, or is a ``dict`` sufficient, 
or a tuple?

To get the mean for a specific subset, then a selection step is sufficient.
To get the average on all nodes for each column but only on step 1: 
```pycon
>>> df.loc[[1, :]].mean(axis=0)
    S1  0.25   
    S2  0.35
    S3  0.45
```

Applying the mean on ``axis=1`` would here get us average for each row:
```pycon
>>> df.mean(axis=1)
  step  node|  
     1     1|   0.3
     1     2|0.3333
     1     3|0.3667
     1     4|   0.4
     2     1|   0.4
     2     2|0.4333
     2     3|0.4667
     2     4|   0.5
```

### Apply

Applying a function on either ``axis``:
```pycon
>>> df.apply(np.mean, axis=1)
  step  node|  
     1     1|   0.3
     1     2|0.3333
     1     3|0.3667
     1     4|   0.4
     2     1|   0.4
     2     2|0.4333
     2     3|0.4667
     2     4|   0.5
```
Is that doable with DPF?

## Plotting

See [here](https://pandas.pydata.org/docs/user_guide/10min.html#plotting).
There is no intrinsic notion of 3D plotting with ``pandas.DataFrame.plot``.
Its default backend is ``matplotlib`` and it is only a helper which exposes usual arguments.

A possibility is to reuse the ``kind`` argument found in ``pandas.DataFrame.plot`` to include a 
``contour`` option, which could even be the default, and which would switch the backend depending 
on the kind of plot requested.
See [here](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.html#pandas-dataframe-plot).

```py
df.plot(
    data,
    x:
    y:
    kind:
)
```

## I/O

```py
df.to_pandas
df.to_ndarray
```