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
and time-like/set dimensions to be able to link data to the underlying simulation.

A value is always indexed using a multi-index composed of at least one spatial entity ID and one
time/set/frequency ID

There will always be a spatial index axis, consisting of one or a combination of 1-based IDs:
node ID, element ID, part ID, 'global', layer ID...
The underlying equivalent PyDPF-Core location being global, nodal, elemental, or elemental nodal.

A DPF ``Dataframe`` wraps a ``FieldsContainer`` (or several?). Operations on the data, and thus on the 
``FieldsContainer`` should be handled by DPF operators and workflows server-side whenever possible.
When working with gRPC however, this would mean going back-and-forth...

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

Sorting values by mesh entity ID:
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
>>> df.sort_values(by='S1', ascending=True)
  step  node|    S1    S2    S3
     2     4|   0.0   0.0   0.1
     1     4|   0.1   0.5   0.6
     1     3|   0.2   0.4   0.5
     1     2|   0.3   0.3   0.4
```

Display the ``DataFrame.index`` or ``DataFrame.columns``:
```pycon
>>> df.index
{"spatial": Index([1, 2, 3, 4, 5], dtype='node'), 
"event": Index([1, 2], dtype='step')}
>>> df.columns
Index(['S1', 'S2', 'S2'], dtype=[<stress>, <stress>, <stress>], units=['Pa', 'Pa', 'Pa'])
```

"Time-like" index naming ideas:
time, event, progress, progression, set, configuration

Or name it depending on the simulation type.

## Select Data

## Retrieve Data

## Interactions between dataframes

## Statistics

## Plotting

## I/O

```py
df.to_pandas
df.to_array
```