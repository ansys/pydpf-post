# DPF Refactoring: Static Mechanical Analysis

Pseudocode for StaticMechanicalSimulation

## Load the results

Load the result files within an instance of `StaticMechanicalSimulation`

The `ansys.dpf.post.load_simulation` function detects by default the type physics
and the type of analysis based on the metadata in the underlying `ansys.dpf.core.Model` instance.
If no physics type or analysis type is recognized automatically, 
it will default to static and mechanical.
You can specify the physics type and analysis type. 

```pycon
>>> static_simulation = dpf.load_simulation(examples.simple_bar)
>>> static_simulation = dpf.load_simulation(examples.simple_bar,
...                                        physics_type=physics_types.mechanical,
...                                        analysis_type=analysis_types.static)
```

To that end, we might want to provide Public Enums for available physics_types and analysis_types
instead of having them private.
```pycon
>>> from ansys.dpf.post import physics_types, analysis_types
>>> # instead of
>>> # from ansys.dpf.post.common import _AnalysisType, _PhysicsType
```

## Explore the simulation object

Print information about the simulation
```pycon
>>> print(static_simulation)
TBD
```

Print the mesh
```pycon
>>> print(static_simulation.mesh)
TBD
```

Print the list of constructed geometries
```pycon
>>> print(static_simulation.geometries)
TBD
```

Print the list of boundary conditions
```pycon
>>> print(static_simulation.boundary_conditions)
TBD
```

Print the list of loads
```pycon
>>> print(static_simulation.loads)
TBD
```

Print the list of available named selections
```pycon
>>> print(static_simulation.named_selections)
TBD
```

Print the list of steps (the time_freq_support)? the list of sub-steps? both? the times?
```pycon
>>> print(static_simulation.time_freq_support)
TBD
>>> print(static_simulation.steps)
TBD
```

Print the list of available results
```pycon
>>> print(static_simulation.results)
TBD
```

General plot of the simulation object with, by default:
- the mesh, not deformed, at step 0
- the geometry
- the boundary conditions
```pycon
>>> static_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)
Image
```

## Create and use geometry

Creation of geometry is for now in PyDPF-Core.
Do we expose those in Post in some way? Simply by importing them in a specific Post module?

### Add geometry to the simulation
```pycon
>>> from ansys.dpf.core.geometry import Line
>>> line = Line(coordinates=[[1.0, 0.0, 1.0], [5.0, 0.0, 1.0]])
>>> static_simulation.geometries.append(line)
>>> static_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)
Image
```

## Create and use boundary conditions

TBD

## Create and use loads

Loads are not read from files yet.

### Create a load
A load would then be derived from a DataObject, which has a mesh support, thus:
```pycon
>>> from ansys.dpf.post.load import Load 
>>> load = Load(data=data_object)
>>> print(load)
TBD
```

### Map loads from one mesh to another 
(which is the same as DataObject.interpolate?)
```pycon
>>> load.interpolate(mesh=another_mesh)
```

This means a load has a support, either nodal, elemental or face.
This means `location.face` must be available.
It can also vary with step/time/frequency/rpm or any dimension present in the DataObject.

The metadata and functionalities unique to a Load would be:
- its naming
- its string representation
- its graphic representation

## Create and use selections

### Create a selection
```pycon
>>> selection = Selection()
>>> selection.select_nodes(nodes=[1, 2, 3])
>>> selection.select_elements(elements=[1, 2, 3])
>>> selection.select_time_freq_indices(time_freq_indices=[1])  # Rework?
>>> # selection.select_steps(steps=[1])
>>> # selection.select_geometry(geometry=Line())
>>> selection.select_named_selection(named_selection=named_selections[0])
```

### Retrieve the spatial part of a selection
```pycon
>>> print(selection.spatial_selection)
TBD
```

### Intersect two spatial selections
```pycon
>>> space_sel_1 = selection_1.spatial_selection
>>> print(space_sel_1)
Spatial selection with 3 node(s):
[1, 2, 3]
>>> space_sel_2 = selection_2.spatial_selection
>>> print(space_sel_2)
Spatial selection with 3 node(s):
[3, 4, 5]
>>> space_sel_1.intersect(space_sel_2)
>>> print(space_sel_1)
Spatial selection with 1 node(s):
[3]
```

### Union of two spatial selections
```pycon
>>> selection_3 = Selection()
>>> selection_3.select_nodes(nodes=[1, 2, 3])
>>> selection_4 = Selection()
>>> selection_4.select_nodes(nodes=[3, 4, 5])
>>> space_sel_3 = selection_3.spatial_selection
>>> space_sel_4 = selection_4.spatial_selection
>>> space_sel_3.union(space_sel_4)
>>> print(space_sel_3)
Spatial selection with 5 node(s):
[1, 2, 3, 4, 5]
```
Do we then need to explicitly update the spatial_selection of selection_1?

### Select based on constructed geometry
```pycon
>>> from ansys.dpf.core.geometry import Line
>>> line = Line(coordinates=[[1.0, 0.0, 1.0], [5.0, 0.0, 1.0]])
>>> selection_5 = Selection()
>>> selection_5.select_geometry(line)
>>> print(selection_5)
Selection with:
- 0 steps

- 100 nodes 
[...]
```
Right now a Line element has an underlying mesh with nodes and elements. We could print out this information.

### Set a selection as default 
for a simulation for further result extractions
```pycon
>>> static_simulation.activate_selection(selection=selection)
>>> print(static_simulation.active_selection)
Selection with:
- 0 steps

- 0 node/elements
 
```

### Remove a selection as default 
for a simulation for further result extractions
```pycon
>>> static_simulation.deactivate_selection()
>>> print(static_simulation.active_selection)
None 
```

### Apply a SpatialSelection to a Simulation 
to get resulting list of entities IDs (nodes or elements)
```pycon
>>> list_of_IDs = selection.spatial_selection.apply_to(static_simulation)
>>> print(list_of_IDs)
[ ... ]
```

## Mesh operations

### Extract the mesh
```pycon
>>> mesh = static_simulation.mesh
>>> print(mesh)
DPF Mesh:
    XXXX nodes
    XXXX elements
    Unit: m
    with solid (3D) elements
```

### Plot the mesh 
with defaults for static:
- not deformed
- at last step
```pycon
>>> mesh.plot(
...     opacity=0.3, title="Static mesh plot", text="defaults to not deformed, last step"
... )
<Image>
```

### Plot elements orientations with triads
```pycon
>>> mesh.plot(triads=True)
<Image>
```

### Map a load onto the mesh
A load is a DataObject, so has a mesh support. Thus:
```pycon
>>> load_2 = mesh.interpolate(load)
```

## Extract specific results

### Extract displacements along X for nodes 1, 2 and 3 at step 1
```pycon
>>> displacement_X = static_simulation.displacement(
...     component="X", nodes=[1, 2, 3], steps=[1]
... )
```

### Extract norm of displacements for nodes 1, 2 and 3 at step 1
```pycon
>>> displacement_norm = static_simulation.displacement(
...     component="N", nodes=[1, 2, 3], steps=[1]
... )
```

### Extract nodal XY stresses for elements 1, 2 and 3 at step 1
```pycon
>>> stress_XY = static_simulation.elemental_stress(
...     component="XY", elements=[1, 2, 3], steps=[1]
... )
```

### Extract first principal nodal stress for a named (elemental or nodal) selection at all steps
```pycon
>>> stress_S1 = static_simulation.nodal_stress(
...     component="S1", named_selection=named_selections[0]
... )
```

### Extract elemental Von Mises stress everywhere at all steps
```pycon
>>> stress_VM = static_simulation.elemental_stress(component="VM")
```

### Extract equivalent elemental nodal elastic strain for a selection at step 1
```pycon
>>> elastic_strain_XY = static_simulation.nodal_elastic_strain(
...     component="XY", selection=selection, steps=[1]
... )
```

### Extract first principal nodal strain for a selection at step 1
```pycon
>>> elastic_strain_E1 = static_simulation.nodal_elastic_strain(
...     component="E1", selection=selection, steps=[1]
... )
```

### Extract nodal plastic strain for a selection at step 1
```pycon
>>> plastic_strain = static_simulation.nodal_plastic_strain(selection=selection, steps=[1])
```

## Manipulate results with DataObject

### Explore a DataObject
Print a DataObject

```pycon
>>> print(stress_S1)
DPF DataObject:
Nodal stress S1
    on:
    - named selection {named_selection_name} # if present
    - N nodes/elements/faces:
        [1, ..., 45]  # summarized list representation if less than X characters long
    - N steps/frequencies/times:
        [...]
    with:
    - S1 stress (Pa)
        [range]
```
Or print-out the DataObject in a dataframe style,
with an index on the main supporting mesh entity id
```pycon
DPF DataObject:
Nodal stress S1 [on *named_selection_name*] # if present
 node ID   step  component   value
       1      1         S1     0.1
       2      1         S1     0.2
       3      1         S1     0.3
       ...
       1      2         S1     0.2
       2      2         S1     0.3
       3      2         S1     NaN
```

Get the name of the DataObject
print(stress_S1.name)
Nodal stress S1

Get the mesh support of the DataObject
print(stress_S1.mesh)

Get the shape of the DataObject
```pycon
>>> print(stress_S1.shape)
list of lengths of each dimension present: [nb_steps, nb_entities,
```

Get specific data from a DataObject (as a new DataObject)
```pycon
>>> stress_S1_step_1 = stress_S1[1]  # Equivalent to stress_S1[0, ...] or stress_S1[0, :, :]
```
or
```pycon
>>> stress_S1_step_1 = stress_S1.step(1)
```

```pycon
>>> print(stress_S1_step_1)
Nodal stress S1
For:
- steps: 1
- named selection: {named_selection[0]}
```

### Get the minimum, maximum, mean of the DataObject 
Use DPF operators for efficiency.

Across all dimensions present in the specific DataObject, unless specified as argument
Take inspiration from the "axis" argument for ``numpy.ndarray.amax``.
Do we want a ``nanmax``? Which propagates NaN values? (Whether we have actual NaNs or not)
These would return DataObjects or a scalar.
```pycon
>>> global_max_stress_S1 = stress_S1.max()
```
```pycon
>>> max_stress_S1_step_1 = stress_S1_step_1.max()
```
or
```pycon
max_stress_S1_step_1 = stress_S1.max(steps=[1])
```

Arguments would include lists of steps, nodes, elements, components... a named selection
```pycon
>>> max_stress_S1_node_1 = stress_S1.max(nodes=[1])
```

Get the mesh support for a specific DataObject (uses the scoping associated to the DataObject)
```pycon
>>> mesh_1 = max_stress_S1_node_1.mesh
>>> print(mesh_1)
Mesh containing:
- 1 nodes: [1]
- 0 elements: []
```

### Plot a DataObject (3D contour plots)
TODO:
plotting and graphing functions can return the Plotter instance

Plot a multistep DataObject - takes an optional "step" argument, defaults to last
```pycon
>>> pl = stress_S1.plot(step=3, return_plotter=True)
```

Plot a one-step DataObject - the "step" argument is not read
```pycon
>>> stress_S1_step_1.plot()
```

Plot principal stresses

TODO: Not really sure how this is different from plotting results


### Create graphs (2D curve plots)
TODO:
the graph method would take optional "x", "y", "names" arguments.

Graph a multistep DataObject
```pycon
>>> max_stress_S1_node_1.graph(x="steps", names="ID")
```
or as proposed in the GitHub discussion #213
```pycon
>>> max_stress_S1_node_1.plot(graph=True, individual_figures=False)
```

Plot waterfall/stagger diagrams -> a waterfall diagram is different from a waterfall plot!

Plot waterfall plots


### Animate a DataObject

Animate a contour in time
```pycon
>>> stress_S1.animate(axis="steps")
```

### Animated graphs

Animate a curve in time
```pycon
>>> max_stress_S1_node_1.animate_graph(axis="steps", x="ID")
```

### Interpolate data from one mesh to another

Starts from the mesh support of the DataObject

Performs interpolation for each step and for each component present by default (for each field)
```pycon
>>> stress_S1_on_another_mesh = stress_S1.interpolate(mesh=another_mesh)
```

### Export to a numpy ndarray
```pycon
>>> displacement_X_arr = displacement_X.as_array()
>>> print(displacement_X_arr)

```

### Export to a pandas Dataframe
```pycon
>>> displacement_X_df = displacement_X.as_data_frame()
>>> print(displacement_X_df)

```
