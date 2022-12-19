# Project: DPF Refactoring

## Aim

Simple, and intuitive, yet powerful, easy-to-use entry level library which remains as an API to ``DPF-Core`` library.

## Requirements

* **Simple**: Focused on exposing data.
* **Intuitive**: Similar to other Python libraries (numpy, pandas, pyvista)
* **Powerful**: Using ``DPF-Core`` under the hood. Everything is just exposing ``DPF-Core`` operators/functions. Using Numpy and Pandas helpers.
* **Easy-to-use**: No more than three lines to retrieve/plot data.


## Current DPF-Core implementation

```py
from ansys.dpf import core as dpf
from ansys.dpf.core import examples

# Load model
model = dpf.Model(examples.static_rst)

# Define scoping for nodes with ID 4 and 6
nodes_scoping = dpf.Scoping(ids=[4, 6], location=dpf.locations.nodal)

# Evaluate displacement for the defined scoping
disp_op = model.results.displacements()
disp_op.inputs.mesh_scoping(nodes_scoping)
disp = disp_op.outputs.fields_container()

# Plot
model.metadata.meshed_region.plot(disp)
```

## Current DPF-Post implementation

```py
from ansys.dpf import post
from ansys.dpf.post import examples

solution = post.load_solution(examples.static_rst)

disp = solution.displacement(node_scoping=[4, 6])
disp.vector.plot_contour()

```


## Proposed DPF-Post implementation
Focused on intuitive manipulation of the data.
```py
from ansys.dpf.post import load_solution
from ansys.dpf.post import examples

solution = load_solution(examples.static_rst)
solution.displacement(nodes=[4,6]).plot()
```

In the future, we are aiming for a more descriptive, and imperative
API such as:

```python
solution.displacement(
          named_selection=["mycm", "myothercm"],
          time_steps = [1, 2],
          component = "X",
        ).plot(
          cmap="viridis",
          show_edges=True,
        )
```

Leveraging ``DPF-Core`` operators.


## Future implementation

### **Selector object**
For spatial, time/frequency and qualifiers (Fluids).

Two ways to apply it.

Creating an object and use it in our results call.

```py
from ansys.dpf.post import select

## Function based
selection = select(
          named_selection=["mycm", "myothercm"],
          time_steps = [1, 2],
)  # Does not change 'simulation' state

solution.displacement(selection=selection)
```

Changing state of the ``Solution`` class.

```py
simulation.activate_selection(
          time_steps = [1, 2],
) # Remains applied in 'simulation' until cleared/deactivated

simulation.displacement().plot()  # plotting only time_steps 1 and 2

simulation.deactivate_selection()
```

You can even reselect from selections (intersections).

```py
selection = select(
          time_steps = [1, 2],
)
selection = selection.apply(nodes=[1, 20])

simulation.displacement(
      selection=selection
    ).plot()
```

### **Mapping results to geometry**

```py
from ansys.dpf.post import tools
from ansys.dpf.post import select

center = (0,0,0)
normal = (1,0,0)

myselection = select(tools.plane(center,normal))
solution.displacement(
    selection=myselection,
    component="X"
    ).plot(cmap="viridis")


```
