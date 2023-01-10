# Project: DPF Refactoring

## Pseudocode for StaticMechanicalSimulation

Load the result files within an instance of `StaticMechanicalSimulation`

The `ansys.dpf.post.load_simulation` function detects by default the type physics
and the type of analysis based on the metadata in the underlying `ansys.dpf.core.Model` instance.
If no physics type or analysis type is recognized automatically, 
it will default to static and mechanical.
You can explicit the physics type and analysis type. 

```py
static_simulation = dpf.load_simulation(examples.simple_bar)
static_simulation = dpf.load_simulation(examples.simple_bar,
                                        physics_type=physics_types.mechanical,
                                        analysis_type=analysis_types.static)
```

To that end, we might want to provide Public Enums for available physics_types and analysis_types
instead of having them private.
```py
from ansys.dpf.post import physics_types, analysis_types
# instead of
# from ansys.dpf.post.common import _AnalysisType, _PhysicsType
```

