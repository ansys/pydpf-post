# Fluid use-cases

## P0: Load a fluid analysis result, print zones, list species and phases (tree)
### Server-side
- Improve print of MeshesContainer to mimic FieldsContainer, of MeshedRegion to print faces
- Add hierarchy between cell to face zones
- Verify vectorization of results


### Pseudo-code
```pycon
>>> from ansys.dpf import post
>>> simulation = post.FluidSimulation(r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj")
>>> print(simulation)
Fluid Simulation.


Data Sources
------------------------------
D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj

DPF Model
------------------------------
DPF Result Info 
  Analysis: unknown 
  Physics Type: unknown 
  Unit system: MKS: m, kg, N, s, V, A, degC 
  Available results: 
    density (Density) :Density 
    enthalpy (Enthalpy) :Enthalpy 
    mass_flow_rate (Mass Flow Rate) :Mass Flow Rate 
    mean_x_velocity (Mean X Velocity) :Mean X Velocity 
    mean_y_velocity (Mean Y Velocity) :Mean Y Velocity 
    mean_z_velocity (Mean Z Velocity) :Mean Z Velocity 
    rms_x_velocity (Rms X Velocity) :RMS X Velocity 
    rms_y_velocity (Rms Y Velocity) :RMS Y Velocity 
    rms_z_velocity (Rms Z Velocity) :RMS Z Velocity 
    static_pressure (Static Pressure) :Static Pressure 
    static_temperature (Static Temperature) :Static Temperature 
    velocity (Velocity) :SV_VELOCITY
  Available Qualifiers Labels:
    - phases: water (1), oil (2)
    - species: O2 (1), H2O (2)
------------------------------
DPF  Meshed Region: 
  16660 nodes 
  13856 elements 
  Unit:  
  With solid (3D) elements
------------------------------
DPF  Time/Freq Support: 
  Number of sets: 3 
Cumulative     Time (s)       LoadStep       Substep         
1              0.009587       23300          1               
2              0.009593       23315          1               
3              0.009600       23334          1               
>>> print(simulation.mesh)
DPF  Mesh: 
  16660 nodes 
  13856 elements 
  Unit:  
  With solid (3D) elements
>>> print(simulation.zones)
Cell zones
  11093 polyhedra cells,  zone id: 283
Face zones
  902 polygonal symmetry faces,  zone id: 217
  34 polygonal velocity-inlet faces,  zone id: 237
  58 polygonal pressure-outlet faces,  zone id: 238
  65 polygonal pressure-outlet faces,  zone id: 239
  65 polygonal pressure-outlet faces,  zone id: 240
  108 polygonal pressure-outlet faces,  zone id: 241
  1119 polygonal wall faces,  zone id: 2
>>> print(simulation.species)
[Species<name: "O2", id: 1>, Species<name: "H2O", id: 2>]
>>> print(simulation.phases)
[Phase<name: "water", id: 1>, Phase<name: "oil", id: 2>]
```

## P0: Extract result data for a given variable and render it
on a whole volume mesh / on a selection of zones for steady and unsteady analyses

### Server-side
- Possible for scalar variables
- For vector by next week -> done
- Tensor results not yet available (long-term)
- MeshInfoProvider:
  - CellToFaceZones
  - Print as a tree with names
  - List of cell zones, list of face zones, node zones with their names
- Support zone label using GenericSupport

### PyDPF-core
- Expose Zone/Part pin
 -Expose faces
- Expose Zones, API to get face zones of a cell zone
- Improve MeshesContainer API for queries on zones...

### PyDPF-post
- FluidSimulation -> done

### Pseudo-code
```pycon
>>> from ansys.dpf import post
>>> simulation = post.FluidSimulation(r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj")
>>> simulation.result_info
DPF Result Info 
  Analysis: unknown 
  Physics Type: unknown 
  Unit system: MKS: m, kg, N, s, V, A, degC 
  Available results: 
    density (Density) :Density 
    enthalpy (Enthalpy) :Enthalpy 
    mass_flow_rate (Mass Flow Rate) :Mass Flow Rate 
    mean_x_velocity (Mean X Velocity) :Mean X Velocity 
    mean_y_velocity (Mean Y Velocity) :Mean Y Velocity 
    mean_z_velocity (Mean Z Velocity) :Mean Z Velocity 
    rms_x_velocity (Rms X Velocity) :RMS X Velocity 
    rms_y_velocity (Rms Y Velocity) :RMS Y Velocity 
    rms_z_velocity (Rms Z Velocity) :RMS Z Velocity 
    static_pressure (Static Pressure) :Static Pressure 
    static_temperature (Static Temperature) :Static Temperature 
    velocity (Velocity) :SV_VELOCITY
  Available Qualifiers Labels:
    - phases: water (1), oil (2)
    - species: O2 (1), H2O (2)
>>> simulation.result_info.density
Available Qualifiers:
  {"phase": 1, "species":1}
  {"phase": 1, "species":2}
  {"phase": 1, "species":3}
  {"phase": 2, "species":1}
  {"phase": 2, "species":2}
>>> # Extract density on all zones for phase 1
>>> density = simulation.density(zone_ids=[-1], phases=[1])
>>> print(density)
>>> # Plotting the Dataframe
>>> density.plot()
>>> # Plotting the result as a volume plot would require extracting the result on the cells
>>> simulation.volume_plot(simulation.density(location="cells"))
```

## P0: Create a cut-plane (section) through zones and plot and animate pressure and velocity on it

### Pseudo-code
```pycon
>>> from ansys.dpf import post
>>> simulation = post.FluidSimulation(r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj")
>>> # Create the section support
>>> plane = Plane(center=[0., 0., 0.], normal=[1., 0., 0.])

>>> # # v1: result extraction on section selection
>>> selection = post.selection.Selection()
>>> selection.select_section(plane)
>>> section_pressure = simulation.pressure(selection=selection)
>>> section_pressure.animate()

>>> # # v2: global result extraction -> mapping on section
>>> # Extract result of interest
>>> pressure = simulation.pressure()
>>> mapped_pressure = post.map(result=pressure, onto=plane)
>>> mapped_pressure.animate()
```

## P0: Create an iso-surface of pressure and plot velocity contour on it

### Pseudo-code
```pycon
>>> from ansys.dpf import post
>>> simulation = post.FluidSimulation(r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj")
>>> # Extract results of interest
>>> pressure = simulation.pressure()
>>> velocity = simulation.velocity()
>>> iso_surface = post.iso_surface(result=pressure)
>>> print(iso_surface)
VizualizationMesh
>>> mapped_velocity = post.map(result=velocity, onto=iso_surface)
>>> mapped_velocity.plot()
```

## P0: Extract result data from a multi-phases / combustion analysis, plot and animate the flame-front as an iso-contour of species/phases
### Server-side
OK

### Pseudo-code
```pycon
```

## P0: For an unsteady analysis get iso-contour over time for a velocity/pressure field and animate it
### Server-side:
OK

### Pseudo-code
```pycon
```

## Extract and plot result on faces, result on cells, result on nodes
### Server-side:
- Allow for faces-only meshed-region

### Pseudo-code
```pycon
```

### PyDPF
- Difference between face/cell-location for cff and face/cell-location for PyDPF:
- A result on a face zone will be exposed as face-located when

### Pseudo-code
```pycon
```

## Compute and plot drag and lift across the body (external aerodynamics)
### Server-side
- update operators to compute integrals using the correct methods
### PyDPF-Post
- Create a specific API for external aerodynamics functions

### Pseudo-code
```pycon
```

## Create a cut-plane through zones and plot and animate a pressure/velocity spatial gradient on it
### Server-side
- Operators for gradient computation from a FieldsContainer

### Pseudo-code
```pycon
```

## Analysis of multi-phases for offshore structures, free-surface analysis, track interface between two phases, different methods to reconstruct the boundary -> iso-contour

### Pseudo-code
```pycon
```

## Extract cell data and average it on nodes (fluent)

### Pseudo-code
```pycon
```

## Extract node data and average it on faces (cfx)

### Pseudo-code
```pycon
```

## For a steady analysis, create a "source" for computing streamlines, plot the streamlines
### Server-side
- Operator for creation of a "particle source"
- Operator for computation of streamlines based on particle-source and FieldsContainer

### Pseudo-code
```pycon
```

## Compute vorticity and plot vortices
### Server-side
- q-criterion operator?
### PyDPF-Post
 -Simulation.vorticity or "compute_vorticity" of a velocity field

### Pseudo-code
```pycon
```
