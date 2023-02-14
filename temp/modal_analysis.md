# DPF Refactoring: Modal Mechanical Analysis

Pseudocode for ModalMechanicalSimulation

## Extract results

### Typical result signature:

```pycon
def result(
        self,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], 
            Tuple[Union[int, List[int], Union[int, List[int]], None] = None,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
    ) -> DataFrame:
```

### Use-cases

#### Extract displacements along X for nodes 1, 2 and 3 at f=0.05Hz

```pycon
displacement_X = harmonic_simulation.displacement(
    components=["X"], node_ids=[1, 2, 3], frequencies=[0.05]
)
```

#### Extract nodal XY stresses for elements 1, 2 and 3 at set 1

```pycon
stress_XY = harmonic_simulation.elemental_stress(
    components=["XY"], element_ids=[1, 2, 3], set_ids=[1]
)
```

#### Extract first principal nodal stress for a named (elemental or nodal) selection at all frequencies

```pycon
stress_S1 = harmonic_simulation.nodal_principal_stress(
    components=["1"], named_selections=named_selections[0]
)
```

#### Extract equivalent elemental strain for a selection at set 1

```pycon
strain_eqv = harmonic_simulation.elemental_eqv_strain(
    selection=selection, set_ids=[1]
)
```