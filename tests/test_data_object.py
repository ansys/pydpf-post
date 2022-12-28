from ansys.dpf.post import examples, load_simulation


def test_get_displacements():
    simulation = load_simulation(examples.download_transient_result())
    disp = simulation.displacement()
    assert len(disp) == 35
    assert disp[0].ndim == 3
    assert disp[0].elementary_data_count == 3820
    disp_x = simulation.displacement(component=0)
    assert len(disp_x) == 35
    assert disp_x[0].ndim == 1
    assert disp[0].elementary_data_count == 3820
    disp_time = simulation.displacement(steps=[1, 25])
    assert len(disp_time) == 2
    assert disp_time[0].elementary_data_count == 3820


def test_get_stresses():
    simulation = load_simulation(examples.download_transient_result())

    # Nodal stress
    stress = simulation.nodal_stress()
    assert len(stress) == 35
    assert stress[0].ndim == 6
    assert stress[0].elementary_data_count == 3812
    stress_xy_nodes = simulation.nodal_stress(
        component="XY", nodes=[100, 101, 102, 103]
    )
    assert len(stress_xy_nodes) == 35
    assert stress_xy_nodes[0].ndim == 1
    assert stress_xy_nodes[0].elementary_data_count == 4

    # Elemental stress
    stress = simulation.elemental_stress(steps=[2, 15])
    assert len(stress) == 2
    assert stress[0].ndim == 6
    assert stress[0].elementary_data_count == 715

    # Raw stress
    stress = simulation.raw_stress()
    assert len(stress) == 35
    assert stress[0].ndim == 6
    assert (
        stress[0].elementary_data_count == 5720
    )  # why the number of entities is now different from data_count?
