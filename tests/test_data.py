from ansys.dpf.post import examples, load_simulation


def test_data():
    simulation = load_simulation(examples.download_transient_result())
    stressObject = simulation.nodal_stress()  # THIS DOES NOT WORK WITH raw_stress()
    stress = stressObject[0]

    assert stress.name == "stress_0.s"
    assert stress.location == "Nodal"
    assert len(stress.data) == stress.n_data
    assert len(stress.ids) == stress.n_data
    assert len(stress.data[0]) == stress.n_dim
    # assert type(stress.time_freq)
    assert stress.unit == "Pa"
    stress.plot()
