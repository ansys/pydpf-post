import ansys.dpf.core as core

import ansys.dpf.post as dpf


def test_dataframe_from_fields_container(simple_bar):
    model = core.Model(simple_bar)
    fc = model.results.displacement().eval()
    df = dpf.DataFrame(data=fc)
    assert df


def test_dataframe_len(multishells, transient_rst):
    # Case of one set with two nodal fields
    model = core.Model(multishells)
    fc = model.results.stress.on_all_time_freqs.on_location("Nodal").eval()
    df = dpf.DataFrame(data=fc)
    assert len(df) == 7079
    # Case of several sets with one field per set
    model = core.Model(transient_rst)
    fc = model.results.displacement.on_all_time_freqs.eval()
    df = dpf.DataFrame(data=fc)
    assert len(df) == 133700


def test_dataframe_str(transient_rst):
    model = core.Model(transient_rst)
    # print(model)
    fc = model.results.displacement.on_all_time_freqs().eval()
    df = dpf.DataFrame(data=fc)
    txt = str(df)
    print()
    print(txt)

    print("\nReference:\n")

    field = fc.get_field_by_time_id(timeid=1)
    print(field.scoping.ids)
    print(field.get_entity_data_by_id(id=525))
    print(field.get_entity_data_by_id(id=534))
    field = fc.get_field_by_time_id(timeid=35)
    print(field.get_entity_data_by_id(id=3817))
    print(field.get_entity_data_by_id(id=3825))


def test_dataframe_head_tail(transient_rst):
    model = core.Model(transient_rst)
    fc = model.results.displacement.on_all_time_freqs().eval()
    df = dpf.DataFrame(data=fc)
    txt = df.head()
    print("\ndf.head()\n")
    print(txt)

    print("\nReference:\n")
    field = fc.get_field_by_time_id(timeid=1)
    ids = field.scoping.ids
    print(ids[:5])
    print(field.get_entity_data_by_id(id=ids[0]))
    print(field.get_entity_data_by_id(id=ids[1]))
    print(field.get_entity_data_by_id(id=ids[2]))
    print(field.get_entity_data_by_id(id=ids[3]))
    print(field.get_entity_data_by_id(id=ids[4]))

    txt = df.head(6)
    print("\ndf.head(6)\n")
    print(txt)

    print("\nReference:\n")
    print(ids[5])
    print(field.get_entity_data_by_id(id=ids[5]))

    txt = df.tail()
    print("\ndf.tail()\n")
    print(txt)

    print("\nReference:\n")
    field = fc.get_field_by_time_id(timeid=35)
    ids = field.scoping.ids
    print(ids[-5:])
    print(field.get_entity_data_by_id(id=ids[-5]))
    print(field.get_entity_data_by_id(id=ids[-4]))
    print(field.get_entity_data_by_id(id=ids[-3]))
    print(field.get_entity_data_by_id(id=ids[-2]))
    print(field.get_entity_data_by_id(id=ids[-1]))

    txt = df.tail(6)
    print("\ndf.tail(6)\n")
    print(txt)

    print("\nReference:\n")
    print(ids[-6])
    print(field.get_entity_data_by_id(id=ids[-6]))
