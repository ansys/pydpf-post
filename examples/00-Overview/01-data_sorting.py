"""
.. _data_sorting:

Data sorting
============
This example shows how to retrieve data sorted by node ID.

"""

###############################################################################
# Imports and load model
# ~~~~~~~~~~~~~~~~~~~~~~
# Import modules
import ansys.dpf.core as core

from ansys.dpf.post import examples, load_solution
from ansys.dpf.post.data_object import DataObject

simulation = load_solution(examples.static_rst, legacy=False)
# simulation = load_solution(examples.download_transient_result(), legacy=False)
# simulation = load_solution(examples.download_all_kinds_of_complexity(), legacy=False)
disp = simulation.displacement()
disp.sort()

###############################################################################
# Create artificial FieldsContainer
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: change this as Post user will not interact with FieldsContainer. Create
# a DataObject from factory? Similar to fieldscontainer_factory in Core?

# Create FieldsContainer from scratch with two different fields.
# Note that the two fields contain different scopings to account for the most
# generic case.
fc = core.FieldsContainer()
field_1 = core.field_from_array([10, 40, 30, 20])
field_1.scoping.ids = [1, 4, 3, 2]
field_2 = core.field_from_array([60, 50, 70, 90, 40])
field_2.scoping.ids = [6, 5, 7, 9, 4]
fc.labels = ["mesh_region"]
fc.add_field({fc.labels[0]: 1}, field_1)
fc.add_field({fc.labels[0]: 2}, field_2)

###############################################################################
# Create DataObject
# ~~~~~~~~~~~~~~~~~
# TODO: This should be the first section of the example, creating a DataObject from arrays/lists

# Create DataObject from FieldsContainer and check that the data is not sorted:
data = DataObject(fc)
data.is_sorted()

###############################################################################
# Sort data in the DataObject by node IDs and check the ``is_sorted`` property:
data.sort()
data.is_sorted()

###############################################################################
# TODO: How will the users interact with the dataobject? How do we expose data, ids, etc?

# Show sorted data
print(f"The first field ids are: {data._fc[0].scoping.ids}")
print(f"The first field data is: {data._fc[0].data}")
print(f"The second field ids are: {data._fc[1].scoping.ids}")
print(f"The second field data is: {data._fc[1].data}")

# TODO: add results from testing performance. How slow is sorting for a mesh of 0.5m nodes? 1m? 5m?
