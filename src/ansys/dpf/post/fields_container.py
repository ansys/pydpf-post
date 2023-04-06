from __future__ import annotations

from ansys import dpf
import ansys.dpf.core as dpf

from ansys.dpf.core.property_field import PropertyField

from typing import Union, Dict

from collections.abc import Sequence
import copy

class LabelSpaceKV:
    def __init__(self, _dict: Dict[str, int], _field):
        self._dict  = _dict
        self._field = _field

    def dict(self):
        return self._dict

    def field(self):
        return self._field
    
    def __str__(self):
        field_str = str(self._field).replace('\n', ' ')
        return f"Label Space: {self._dict} with field {field_str}"

class PropertyFieldsContainer(Sequence):
    def __init__(self, fields_container=None, server=None):
        # default constructor
        self._labels             = []    # used by Dataframe
        self.scopings           = []
        self.server             = None  # used by Dataframe
        
        self.label_spaces       = []
        self.ids                = []
        
        # fields_container copy
        if fields_container is not None:
            self._labels         = copy.deepcopy(fields_container.labels)
            self.scopings       = copy.deepcopy(fields_container.scopings)
            self.server         = copy.deepcopy(fields_container.server)
            
            self.label_spaces   = copy.deepcopy(fields_container.label_spaces)
            self.ids            = copy.deepcopy(fields_container.ids)
        # server copy
        if server is not None:
            self.server = server
        
    # Collection
    def __str__(self):
        str = f"DPF PropertyFieldsContainer with {len(self)} fields\n"
        for idx, ls in enumerate(self.label_spaces):
            str += f"\t {idx}: {ls}\n"
        
        return str

    @property
    def labels(self):
        return self._labels
    
    @labels.setter
    def labels(self, vals):
        self.set_labels(vals)

    def set_labels(self, labels):
        if len(self._labels) != 0:
            raise ValueError("labels already set")

        for l in labels:
            self.add_label(l)
    
    def add_label(self, label):
        if label not in self._labels:
            self._labels.append(label)
            self.scopings.append([])
    
    def has_label(self, label):
        return label in self.labels
    
    # used by Dataframe
    def get_label_space(self, idx):
        return self.label_spaces[idx]._dict
    
    # used by Dataframe
    def get_label_scoping(self, label="time"):
        if label in self.labels:
            scoping_ids = self.scopings[self.labels.index(label)]
            return dpf.Scoping(ids=scoping_ids, location="")
        raise KeyError("label {label} not found")
    
    def add_entry(self, label_space: Dict[str, int], value):
        new_id = self._new_id()
        
        if hasattr(value, "_server"):
            self.server = value._server
        
        # add Label Space
        self.label_spaces.append(LabelSpaceKV(label_space, value))

        # Update IDs
        self.ids.append(new_id)

        # Update Scopings
        for label in label_space.keys():
            label_idx = self.labels.index(label)
            self.scopings[label_idx].append(new_id)
    
    def get_entries(self, label_space_or_index):
        if isinstance(label_space_or_index, int):
            idx: int = label_space_or_index
            return [self.label_spaces[idx].field()]
        else:
            _dict: Dict[str, int] = label_space_or_index
            are_keys_in_labels = [key in self.labels for key in _dict.keys()]
            if all(are_keys_in_labels):
                remaining = set(range(len(self.label_spaces)))
                for key in _dict.keys():
                    val = _dict[key]
                    to_remove = set()
                    for idx in remaining:
                        ls = self.label_spaces[idx]
                        if ls.dict()[key] != val:
                            to_remove.add(idx)
                    remaining = remaining.difference(to_remove)
                
                idx_to_field = lambda idx: self.label_spaces[idx].field()
                return list(map(idx_to_field, remaining))
            else:
                bad_idx = are_keys_in_labels.index(False)
                bad_key = _dict.keys()[bad_idx]
                raise KeyError(f"Key {bad_key} is not in labels: {self.labels}")

    def get_entry(self, label_space_or_index):
        ret = self.get_entries(label_space_or_index)

        if len(ret) != 0:
            return ret[0]

        raise IndexError("Could not find corresponding entry")

    def _new_id(self):
        if len(self.ids) == 0:
            self.last_id = 1
            return self.last_id
        else:
            self.last_id += 1
            return self.last_id

    # FieldsContainer
    def create_subtype(self, obj_by_copy):
        return PropertyField(property_field=obj_by_copy, server=self.server)

    def get_fields_by_time_complex_ids(self, timeid=None, complexid=None):
        label_space = {"time": timeid, "complex":complexid}
        return self.get_fields(label_space)

    def get_field_by_time_complex_ids(self, timeid=None, complexid=None):
        label_space = {"time":timeid, "complex": complexid}
        return self.get_field(label_space)

    def __time_complex_label_space__(self, timeid=None, complexid=None):
        raise NotImplementedError

    # used by Dataframe
    def get_fields(self, label_space):
        return self.get_entries(label_space)

    def get_field(self, label_space_or_index):
        return self.get_entry(label_space_or_index)

    def get_field_by_time_id(self, timeid=None):
        label_space = {"time":timeid}
        if self.has_label("complex"):
            label_space["complex"] = 1
            return self.get_field(label_space)

    def get_imaginary_fields(self, timeid=None):
        label_space = {"time":timeid, "complex":1}
        return self.get_fields(label_space)

    def get_imaginary_field(self, timeid=None):
        label_space = {"time":timeid, "complex":1}
        return self.get_field(label_space)

    # used by Dataframe
    def __getitem__(self, key):
        return self.get_field(key)

    def __len__(self):
        return len(self.label_spaces)

    def add_field(self, label_space, field):
        self.add_entry(label_space, field)

    def add_field_by_time_id(self, field, timeid=1):
        if not self.has_label("time"):
            self.add_label("time")

        label_space = {"time":timeid}

        if self.has_label("complex"):
            label_space["complex"] = 0
        
        self.add_field(label_space, field)

    def add_imaginary_field(self, field, timeid=1):
        if not self.has_label("time"):
            self.add_label("time")
        if not self.has_label("complex"):
            self.add_label("complex")
        
        label_space = {"time":timeid, "complex":1}
        self.add_field(label_space, field)

    def select_component(self, index):
        raise NotImplementedError

    @property
    def time_freq_support(self):
        raise NotImplementedError

    @time_freq_support.setter
    def time_freq_support(self, value):
        raise NotImplementedError

    def deep_copy(self, server=None):
        raise NotImplementedError

    def get_time_scoping(self):
        return self.get_label_scoping("time")

    def animate(self, save_as=None, deform_by=None, scale_factor=1.0, **kwargs):
        raise NotImplementedError

    def __sub__(self, fields_b):
        raise NotImplementedError

    def __pow__(self, value):
        raise NotImplementedError

    def __mul__(self, value):
        raise NotImplementedError
