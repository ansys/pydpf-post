
from __future__ import annotations

import numpy as np
import index as post_idx

class SeriesILoc:
    def __init__(self, series: Series):
        self._series = series
    
    def __getitem__(self, key):
        indices = []

        if isinstance(key, int):
            idx = key
            indices = [idx]
        elif isinstance(key, list):
            l = key
            if len(l) != 0:
                first_elem = l[0]
                if isinstance(first_elem, bool):
                    list_bool = l
                    indices = self._get_indices_from_bool_list(list_bool)
                elif isinstance(first_elem, int):
                    list_int = l
                    indices = list_int
        elif isinstance(key, slice):
            sl = key
            start, stop, step = sl.indices(len(self._series._data))
            indices = np.arange(start, stop, step)

        # get indices to keep
        new_data = self._series._data.take(indices)
        new_index = self._series._index.take(indices)
        new_name = self._series.name
        return Series(new_data, new_index, new_name)


    def _get_indices_from_bool_list(self, list_bool):
        return np.arange(len(list_bool))[list_bool]
    

class Series:
    def __init__(self, data, index=None, name=None):
        if index == None:
            index = post_idx.RangeIndex(0,len(data), 1, name)
        
        self._data = data
        self._index = index
    
    def get(self, k):
        idx = self._index._mapping[k]
        return self._data[idx]

    def __getitem__(self, k):
        return self.iloc[k]

    def __len__(self):
        return len(self._data)

    @property
    def index(self) -> post_idx.Index:
        return self._index
    
    @property
    def name(self) -> str:
        return self._index.name
    
    @property
    def iloc(self):
        return SeriesILoc(self)
    
    def __repr__(self):
        ret_str = ""
        if self.name:
            ret_str += f"Name: {self.name}\n"
        
        for idx in self._index._data:
            val = self.get(idx)
            ret_str += f"{idx}\t{val}\n"
        return ret_str