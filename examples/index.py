from __future__ import annotations

from numpy.typing import NDArray

from typing import Union, Any, List
import numpy as np


class Index:
    def __init__(self, data: Union[Index, NDArray], name: str = "", mapping_to_idx = None):
        # copy constructor
        if isinstance(data, Index):
            other = data
            self._data      = other._data
            self._name      = other._name
            self._mapping   = other._mapping
        else:
            self._data = data
            self._name = name

            # TODO: not optimal
            if mapping_to_idx == None:
                mapping_to_idx = dict(zip(self._data, range(len(self._data))))

            self._mapping = mapping_to_idx

    def __len__(self):
        return len(self.array)

    @property
    def array(self):
        return self._data
    
    @property
    def name(self) -> str:
        return self._name
    
    def get_loc(self, key):
        return self._mapping[key]
    
    def union(self, index):
        new_data = np.append(self._data, index._data)
        new_name = self._name
        return Index(new_data, new_name)

    def intersection(self, index):
        new_data = np.intersect1d(self._data, index._data)
        new_name = self._name
        return Index(new_data, new_name)
    
    def take(self, indices):
        new_data = np.take(self._data, indices)
        new_name = self._name
        return Index(new_data, new_name)
    
    def __repr__(self):
        return f"Index({self.array}, name=\"{self.name}\")"
    

class RangeIndex(Index):
    def __init__(self, start: Union[Index, RangeIndex, int, None] = None, stop=None, step=None, name: str = ""):
        if isinstance(start, RangeIndex): # copy
            other = start
            self.__init__(other._start, other._stop, other._step, other._name)
        elif isinstance(start, Index): # copy
            raise NotImplementedError
        else:
            if start == None:
                start = 0
            if stop == None:
                stop = 0
            if step == None:
                step = 1
            super().__init__(np.arange(start, stop, step), name)
            self._start = start
            self._stop = stop
            self._step = step
    
    def union(self, index: RangeIndex):
        if isinstance(index, RangeIndex) and index._step == self._step:
            step = self._step
            name = self._name
            # self is contained in index
            if self._start >= index._start and self._stop <= index._stop:
                start = index._start
                stop = index._stop
                return RangeIndex(start, stop, step, name)
            # self contains index
            if self._start <= index._start and self._stop >= index._stop:
                start = self._start
                stop = self._stop
                return RangeIndex(start, stop, step, name)
            # self overlaps index
            if self._stop >= index._start:
                start = self._start
                stop = index._stop
                return RangeIndex(start, stop, step, name)
            # self gets overlapped by index
            if self._start <= index._stop:
                start = index._start
                stop = self._stop
                return RangeIndex(start, stop, step, name)
        
        return super().union(index)
    
    def intersection(self, index: RangeIndex):
        if isinstance(index, RangeIndex) and index._step == self._step:
            step = self._step
            name = self._name
            # self is contained in index
            if self._start >= index._start and self._stop <= index._stop:
                start = self._start
                stop = self._stop
                return RangeIndex(start, stop, step, name)
            # self contains index
            if self._start <= index._start and self._stop >= index._stop:
                start = index._start
                stop = index._stop
                return RangeIndex(start, stop, step, name)
            # self overlaps index
            if self._stop >= index._start:
                start = index._start
                stop = self._stop
                return RangeIndex(start, stop, step, name)
            # self gets overlapped by index
            if self._start <= index._stop:
                start = index._stop
                stop = self._start
                return RangeIndex(start, stop, step, name)
        
        return super().intersection(index)

    def __repr__(self):
        return f"RangeIndex(start={self._start}, stop={self._stop}, step={self._step}, name={self._name})"

class CategoricalIndex(Index):
    def __init__(self, data: Union[Index, CategoricalIndex, List[Any], NDArray[Any]], categories=[], name: str = ""):
        if isinstance(data, CategoricalIndex): # copy
            other = data
            self.__init__(other._data, other._categories, other._name)
        elif isinstance(data, Index): # copy
            other = data
            self.__init__(data=other._data, name=other._name)
        else: # base constructor
            data = np.array(data)
            if categories == None or categories==[]:
                categories = list(set(data))

            self._categories = categories

            mapping = [np.where(data == cat) for cat in self._categories]

            super().__init__(data, name, mapping)

    @property
    def categories(self):
        return self._categories


    def __repr__(self):
        return f"CategoricalIndex(data={self.array}, categories={self._categories}, name={self._name})"
    
class MultiIndex(Index):
    def __init__(self, levels: Union[MultiIndex, Index, List[NDArray[Any]]] = [], names=[]):
        if isinstance(levels, MultiIndex):
            other = levels
            self.__init__(other._levels, other._names)
        elif isinstance(levels, Index):
            other = levels
            _levels = [other._data]
            _names = [other._name]
            self.__init__(_levels, _names)
        else:
            self._levels = levels
            if len(names) != len(levels):
                missing = len(levels)-len(names)
                names += [""] * missing
            self._names = names

            # workaround to assign a list of tuples to np array
            unpacked = list(zip(*self._levels))
            np_arr = np.empty(len(unpacked), dtype=object)
            np_arr[:] = unpacked
            super().__init__(data=np_arr)

    @staticmethod
    def from_arrays(levels=[], names=[]):
        return MultiIndex(levels, names)
    
    @staticmethod
    def from_product(array_list, names=[]):
        arr_l = np.meshgrid(*array_list)
        arrays = [arr.ravel() for arr in arr_l]
        return MultiIndex(arrays, names)

    def to_flat_index(self):
        return Index(data=self._data)

    def __repr__(self):
        return f"MultiIndex({self.array}, names={self._names})"