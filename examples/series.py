
from __future__ import annotations

from typing import Tuple

from collections.abc import Sequence, Iterator

import math
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
            return self._series._data[idx]
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

#                   LEFT             RIGHT
########### +-------------------+---------------+ 
#           | column index name | columns labels|   <-- only for DataFrame
#  HEADER   +-------------------+---------------+
#           | row index name    |
########### +-------------------+---------------+
#           | row labels        | Values        |
#           |                   |               |
#           |                   |               |
#  VALUES   |                   |               |
#           |                   |               |
#           |                   |               |
#           |                   |               |
########### +-------------------+---------------+
#  FOOTER       Shape [x rows]
###########

class SeriesFormatter:
    def __init__(self, series: Series, row_end: int, precision: int = 6):
        self._series = series

        # avoid wasting processing time
        self._row_end = row_end

        # precompute lengths
        self._precision = precision
        self._left_width  = self._compute_left_width()
        self._right_width = self._compute_right_width()
        self._value_width = self._compute_value_width()
    
    def _compute_left_width(self):
        #name_size = len(self._series.name)
        label_size = len(self._series.index.name)
        max_size = label_size #max(name_size, label_size)

        if isinstance(self._series.index,post_idx.MultiIndex):
            max_size = self._compute_label_width() * len(self._series.index._levels)

        return self._next_multiple(max_size, 4)

    def _compute_right_width(self):
        value_width = self._compute_value_width()

        return self._next_multiple(value_width, 4)

    def _compute_value_width(self):
        
        value_width = int(math.log10(len(self._series))) + 1
        if np.issubdtype(self._series._data.dtype, np.floating):
            value_width = self._precision + 3 + 1# integer part, dot, minus sign, + space
        
        return self._next_multiple(value_width, 4)

    def _compute_label_width(self):
        label_size = len(self._series.index.name)
        if isinstance(self._series.index,post_idx.MultiIndex):
            label_size = max( [
                                max([len(str(e)) for e in self._series.index._data[i]])
                                for i in range(min(self._row_end,len(self._series)))
                            ]
                           )
            names_size = max(map(len,self._series.index._names))
            label_size = max(label_size, names_size)
        return self._next_multiple(label_size, 4)

    def _get_value_fmt(self):
        value_fmt = "{: " + f"{self._value_width}" + "}"
        if np.issubdtype(self._series._data.dtype, np.floating):
            value_fmt = "{: " + f"{self._precision}e" + "}"
        return value_fmt

    def _next_multiple(self, x, n):
        if (x % n) != 0:
            return (int(x/n)+1)*n
        return x

    def generate_string(self):
        ret = ""
        ret += self._generate_header()
        ret += self._generate_values()
        ret += self._generate_footer()
        return ret
    
    def _generate_header(self):
        label_fmt = "{:<" + f"{self._compute_label_width()}" + "}"
        
        ret = f"Name: {self._series._name}\n"

        if isinstance(self._series._index,post_idx.MultiIndex):
            for n in self._series._index._names:
                ret += label_fmt.format(n)
        else:
            ret += label_fmt.format(self._series.index.name)
        
        ret += "\n"
        return ret

    def _generate_values(self):
        ret = ""

        label_fmt = "{:<" + f"{self._compute_label_width()}" + "}"
        value_fmt = self._get_value_fmt()
        
        for zidx, val in enumerate(self._series[:self._row_end]):
            id  = self._series.index._data[zidx]
            if not isinstance(id, tuple):
                id = (id,)
            for e in id:
                ret += label_fmt.format(e)
            ret += value_fmt.format(val)
            ret += "\n"
        return ret
    
    def _generate_footer(self):
        return f"[{len(self._series)} rows]"

class SeriesIterator(Iterator):
    def __init__(self, series):
        self._series = series
        self.idx = 0

    def __next__(self):
        if self.idx < len(self._series):
            ret = self._series.iloc[self.idx]
            self.idx += 1
            return ret
        raise StopIteration
    
    def __iter__(self):
        return SeriesIterator(self._series)

class Series(Sequence):
    def __init__(self, data, index=None, name=""):
        if index == None:
            index = post_idx.RangeIndex(0,len(data), 1)
        
        self._data = data
        self._index = index
        self._name = name
    
    def get(self, k):
        idx = self._index._mapping[k]
        return self._data[idx]

    def __iter__(self):
        return SeriesIterator(self)

    def __getitem__(self, k):
        return self.iloc[k]

    def __len__(self):
        return len(self._index)

    @property
    def index(self) -> post_idx.Index:
        return self._index
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def iloc(self):
        return SeriesILoc(self)
    
    def __repr__(self):

        return SeriesFormatter(self, 10, 10).generate_string()


        ret_str = ""
        if self.name:
            ret_str += f"Name: {self.name}\n"
        
        for idx in self._index._data:
            val = self.get(idx)
            ret_str += f"{idx}\t{val}\n"
        return ret_str