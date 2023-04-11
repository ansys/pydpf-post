
import pandas as pd
import numpy as np
import index
import series

import random


def random_int_list(size, _min=0, _max=100):
    return np.array([random.randrange(_min,_max) for i in range(size)])

def random_char_list(size):
    return np.array(map(lambda i: chr(i),random_int_list(size, ord('A'), ord('Z'))))

def test_range_index(range_index_type):
    range_index1 = range_index_type(0, 8, name="test_range1")
    range_index2 = range_index_type(4, 10)

    union_index = range_index1.union(range_index2)
    intersection_index = range_index1.intersection(range_index2)

    extracted_index = range_index2.take([0, 1, 2])
    print(extracted_index)
    print(union_index)
    print(intersection_index)

def test_cat_index(cat_index_type):
    cat_index1 = cat_index_type(data=['a','a','b','c','a'])
    cat_index2 = cat_index_type(data=['f','f','e','g','f'])

    print(cat_index1)
    print(cat_index2)
    union_index = cat_index1.union(cat_index2)
    intersection_index = cat_index1.intersection(cat_index2)
    print(union_index)
    print(intersection_index)


def test_multi_index(multi_index_type):
    multi_index1 = multi_index_type.from_arrays([["a","b","c"], [1,2,3]])
    multi_index2 = multi_index_type.from_product([[4,5,6], ["a","b"], ["oui", "non"]])

    print(multi_index1)
    print(multi_index2)

def test_series_iloc(series):
    
    elem_0      = series.iloc[0]
    elems_list  = series.iloc[[0,3,4]]
    elems_slice = series.iloc[0:5]
    elems_bool  = series.iloc[[True, False, False, True, True]]

    print(series)

    print(elem_0)
    print(elems_list)
    print(elems_slice)
    print(elems_bool)
    print(elems_bool.index)

def test_series_multi_idx(series):
    print(series)
    print(series["one"])
    print(series["one"]["foo"])
    print(series["one"]["foo"]["X"])

    print(series.loc[("one","foo","X")])
    print(series.xs(("foo"), level=("name")))
    print(series.xs(("bar","Z"), level=("name","component")))

print("===Range Indexes===")
print("----pandas----")
test_range_index(pd.RangeIndex)
print("----post_index----")
test_range_index(index.RangeIndex)

print("===Categorical Indexes===")
print("----pandas----")
test_cat_index(pd.CategoricalIndex)
print("----post_index----")
test_cat_index(index.CategoricalIndex)

print("====MultiIndex====")
print("----pandas----")
test_multi_index(pd.MultiIndex)
print("----post_index----")
test_multi_index(index.MultiIndex)


sample_intl = random_int_list(5)

pd_series = pd.Series(data=sample_intl,name="test_pd_series", index=pd.RangeIndex(0,5,1, name="X"))
post_series = series.Series(data=sample_intl, name="test_post_series")
print("====Series====")
print("----pandas----")
test_series_iloc(pd_series)
print("----post_series----")
test_series_iloc(post_series)


rand_data = np.random.randn(3*2*3)

multi_pd_index = pd.MultiIndex.from_product([["one","two","three"],["foo","bar"], ["X","Y","Z"]], names=["num","name","component"])
multi_pd_series = pd.Series(data=rand_data, name="multi_pd_series",index=multi_pd_index)
test_series_multi_idx(multi_pd_series)