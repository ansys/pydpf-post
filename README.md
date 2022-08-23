# DPF-Post - Ansys Data Post-Processing Framework
[![PyAnsys](https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC)](https://docs.pyansys.com/)
[![Python](https://img.shields.io/badge/Python-3.8-blue)](https://pypi.org/project/ansys-dpf-post/)
[![pypi](https://badge.fury.io/py/ansys-dpf-post.svg?logo=python&logoColor=white)](https://pypi.org/project/ansys-dpf-post)
[![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The Data Processing Framework (DPF) is designed to provide numerical
simulation users/engineers with a toolbox for accessing and
transforming simulation data. DPF can access data from solver result
files as well as several neutral formats (csv, hdf5, vtk,
etc.). Various operators are available allowing the manipulation and
the transformation of this data.

The Python `ansys-dpf-post` package provides a simplified Python
interface to DPF, thus enabling rapid postprocessing without ever
leaving a Python environment. 

This module leverages the DPF-Core project's ``ansys-dpf-core`` package and can
be found by visiting [PyDPF-Core
GitHub](https://github.com/pyansys/pydpf-core).  Use ``ansys-dpf-core`` for
building more advanced and customized workflows using Ansys's DPF.

Visit the [DPF-Post Documentation](https://postdocs.pyansys.com) for a
detailed description of the package, or see the [Examples
Gallery](https://postdocs.pyansys.com/examples/index.html) for more
detailed examples.


### Installation

Install this repository with:

```
pip install ansys-dpf-post
```

You can also clone and install this repository with:

```
git clone https://github.com/pyansys/pydpf-post
cd pydpf-post
pip install . --user
```

### Running DPF-Post
Provided you have ANSYS 2021R1 installed, a DPF server will start
automatically once you start using DPF-Post.  Should you wish to use
DPF-Post without 2020R1, see the [DPF Docker](https://dpfdocs.pyansys.com/getting_started/docker.html) documentation.

Opening and plotting a result file generated from Ansys workbench or
MAPDL is as easy as:

```python
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> solution = post.load_solution(examples.multishells_rst)
>>> stress = solution.stress()
>>> stress.xx.plot_contour(show_edges=False)
```

![Example Stress Plot](https://github.com/pyansys/dpf-post/raw/master/docs/source/images/main_example.png)


Or extract the raw data as a `numpy` array with:

```python
>>> stress.xx.get_data_at_field(0)
array([-3.37871094e+10, -4.42471752e+10, -4.13249463e+10, ...,
        3.66408342e+10,  1.40736914e+11,  1.38633557e+11])
```

### Key Features


**Computational Efficiency**

The DPF-Post module is based on DPF Framework that been developed with
a data framework that localizes the loading and post-processing within
the DPF server, enabling rapid post-processing workflows as this is
written in C and FORTRAN.  At the same time, the DPF-Post Python
module presents the result in Pythonic manner, allowing for the rapid
development of simple or complex post-processing scripts.


**Easy to use**

The API of DPF-Post module has been developed in order to make easy
post-processing steps easier by automating the use of DPF's chained
operators.  This allows for fast post-processing of potentially
multi-gigabyte models in a short script.  DPF-Post also details the
usage of the operators used when computing the results so you can also
build your own custom, low level scripts using the
[DPF-Core](https://github.com/pyansys/pydpf-core) module.


### License

``PyDPF-Post`` is licensed under the MIT license.  For more information, see the
[LICENSE](https://github.com/pyansys/dpf-post/raw/master/LICENSE).
