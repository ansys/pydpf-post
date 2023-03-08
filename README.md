# PyDPF-Post - Ansys Data Post-Processing Framework
[![PyAnsys](https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC)](https://docs.pyansys.com/)
[![Python](https://img.shields.io/pypi/pyversions/ansys-dpf-post?logo=pypi)](https://pypi.org/project/ansys-dpf-post/)
[![pypi](https://badge.fury.io/py/ansys-dpf-post.svg?logo=python&logoColor=white)](https://pypi.org/project/ansys-dpf-post)
[![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The Data Processing Framework (DPF) is designed to provide numerical
simulation users/engineers with a toolbox for accessing and
transforming simulation data.

The Python `ansys-dpf-post` package provides a high level, physics oriented API for postprocessing.
Loading a simulation (defined by its result files) allows you to extract simulation metadata as well
as results and apply postprocessing operations on it.

This module leverages the PyDPF-Core project's ``ansys-dpf-core`` package and can
be found by visiting [PyDPF-Core
GitHub](https://github.com/pyansys/pydpf-core).  Use ``ansys-dpf-core`` for
building more advanced and customized workflows using Ansys DPF.

## Documentation

Visit the [PyDPF-Post Documentation](https://postdocs.pyansys.com) for a
detailed description of the package, or see the [Examples
Gallery](https://postdocs.pyansys.com/examples/index.html) for more
detailed examples.

## Installation

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

## Brief Demo

Provided you have ANSYS 2023 R1 installed, a DPF server starts
automatically once you start using PyDPF-Post.
Loading a simulation to extract and post-process results:

```pycon
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> simulation = post.load_simulation(examples.download_crankshaft())
>>> displacement = simulation.displacement()
>>> print(displacement)
```
```pycon
             results         U
              set_id         3
      node      comp          
      4872         X -3.41e-05
                   Y  1.54e-03
                   Z -2.64e-06
      9005         X -5.56e-05
                   Y  1.44e-03
                   Z  5.31e-06
       ...
```
```pycon
>>> displacement.plot()
```
![Example Displacement plot Crankshaft](https://github.com/pyansys/dpf-post/raw/master/docs/source/images/crankshaft_disp.png)
```pycon
>>> stress_eqv = simulation.stress_eqv_von_mises_nodal()
>>> stress_eqv.plot()
```
![Example Stress plot Crankshaft](https://github.com/pyansys/dpf-post/raw/master/docs/source/images/crankshaft_stress.png)

To run PyDPF-Post with Ansys versions starting from 2021 R1 to 2022 R2, use the following legacy PyDPF-Post 
tools:

```pycon
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> solution = post.load_solution(examples.download_crankshaft())
>>> stress = solution.stress()
>>> stress.eqv.plot_contour(show_edges=False)
```
![Example Stress plot Crankshaft](https://github.com/pyansys/dpf-post/raw/master/docs/source/images/crankshaft_stress.png)

## License

``PyDPF-Post`` is licensed under the MIT license.  For more information, see the
[LICENSE](https://github.com/pyansys/dpf-post/raw/master/LICENSE).
