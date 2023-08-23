# PyDPF-Post - Ansys Data PostProcessing Framework
[![PyAnsys](https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC)](https://docs.pyansys.com/)
[![Python](https://img.shields.io/pypi/pyversions/ansys-dpf-post?logo=pypi)](https://pypi.org/project/ansys-dpf-post/)
[![pypi](https://badge.fury.io/py/ansys-dpf-post.svg?logo=python&logoColor=white)](https://pypi.org/project/ansys-dpf-post)
[![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Ansys Data Processing Framework (DPF) provides numerical simulation 
users and engineers with a toolbox for accessing and transforming simulation 
data. With DPF, you can perform complex preprocessing or postprocessing of
large amounts of simulation data within a simulation workflow.

The Python `ansys-dpf-post` package provides a high-level, physics-oriented
API for postprocessing. Loading a simulation (defined by its results files)
allows you to extract simulation metadata and results and then apply
postprocessing operations on them.

The latest version of DPF supports Ansys solver results files for:

  - Mechanical APDL (`.rst`, `.mode`, `.rfrq`, `.rdsp`)
  - LS-DYNA (`.d3plot`, `.binout`)
  - Fluent (`.cas/dat.h5`, `.flprj`)
  - CFX (`.cad/dat.cff`, `.flprj`)

For more information on file support, see the [main page](https://dpf.docs.pyansys.com/version/stable/index.html)
in the PyDPF-Core documentation.

PyDPF-Post leverages the PyDPF-Core project's ``ansys-dpf-core`` package, which is
available at [PyDPF-Core GitHub](https://github.com/ansys/pydpf-core).
Use the ``ansys-dpf-core`` package for building more advanced and customized workflows
using Ansys DPF.

## Documentation and issues

Documentation for the latest stable release of PyPDF-Post is hosted at
[PyDPF-Post documentation](https://post.docs.pyansys.com/version/stable/).

In the upper right corner of the documentation's title bar, there is an option for switching from
viewing the documentation for the latest stable release to viewing the documentation for the
development version or previously released versions.

You can also [view](https://cheatsheets.docs.pyansys.com/pydpf-post_cheat_sheet.png) or
[download](https://cheatsheets.docs.pyansys.com/pydpf-post_cheat_sheet.pdf) the
PyDPF-Post cheat sheet. This one-page reference provides syntax rules and commands
for using PyDPF-Post.

On the [PyDPF-Post Issues](https://github.com/ansys/pydpf-post/issues) page,
you can create issues to report bugs and request new features. On the
[PyDPF-Post Discussions](https://github.com/ansys/pydpf-post/discussions) page or the [Discussions](https://discuss.ansys.com/)
page on the Ansys Developer portal, you can post questions, share ideas, and get community feedback. 

To reach the project support team, email [pyansys.core@ansys.com](mailto:pyansys.core@ansys.com).

## Installation

To install this package, run this command:

```
pip install ansys-dpf-post
```

You can also clone and install this package with these commands:

```
git clone https://github.com/ansys/pydpf-post
cd pydpf-post
pip install . --user
```

## Brief demo

Provided you have Ansys 2023 R1 or later installed, a DPF server automatically starts
once you start using PyDPF-Post.

To load a simulation for a MAPDL result file to extract and
postprocess results, use this code:

```pycon
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> simulation = post.load_simulation(examples.download_crankshaft())
>>> displacement = simulation.displacement()
>>> print(displacement)
```
```pycon
             results       U (m)
             set_ids           3
 node_ids components            
     4872          X -3.4137e-05
                   Y  1.5417e-03
                   Z -2.6398e-06
     9005          X -5.5625e-05
                   Y  1.4448e-03
                   Z  5.3134e-06
      ...        ...         ...
```
```pycon
>>> displacement.plot()
```
![Example Displacement plot Crankshaft](https://github.com/ansys/pydpf-post/raw/master/docs/source/images/crankshaft_disp.png)
```pycon
>>> stress_eqv = simulation.stress_eqv_von_mises_nodal()
>>> stress_eqv.plot()
```
![Example Stress plot Crankshaft](https://github.com/ansys/pydpf-post/raw/master/docs/source/images/crankshaft_stress.png)

To run PyDPF-Post with Ansys 2021 R1 through 2022 R2, use this code to
start the legacy PyDPF-Post tools:

```pycon
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> solution = post.load_solution(examples.download_crankshaft())
>>> stress = solution.stress()
>>> stress.eqv.plot_contour(show_edges=False)
```
![Example Stress plot Crankshaft](https://github.com/ansys/pydpf-post/raw/master/docs/source/images/crankshaft_stress.png)

## License and acknowledgements

PyDPF-Post is licensed under the MIT license. For more information, see the
[LICENSE](https://github.com/ansys/pydpf-post/raw/master/LICENSE) file.

PyDPF-Post makes no commercial claim over Ansys whatsoever. This library
extends the functionality of Ansys DPF by adding a Python interface
to DPF without changing the core behavior or license of the original
software.