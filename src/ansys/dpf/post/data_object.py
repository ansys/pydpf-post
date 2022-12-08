import pandas as pd    #Make sure that pandas is in the requirements
import numpy as np     #Make sure that numpy is in the requirements
from ansys.dpf.core.fields_container import FieldsContainer

class DataObject(FieldsContainer):
    """Exposes the data generated by a field container.

    Parameters
    ----------
    fields_container : class:`fieldc_from_dpf`
        fields_container
    server: str
        Name of the server.
    """

    def __init__(self, fields_container=None, server=None):
        super().__init__(fields_container._internal_obj, server)

    def __min__(self, *args, **kwargs):
        self.as_numpy_array().min()

    def __max__(self, *args, **kwargs):
        self.as_numpy_array().max()

    def as_data_frame(self, *args):
        """Returns the data from the field container as a pandas data_frame.

        Parameters
        ----------
        data : class:`ansys.dpf.core.fields_container.FieldsContainer`
            Fields container exposing the data.

        Returns
        ----------
        class:`pandas.core.frame.DataFrame`
            

        Examples
        --------
        >>> import pandas as pd
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.multishells_rst)
        >>> stress = solution.stress()
        >>> stress.xx.plot_contour(show_edges=False)
        # Here the example will have to be updated after refactoring.
        >>> data = stress.xx.get_data_at_field(0)
        >>> data_frame = pd.DataFrame(data)

        >>> data2 = [[0,1,2], [3,4,5], [6,7,8], [9,10,11]]
        >>> columns = ['x', 'y', 'z']
        >>> df = pd.DataFrame(data2, columns)
        >>> df = df.transpose()
        """

        columns = None
        for arg in args:
            columns.appends(arg)

        #load data into a DataFrame object:
        data_frame = pd.DataFrame(self.field(0), columns)
        transposed = data_frame.transpose()
        return transposed

    def as_numpy_array(self):
        return np.ndarray(self.get)

    def plot(self, **kwargs):
        self[-1].plot(**kwargs)