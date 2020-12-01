import os

from ansys.dpf.post.common import Grouping as grouping
from ansys.dpf.post.common import _AvailableKeywords as available_keywords
from ansys.dpf.post.post_utility import result
from ansys import dpf

"""Post-processing module. Using Data Processing Framework.
Allow to create a result object, then use it to get wanted results.

Example
-----
from ansys.dpf import post
result = post.result("file.rst")
disp = result.nodal_displacement() 

"""

dpf.core.start_local_server()

def build_docs(path=None):
    """Build HTML documentation.  This outputs all available
    operator types for the loaded operators.

    HTML is saved in the directory mentionned as parameters.
    
    Example
    -----
    from ansys.dpf import post
    post.build_doc("d:/temp/dpf_doc.html")
    
    Parameters
    -----
    str: output path for the documentation. Default is current directory.
    """
    if path is None:
        path = os.getcwd()
        path += "/dataProcessingDoc.html"
    doc_op = dpf.core.Operator('html_doc')
    doc_op.inputs.output_path.connect(path)
    doc_op.run()
    



