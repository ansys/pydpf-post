"""Module containing the ``TransientMechanicalSimulation`` class."""
from typing import List, Union
import warnings

from ansys.dpf import core
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import MechanicalSimulation, ResultCategory


class TransientMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical transient simulations."""

    def _build_time_freq_scoping(
        self,
        selection: Union[Selection, None],
        times: Union[float, List[float], None],
        time_step_ids: Union[int, List[int], None],
    ) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        # create from selection in priority
        if selection:
            return selection.time_freq_selection._evaluate_on(simulation=self)
        # else from time_step_ids
        if time_step_ids:
            if isinstance(time_step_ids, int):
                time_step_ids = [time_step_ids]
            return core.time_freq_scoping_factory.scoping_by_sets(
                cumulative_sets=time_step_ids, server=self._model._server
            )
        # else from times
        if times:
            if isinstance(times, float):
                times = [times]
            raise NotImplementedError
        # Otherwise, no argument was given, create a time_freq_scoping of the whole results
        return core.time_freq_scoping_factory.scoping_on_all_time_freqs(self._model)

    def _get_result(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            base_name:
                Base name for the requested result.
            location:
                Location requested.
            category:
                Type of result requested. See the :class:`ResultCategory` class.
            components:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        # Build the targeted time scoping
        time_scoping = self._build_time_freq_scoping(selection, times, time_step_ids)

        # Build the targeted mesh scoping
        mesh_scoping = self._build_mesh_scoping(
            selection,
            nodes,
            elements,
            named_selection,
            location=location,
        )

        # Build the list of requested results
        if category in [ResultCategory.scalar, ResultCategory.equivalent]:
            # A scalar or equivalent result has no components
            to_extract = None
            columns = [base_name]
        elif category in [ResultCategory.vector, ResultCategory.matrix]:
            # A matrix or vector result can have components selected
            to_extract, columns = self._build_components_from_components(
                base_name=base_name, category=category, components=components
            )
        elif category == ResultCategory.principal:
            # A principal type of result can have components selected
            to_extract, columns = self._build_components_from_principal(
                base_name=base_name, components=components
            )
        else:
            raise ValueError(f"'{category}' is not a valid category value.")

        # Initialize a workflow
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False

        # Instantiate the main result operator
        result_op = self._build_result_operator(
            name=base_name,
            time_scoping=time_scoping,
            mesh_scoping=mesh_scoping,
            location=location,
        )
        # Its output is selected as future workflow output for now
        out = result_op.outputs.fields_container

        # Add a step to compute principal invariants if result is principal
        if category == ResultCategory.principal:
            # Instantiate the required operator
            principal_op = self._model.operator(name="eig_values_fc")
            principal_op.connect(0, out)
            wf.add_operator(operator=principal_op)
            # Set as future output of the workflow
            out = principal_op.outputs.fields_container

        # Add a step to compute equivalent if result is equivalent
        elif category == ResultCategory.equivalent:
            # If a stress result, use one operator
            if base_name[0] == "S":
                equivalent_op = self._model.operator(name="segalmaneqv_fc")
            # If a strain result, use another
            elif base_name[0] == "E":
                equivalent_op = self._model.operator(name="eqv_fc")
            # Throw otherwise
            else:
                raise ValueError(
                    f"Category {ResultCategory.equivalent} "
                    "is only available for stress or strain results."
                )
            equivalent_op.connect(0, out)
            wf.add_operator(operator=equivalent_op)
            # Set as future output of the workflow
            out = equivalent_op.outputs.fields_container

        # Add an optional component selection step if result is vector, matrix, or principal
        if (
            category
            in [ResultCategory.vector, ResultCategory.matrix, ResultCategory.principal]
        ) and (to_extract is not None):
            # Instantiate a component selector operator
            extract_op = self._model.operator(name="component_selector_fc")
            # Feed it the current workflow output
            extract_op.connect(0, out)
            # Feed it the requested components
            extract_op.connect(1, to_extract)
            wf.add_operator(operator=extract_op)
            # Set as future output of the workflow
            out = extract_op.outputs.fields_container

        # Add an optional norm operation if requested
        if norm:
            norm_op = self._model.operator(name="norm_fc")
            norm_op.connect(0, out)
            wf.add_operator(operator=norm_op)
            out = norm_op.outputs.fields_container

        # Set the workflow output
        wf.set_output_name("out", out)
        # Evaluate  the workflow
        fc = wf.get_output("out", core.types.fields_container)

        # Test for empty results
        if (len(fc) == 0) or all([len(f) == 0 for f in fc]):
            warnings.warn(
                message=f"Returned Dataframe with columns {columns} is empty.",
                category=UserWarning,
            )
        # Return the result wrapped in a DPF_Dataframe
        return DataObject(
            fields_container=fc,
            columns=columns,
            mesh_scoping=mesh_scoping,
        )

    def displacement(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract displacement results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                Times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements whose nodes to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="U",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def velocity(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract velocity results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                Times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements whose nodes to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="V",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def acceleration(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract acceleration results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                Times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements whose nodes to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="A",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal stress results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental_nodal,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental stress results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal stress results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.nodal,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_principal(
        self,
        component_ids: Union[List[str], List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal principal stress results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental_nodal,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_principal_elemental(
        self,
        component_ids: Union[List[str], List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental principal stress results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_principal_nodal(
        self,
        component_ids: Union[List[str], List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal principal stress results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.nodal,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_eqv_von_mises(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal equivalent Von Mises stress results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental_nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_eqv_von_mises_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent Von Mises stress results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_eqv_von_mises_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent Von Mises stress results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=core.locations.elemental_nodal,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=core.locations.nodal,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=core.locations.elemental,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_principal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal principal elastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=core.locations.elemental_nodal,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_principal_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal principal elastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=core.locations.nodal,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_principal_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental principal elastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=core.locations.elemental,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_state_variable(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal plastic state variable results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_PSV",
            location=core.locations.elemental_nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_state_variable_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental plastic state variable results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_PSV",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_state_variable_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal plastic state variable results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_PSV",
            location=core.locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal plastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.elemental_nodal,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal plastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.nodal,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental plastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z", "XX", "XY",
                "XZ", or their respective equivalents 1, 2, 3, 4, 5, 6.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.elemental,
            category=ResultCategory.matrix,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_principal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal principal plastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.elemental_nodal,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_principal_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal principal plastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.nodal,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_principal_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental principal plastic strain results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.elemental,
            category=ResultCategory.principal,
            components=component_ids,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_eqv(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal equivalent plastic strain results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.elemental_nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_eqv_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent plastic strain results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_eqv_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent plastic strain results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=core.locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def reaction_force(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract reaction force results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="RF",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_volume(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental volume results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENG_VOL",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_mass(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental mass results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ElementalMass",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_heat_generation(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental heat generation results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EHC",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def element_centroids(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element centroids results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="centroids",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def thickness(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element thickness results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="thickness",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def element_orientations(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal element orientations results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=core.locations.elemental_nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def element_orientations_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental element orientations results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def element_orientations_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal element orientations results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=core.locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def artificial_hourglass_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract artificial hourglass energy results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENG_AHO",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def thermal_dissipation_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract thermal dissipation energy results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENG_TH",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def kinetic_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract kinetic energy results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENG_KE",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def hydrostatic_pressure(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract hydrostatic pressure element nodal results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_HPRES",
            location=core.locations.elemental_nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def hydrostatic_pressure_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract hydrostatic pressure nodal results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_HPRES",
            location=core.locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def hydrostatic_pressure_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract hydrostatic pressure elemental results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_HPRES",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def structural_temperature(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract structural temperature element nodal results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="BFE",
            location=core.locations.elemental_nodal,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def structural_temperature_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract structural temperature nodal results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="BFE",
            location=core.locations.nodal,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def structural_temperature_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract structural temperature elemental results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="BFE",
            location=core.locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def element_nodal_forces(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element nodal forces results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=core.locations.elemental_nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def element_nodal_forces_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element nodal forces nodal results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def element_nodal_forces_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element nodal forces elemental results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            time_steps_ids:
                List of time steps IDs to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=core.locations.elemental,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def nodal_force(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal force results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                Times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                Load steps to get results for.
            sub_steps:
                Sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements whose nodes to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="F",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def nodal_moment(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        time_step_ids: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal moment results from the simulation.

        Args:
            component_ids:
                Components to get results for. Can be one or several of "X", "Y", "Z",
                or their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                Times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                Load steps to get results for.
            sub_steps:
                Sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements whose nodes to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="M",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            time_step_ids=time_step_ids,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )
