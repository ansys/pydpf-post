from ansys.dpf.post.enums import ResultCategory

component_label_to_index = {
    "1": 0,
    "2": 1,
    "3": 2,
    "4": 3,
    "5": 4,
    "6": 5,
    "X": 0,
    "Y": 1,
    "Z": 2,
    "XX": 0,
    "YY": 1,
    "ZZ": 2,
    "XY": 3,
    "YZ": 4,
    "XZ": 5,
}

vector_component_names = ["X", "Y", "Z"]
matrix_component_names = ["XX", "YY", "ZZ", "XY", "YZ", "XZ"]
principal_names = ["1", "2", "3"]


def build_components_for_vector(base_name, components):
    out, columns = build_components(base_name, components, vector_component_names)
    return out, columns


def build_components_for_matrix(base_name, components):
    out, columns = build_components(base_name, components, matrix_component_names)
    return out, columns


def build_components(base_name, components, component_names):
    # Create operator internal names based on components
    out = []
    if components is None:
        out = None
    else:
        if isinstance(components, int) or isinstance(components, str):
            components = [components]
        if not isinstance(components, list):
            raise ValueError(
                "Argument 'components' must be an int, a str, or a list of either."
            )
        for comp in components:
            if not (isinstance(comp, str) or isinstance(comp, int)):
                raise ValueError(
                    "Argument 'components' can only contain integers and/or strings.\n"
                    f"The provided component '{comp}' is not valid."
                )
            if isinstance(comp, int):
                comp = str(comp)
            if comp not in component_label_to_index.keys():
                raise ValueError(
                    f"Component {comp} is not valid. Please use one of: "
                    f"{list(component_label_to_index.keys())}."
                )
            out.append(component_label_to_index[comp])

    # Take unique values and build names list
    if out is None:
        columns = [base_name + comp for comp in component_names]
    else:
        out = list(set(out))
        columns = [base_name + component_names[i] for i in out]
    return out, columns


def build_components_for_principal(base_name, components):
    # Create operator internal names based on principal components
    out = []
    if components is None:
        components = [1]

    if isinstance(components, int) or isinstance(components, str):
        components = [components]
    if not isinstance(components, list):
        raise ValueError(
            "Argument 'components' must be an int, a str, or a list of either."
        )
    for comp in components:
        if not (isinstance(comp, str) or isinstance(comp, int)):
            raise ValueError(
                "Argument 'components' can only contain integers and/or strings."
            )
        if str(comp) not in principal_names:
            raise ValueError(
                "A principal component ID must be one of: " f"{principal_names}."
            )
        out.append(int(comp) - 1)

    # Take unique values
    if out is not None:
        out = list(set(out))
    # Build columns names
    if out is None:
        columns = [base_name + str(comp) for comp in principal_names]
    else:
        columns = [base_name + principal_names[i] for i in out]
    return out, columns


def create_components(base_name: str, category: ResultCategory, components):
    comp = None
    # Build the list of requested results
    if category in [ResultCategory.scalar, ResultCategory.equivalent]:
        # A scalar or equivalent result has no components
        to_extract = None
        columns = [base_name]
    elif category == ResultCategory.vector:
        # A vector result can have components selected
        to_extract, columns = build_components_for_vector(
            base_name=base_name, components=components
        )
        if to_extract is not None:
            comp = [vector_component_names[i] for i in to_extract]
        else:
            comp = vector_component_names
    elif category == ResultCategory.matrix:
        # A vector result can have components selected
        to_extract, columns = build_components_for_matrix(
            base_name=base_name, components=components
        )
        if to_extract is not None:
            comp = [matrix_component_names[i] for i in to_extract]
        else:
            comp = matrix_component_names
    elif category == ResultCategory.principal:
        # A principal type of result can have components selected
        to_extract, columns = build_components_for_principal(
            base_name=base_name, components=components
        )
        comp = [principal_names[i] for i in to_extract]
    else:
        raise ValueError(f"'{category}' is not a valid category value.")
    return comp, to_extract, columns
