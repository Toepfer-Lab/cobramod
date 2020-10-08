from typing import Generator, Union, Dict, NamedTuple


class MetaboliteTuple(NamedTuple):
    identifier: str
    coefficient: float


def _get_keys(raw: str) -> Generator:
    """
    Returns as Generator the keys for KEGGs data. Input is raw text
    """
    lines = (line for line in raw.split(sep="\n"))
    for line in lines:
        segment = line.split(" ")
        if segment[0] == "///":
            continue
        if segment[0] != "":
            yield segment[0]


def _find_key(line: str, keys: set) -> Union[str, None]:
    """
    Return key if found in given line. Else, None is returned.
    """
    for key in keys:
        if key in line:
            return key
    return None


def _create_dict(raw: str) -> dict:
    """
    Formats most of the keys for Keggs data.
    """
    lines = (line for line in raw.split("\n"))
    keys = set(_get_keys(raw=raw))
    actual_key = str()
    kegg_dict: Dict[str, list] = {"": []}
    for line in lines:
        key = _find_key(line=line, keys=keys)
        if key is None:
            # Append to older key.
            kegg_dict[actual_key].append(line.strip().rstrip())
        else:
            actual_key = key
            line = line[line.find(key)+len(key):].strip().rstrip()
            kegg_dict[actual_key] = [line]
    del kegg_dict[""]
    return kegg_dict


def _get_coefficient(line: str, SIDE: int) -> MetaboliteTuple:
    """
    Returns a NamedTuple with the metabolite identifier and coefficient
    that appears in given line.

    Args:
        line (str): string with information
        SIDE (int): Constant to multiply depending on the side

    Returns:
        NamedTuple: tupple with identifier and coefficient
    """
    # " 4 H+ "
    line = line.strip().rstrip().split(" ")
    try:
        return MetaboliteTuple(
            identifier=line[1], coefficient=float(line[0]) * SIDE)
    except IndexError:
        return MetaboliteTuple(identifier=line[0], coefficient=1.0 * SIDE)


def _give_metabolites(line: str) -> dict:
    """
    Returns dictionary with metabolites identifiers and corresponding
    coeffient
    """
    middle = line.find("=")
    temp_dict = dict()
    reactants = line[:middle-1]
    products = line[middle+2:]
    for item in reactants.split(" + "):
        MetaTuple = _get_coefficient(line=item, SIDE=-1)
        temp_dict[MetaTuple.identifier] = MetaTuple.coefficient
    for item in products.split(" + "):
        MetaTuple = _get_coefficient(line=item, SIDE=1)
        temp_dict[MetaTuple.identifier] = MetaTuple.coefficient
    return temp_dict


def _get_reversibility(line: str) -> tuple:
    """
    Returns the bounds for reaction depending of the string
    """
    # FIXME: Direction depends also from extra keys
    middle = line.find("=")
    line = line[middle-1:middle+2]
    if line == "<=>":
        bounds = (-1000, 1000)
    elif line == "==>":
        bounds = (0, 1000)
    elif line == "<==":
        bounds = (-1000, 0)
    else:
        raise Warning(
            f'"Reversibility for "{line}" could not be found')
    return bounds
