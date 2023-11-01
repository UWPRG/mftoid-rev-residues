import mdtraj as md
import numpy as np
import networkx as nx


def make_bondgraph(top):
    """Returns a bond graph from topology"""
    G = nx.Graph()
    G.add_nodes_from(top.atoms)
    G.add_edges_from(top.bonds)
    element_dict = {atom: atom.element.symbol for atom in G.nodes}
    nx.set_node_attributes(G, element_dict, "element")
    return G


def get_pdb_compare_dict(pdb_fname):
    """Parses PDB file"""
    compare = dict()
    idx_counter = 0
    with open(pdb_fname) as f:
        for line in f:
            split_line = line.split()
            if (split_line[0] == "HETATM") or (split_line[0] == "ATOM"):
                if "C:" in split_line[-1]:
                    connection = split_line[-1].split("C:")[-1]
                else:
                    connection = None
                compare[idx_counter] = connection
                idx_counter += 1
    return compare


def parse_pdb(pdb_fname):
    """Returns a graph object with available connectivity information"""
    trj = md.load(pdb_fname).center_coordinates()
    G = make_bondgraph(trj.top)
    compare_dict = get_pdb_compare_dict(pdb_fname)

    if len(compare_dict) != trj.top.n_atoms:
        raise ValueError(
            f"Error reading {pdb_fname}. Connect dict ({len(compare_dict)}) != {trj.top.n_atoms} atoms"
        )

    compare_dict = {trj.top.atom(k): v for k, v in compare_dict.items()}
    nx.set_node_attributes(G, compare_dict, "compare")
    return G


def rotate_vec(a, b):
    # rotate a onto b
    def skew(x):
        return np.array([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])

    v = np.cross(a, b)
    c = np.dot(a, b)
    v_skew = skew(v)

    R = np.eye(3) + v_skew + np.dot(v_skew, v_skew) / (1.0 + c)

    return R.T


def subgraph_match(G1, G2):
    "Check if G2 is in G1"

    def element_match(n1, n2):
        if n1["element"] == n2["element"]:
            return True
        return False

    GM = nx.algorithms.isomorphism.GraphMatcher(G2, G1, element_match)

    if GM.subgraph_is_isomorphic():
        return list(GM.subgraph_isomorphisms_iter())
    else:
        raise ValueError("No matching subgraphs")


def new_line(line, idx):
    l = line.split()
    l.append((str(idx)))
    return "".join(i + "   " for i in l) + "\n"
