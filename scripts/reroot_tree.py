#!/usr/bin/env python3
"""
Reroot a Newick tree on a named outgroup leaf.

FastTree always outputs an unrooted tree (drawn with an arbitrary
trifurcating "root" that has no biological meaning); IQ-TREE can be told
an outgroup directly (-o flag) and roots its own output accordingly, but
FastTree has no equivalent option. This script does the same job for a
FastTree (or any) Newick file: standard "outgroup rooting" -- insert a new
root on the branch leading to the named leaf, splitting that branch's
length in half, and reverse the parent/child direction of every edge on
the path back to the original (arbitrary) root so the rest of the tree
hangs correctly off the new root.

Internal node labels (branch support values, e.g. FastTree's SH-like
support or a bootstrap percentage) are preserved: each label stays
attached to the same edge/bipartition it described before rerooting,
since rerooting an unrooted tree doesn't change what any internal edge
separates -- only which end is drawn as "up".

Usage: reroot_tree.py in.nwk out.nwk <outgroup_substring>
  <outgroup_substring> matches against leaf names (e.g. an accession like
  "P46313" is enough; it doesn't need to be the full header).
"""
import sys
import re


class Node:
    def __init__(self, name=None, length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None

    def is_leaf(self):
        return len(self.children) == 0


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if s[pos] == ",":
                    pos += 1
                    continue
                elif s[pos] == ")":
                    pos += 1
                    break
        m = re.match(r"[^,:()]*", s[pos:])
        if m:
            name = m.group(0)
            pos += len(name)
            if name:
                node.name = name
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[-0-9.eE+]+", s[pos:])
            length = float(m.group(0))
            pos += len(m.group(0))
            node.length = length
        return node

    return parse_node()


def to_newick(node, is_root=True):
    if node.is_leaf():
        s = node.name or ""
    else:
        s = "(" + ",".join(to_newick(c, is_root=False) for c in node.children) + ")"
        if node.name:
            s += node.name
    if not is_root:
        s += f":{node.length:.6f}"
    return s


def collect_nodes(node, out):
    out.append(node)
    for c in node.children:
        collect_nodes(c, out)


def reroot(root0, outgroup_substr):
    all_nodes = []
    collect_nodes(root0, all_nodes)

    # id -> label (leaf name, or internal support label if any)
    label = {id(n): n.name for n in all_nodes}
    is_leaf = {id(n): n.is_leaf() for n in all_nodes}

    # Undirected adjacency built from the original parent/child edges.
    adj = {id(n): [] for n in all_nodes}
    for n in all_nodes:
        if n.parent is not None:
            adj[id(n)].append((id(n.parent), n.length))
            adj[id(n.parent)].append((id(n), n.length))

    id_to_node = {id(n): n for n in all_nodes}

    matches = [n for n in all_nodes if n.is_leaf() and n.name and outgroup_substr in n.name]
    if not matches:
        sys.exit(f"ERROR: no leaf name contains '{outgroup_substr}'")
    if len(matches) > 1:
        sys.exit(f"ERROR: '{outgroup_substr}' matches multiple leaves: "
                  f"{[m.name for m in matches]} -- use a more specific substring")
    og = matches[0]
    og_id = id(og)

    neighbors = adj[og_id]
    if len(neighbors) != 1:
        sys.exit(f"ERROR: outgroup leaf '{og.name}' does not have exactly one neighbor "
                  f"(got {len(neighbors)}) -- tree may be malformed")
    p_id, e_len = neighbors[0]

    new_root = Node(name=None, length=0.0)
    og_new = Node(name=label[og_id], length=e_len / 2.0)
    p_new = Node(name=label[p_id] if not is_leaf[p_id] else label[p_id], length=e_len / 2.0)
    new_root.children = [og_new, p_new]
    og_new.parent = new_root
    p_new.parent = new_root

    def build(node_obj, current_id, came_from_id):
        for nbr_id, length in adj[current_id]:
            if nbr_id == came_from_id:
                continue
            child_obj = Node(name=label[nbr_id], length=length)
            child_obj.parent = node_obj
            node_obj.children.append(child_obj)
            build(child_obj, nbr_id, current_id)

    build(p_new, p_id, og_id)
    return new_root


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: reroot_tree.py in.nwk out.nwk <outgroup_substring>")
    in_path, out_path, outgroup_substr = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(in_path) as fh:
        newick = fh.read()
    original_root = parse_newick(newick)
    rerooted = reroot(original_root, outgroup_substr)
    with open(out_path, "w") as fh:
        fh.write(to_newick(rerooted, is_root=True) + ";\n")
    print(f"wrote {out_path} (rerooted on outgroup matching '{outgroup_substr}')")
