#!/usr/bin/env python3
"""
Minimal Newick parser + simple cladogram/phylogram plotter using only
matplotlib (no ete3/biopython dependency).

Usage: python3 plot_tree.py fad2_tree.nwk fad2_tree.png
"""
import sys
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
        # optional name
        m = re.match(r"[^,:()]*", s[pos:])
        if m:
            name = m.group(0)
            pos += len(name)
            if name:
                node.name = name
        # optional branch length
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[-0-9.eE+]+", s[pos:])
            length = float(m.group(0))
            pos += len(m.group(0))
            node.length = length
        return node

    root = parse_node()
    return root


def get_leaves(node):
    if node.is_leaf():
        return [node]
    leaves = []
    for c in node.children:
        leaves.extend(get_leaves(c))
    return leaves


def assign_positions(node, depth_acc, y_counter, y_positions):
    """Assign x = cumulative branch length from root, y = leaf order."""
    x = depth_acc + node.length
    if node.is_leaf():
        y = y_counter[0]
        y_counter[0] += 1
        y_positions[id(node)] = (x, y)
        return y
    ys = []
    for c in node.children:
        ys.append(assign_positions(c, x, y_counter, y_positions))
    y = sum(ys) / len(ys)
    y_positions[id(node)] = (x, y)
    return y


def plot(root, out_png, title="FAD2 protein tree (FastTree, JTT+CAT)"):
    y_positions = {}
    assign_positions(root, 0.0, [0], y_positions)

    fig, ax = plt.subplots(figsize=(7, 3.5))

    def draw(node):
        x, y = y_positions[id(node)]
        if node.parent is not None:
            px, py = y_positions[id(node.parent)]
            # horizontal branch
            ax.plot([px, x], [y, y], color="black", lw=1.5)
        if not node.is_leaf():
            child_ys = [y_positions[id(c)][1] for c in node.children]
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="black", lw=1.5)
            for c in node.children:
                draw(c)
        else:
            ax.text(x + 0.005, y, node.name, va="center", ha="left", fontsize=9)

    draw(root)
    ax.set_xlabel("Substitutions per site")
    ax.set_yticks([])
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: plot_tree.py tree.nwk out.png")
    with open(sys.argv[1]) as fh:
        newick = fh.read()
    root = parse_newick(newick)
    plot(root, sys.argv[2])
