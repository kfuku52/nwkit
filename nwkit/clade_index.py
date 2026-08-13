import csv
import hashlib
from io import StringIO

CLADE_ID_PREFIX = "clade-sha256:"


class CladeIndex:
    """Compact descendant-tip index with child-order-independent clade IDs."""

    def __init__(self, tree):
        self.names = tuple(sorted(str(leaf.name) for leaf in tree.leaves()))
        bit_by_name = {name: 1 << index for index, name in enumerate(self.names)}
        self.mask_by_node = dict()
        for node in tree.traverse(strategy="postorder"):
            if node.is_leaf:
                self.mask_by_node[node] = bit_by_name[str(node.name)]
            else:
                mask = 0
                for child in node.children:
                    mask |= self.mask_by_node[child]
                self.mask_by_node[node] = mask
        self._clade_id_by_node = dict()

    def names_for_mask(self, mask):
        names = []
        remaining = int(mask) & ((1 << len(self.names)) - 1)
        while remaining:
            lowest_bit = remaining & -remaining
            names.append(self.names[lowest_bit.bit_length() - 1])
            remaining ^= lowest_bit
        return tuple(names)

    def names_for_node(self, node):
        return self.names_for_mask(self.mask_by_node[node])

    def csv_for_mask(self, mask):
        output = StringIO()
        csv.writer(output, lineterminator="").writerow(self.names_for_mask(mask))
        return output.getvalue()

    def csv_for_node(self, node):
        return self.csv_for_mask(self.mask_by_node[node])

    def count_for_mask(self, mask):
        return int(mask).bit_count()

    def count_for_node(self, node):
        return self.count_for_mask(self.mask_by_node[node])

    def clade_id_for_node(self, node):
        if node not in self._clade_id_by_node:
            digest = hashlib.sha256(b"nwkit-clade-v1\0")
            for name in self.names_for_node(node):
                encoded = name.encode("utf-8")
                digest.update(len(encoded).to_bytes(8, byteorder="big"))
                digest.update(encoded)
            self._clade_id_by_node[node] = CLADE_ID_PREFIX + digest.hexdigest()
        return self._clade_id_by_node[node]


class LcaIndex:
    """Binary-lifting LCA index for repeated reconciliation queries."""

    def __init__(self, tree):
        self.nodes = tuple(tree.traverse(strategy="preorder"))
        self.index_by_node = {node: index for index, node in enumerate(self.nodes)}
        parent = list()
        depth = list()
        for node in self.nodes:
            if node.is_root:
                parent.append(self.index_by_node[node])
                depth.append(0)
            else:
                parent.append(self.index_by_node[node.up])
                depth.append(depth[self.index_by_node[node.up]] + 1)
        self.depth = depth
        self.ancestors = [parent]
        max_levels = max(1, max(depth, default=0).bit_length())
        for _ in range(1, max_levels):
            previous = self.ancestors[-1]
            self.ancestors.append([previous[value] for value in previous])

    def common_ancestor(self, node1, node2):
        index1 = self.index_by_node[node1]
        index2 = self.index_by_node[node2]
        if self.depth[index1] < self.depth[index2]:
            index1, index2 = index2, index1
        depth_difference = self.depth[index1] - self.depth[index2]
        for level, ancestors in enumerate(self.ancestors):
            if depth_difference & (1 << level):
                index1 = ancestors[index1]
        if index1 == index2:
            return self.nodes[index1]
        for ancestors in reversed(self.ancestors):
            if ancestors[index1] != ancestors[index2]:
                index1 = ancestors[index1]
                index2 = ancestors[index2]
        return self.nodes[self.ancestors[0][index1]]
