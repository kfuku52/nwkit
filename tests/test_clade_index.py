from ete4 import Tree

from nwkit.clade_index import CladeIndex


def test_names_for_mask_visits_only_selected_bits_in_name_order():
    index = CladeIndex(Tree("((D,C),(B,A));", parser=1))

    assert index.names_for_mask(0b1010) == ("B", "D")


def test_names_for_mask_ignores_bits_outside_the_index():
    index = CladeIndex(Tree("(A,B,C);", parser=1))

    assert index.names_for_mask(0b10101) == ("A", "C")
    assert index.names_for_mask(-2) == ("B", "C")


def test_names_for_wide_mask_preserves_order_and_ignores_outside_bits():
    names = tuple(f"T{i:03d}" for i in range(600))
    index = CladeIndex(Tree(f"({','.join(reversed(names))});", parser=1))
    selected = tuple(i for i in range(len(names)) if i % 7 < 3)
    mask = sum(1 << i for i in selected) | (1 << 900)
    assert index.names_for_mask(mask) == tuple(names[i] for i in selected)
    assert index.names_for_mask(-2) == names[1:]
