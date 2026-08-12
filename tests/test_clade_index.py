from ete4 import Tree

from nwkit.clade_index import CladeIndex


def test_names_for_mask_visits_only_selected_bits_in_name_order():
    index = CladeIndex(Tree("((D,C),(B,A));", parser=1))

    assert index.names_for_mask(0b1010) == ("B", "D")


def test_names_for_mask_ignores_bits_outside_the_index():
    index = CladeIndex(Tree("(A,B,C);", parser=1))

    assert index.names_for_mask(0b10101) == ("A", "C")
    assert index.names_for_mask(-2) == ("B", "C")
