from unittest import TestCase
from hamcrest import assert_that, equal_to, close_to
from balancer import atoms_in, molecular_weight



class AtomTest(TestCase):
    def test_parses_compound(self):
        self.check_composition('H2SO4', {'H': 2, 'S': 1, 'O': 4})
        self.check_composition('SeS2', {'Se': 1, 'S': 2})
        self.check_composition('Ca(C2H3O2)2', {'Ca': 1, 'C': 4, 'H': 6, 'O': 4})

    def test_can_calculate_molecular_weight(self):
        assert_that(molecular_weight('H2SO4'), is_about(98.072))
        assert_that(molecular_weight('Ca(C2H3O2)2'), is_about(158.166))

    def check_composition(self, compound, composition):
        assert_that(dict(atoms_in(compound)), equal_to(composition, ))


def is_about(value):
    return close_to(value, value * 0.0001)
