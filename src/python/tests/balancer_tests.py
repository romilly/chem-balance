from unittest import TestCase
from hamcrest import assert_that, equal_to, close_to
from chemistry.balancer import Compound, balance


class CompoundTest(TestCase):
    def test_parses_compound(self):
        self.check_composition('H2SO4', {'H': 2, 'S': 1, 'O': 4})
        self.check_composition('SeS2', {'Se': 1, 'S': 2})
        self.check_composition('Ca(C2H3O2)2', {'Ca': 1, 'C': 4, 'H': 6, 'O': 4})

    def test_can_calculate_molecular_weight(self):
        assert_that(Compound('H2SO4').molecular_weight(), is_about(98.072))
        assert_that(Compound('Ca(C2H3O2)2').molecular_weight(), is_about(158.166))

    def check_composition(self, formula, composition):
        assert_that(Compound(formula).composition(), equal_to(composition))

class BalancerTest(TestCase):
    def test_balances_feasible_equation(self):
        # assert_that(balance('Fe+Cl2=FeCl3'), equal_to('2Fe+3Cl2=2FeCl3'))
        assert_that(balance('C2H4+O2=CO2+H2O'), equal_to('C2H4+3O2=2CO2+2H2O'))


def is_about(value):
    return close_to(value, value * 0.0001)
