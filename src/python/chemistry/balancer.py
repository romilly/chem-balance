"""
Determine the composition of a compound from its molecular formula, calculate its molecular weight and balance chemical equations.

This module implements:
Compound - a class that is created from a formula and knows its constituent elements, their quantities and its molecular weight.
amu - a function that returns the atomic mass of a given element
balance - a function that returns the balanced version of a chemical equation
"""
import re
from collections import defaultdict
from fractions import Fraction
from math import gcd
from functools import reduce

"""
regular expressions used in parsing chemical formulae and equations
"""
eq = '[A-Z][a-z]?[\d]*'
expr_form = re.compile('(\((?:%s)*\)[\d]*|%s)' % (eq,eq))
bracketed_form = re.compile('\(((?:%s)*)\)([\d]*)' % eq)
compound_form = re.compile('([A-Z][a-z]?)([\d]*)')

"""
Compound - a class that is created from a formula and knows its constituent elements, their quantities and its molecular weight.

The formula should consist of a series of elements and quantities, as in H2SO4, CH4, Ca2CO3.
The formula can also include bracketed items: CH3(CH2)50CH3
"""


class Compound():
    def __init__(self, formula):
        self._formula = formula
        self._elements = self._atoms_in(formula)

    def formula(self):
        return self._formula

    def elements(self):
        return self._elements.keys()

    def __getitem__(self, item):
        return self._elements[item]

    def __setitem__(self, key, value):
        self._elements[key] = value

    def _atoms_in(self, formula):
        dict = defaultdict(int)
        for expr in expr_form.findall(formula):
            if expr.startswith('('):
                fm = bracketed_form.match(expr)
                chain = fm.group(1)
                factor = count(fm.group(2))
                m = Compound(chain)
                for element in m.elements():
                    dict[element] += factor *m[element]
            else:
                for (element, qty) in compound_form.findall(expr):
                    dict[element] += count(qty)
        return dict

    def composition(self):
        return self._elements

    def molecular_weight(self):
        composition = self.composition()
        return sum([amu(element) * composition[element] for element in self.elements()])


def count(string):
    return int(string) if string else 1


def amu(element):
    return edict[element][2]


def all_zero(row):
    for item in row:
        if item != 0:
            return False
    return True

# lcm modified from https://gist.github.com/endolith/114336


def lcm(*numbers):
    """Return lowest common multiple."""
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, *numbers)


class BalancerException(Exception):
    pass


def balance(equation):
    compounds = compounds_in(equation)
    elements = elements_in(compounds)
    m = composition_matrix(compounds, elements)
    m_factors, p_factors = find_balancing_factors(m)
    reagents, products = rp(equation)
    br = balanced(p_factors, reagents)
    pr = balanced(m_factors, products)
    return br+'='+pr


def rp(equation):
    reagents, products = equation.split('=')
    reagents = reagents.split('+')
    products = products.split('+')
    return reagents, products


def find_balancing_factors( m):
    me = catenate(transpose(m), identity(len(m[0])))
    ref = reduced_row_echelon_form(me)
    zrows = zeroish_rows(len(m), ref)
    if len(zrows) == 0:
        raise BalancerException('No solutions')
    if len(zrows) > 1:
        raise BalancerException('Multiple Solutions')
    m_factors, p_factors = find_factors(zrows[0])
    return m_factors, p_factors


def zeroish_rows(n, rows):
    return [row[n:] for row in rows if all_zero(row[:n])]


def compounds_in(equation):
    return [Compound(formula) for formula in re.split('\+|=', equation)]


def elements_in(compounds):
    elements = set()
    for compound in compounds:
        for element in compound.elements():
            elements.add(element)
    return sorted(elements)


def transpose(m):
    return [list(row) for row in zip(*m)]


def composition_matrix(compounds, elements):
    return [[Fraction(compound[element]) for compound in compounds] for element in elements]


def find_factors(row):
    d = lcm([factor.denominator for factor in row])
    i_f = [factor*d for factor in row]
    return (negated(i_f)), (positive(i_f))


def negated(f):
    return [abs(factor) for factor in f if factor < 0]


def positive(f):
    return [factor for factor in f if factor > 0]


def balanced(factors, compounds):
    return '+'.join([prefix(factor) + compound for (factor, compound) in
                     zip(factors, compounds)])

def prefix(factor):
    return str(factor) if factor > 1 else ''


# from https://rosettacode.org/wiki/Reduced_row_echelon_form#Python
# adapted and additional functions added


def reduced_row_echelon_form(m):
    if not m: return
    reduced_matrix = [[m[i][j] for j in range(len(m[0]))] for i in range(len(m))]
    lead = 0
    rowCount = len(reduced_matrix)
    columnCount = len(reduced_matrix[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return reduced_matrix
        i = r
        while reduced_matrix[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return reduced_matrix
        reduced_matrix[i], reduced_matrix[r] = reduced_matrix[r], reduced_matrix[i]
        lv = reduced_matrix[r][lead]
        reduced_matrix[r] = [mrx / lv for mrx in reduced_matrix[r]]
        for i in range(rowCount):
            if i != r:
                lv = reduced_matrix[i][lead]
                reduced_matrix[i] = [iv - lv * rv for rv, iv in zip(reduced_matrix[r], reduced_matrix[i])]
        lead += 1
    return reduced_matrix

def identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def catenate(m, n):
    return [m[i]+n[i] for i in range(len(m))]

# taken from https://en.wikipedia.org/wiki/Standard_atomic_weight

edict = {
    'H' : (1, 'hydrogen', 1.008),
    'He' : (2, 'helium', 4.0026),
    'Li' : (3, 'lithium', 6.94),
    'Be' : (4, 'beryllium', 9.0122),
    'B' : (5, 'boron', 10.81),
    'C' : (6, 'carbon', 12.011),
    'N' : (7, 'nitrogen', 14.007),
    'O' : (8, 'oxygen', 15.999),
    'F' : (9, 'fluorine', 18.998),
    'Ne' : (10, 'neon', 20.180),
    'Na' : (11, 'sodium', 22.990),
    'Mg' : (12, 'magnesium', 24.305),
    'Al' : (13, 'aluminium', 26.982),
    'Si' : (14, 'silicon', 28.085),
    'P' : (15, 'phosphorus', 30.974),
    'S' : (16, 'sulfur', 32.06),
    'Cl' : (17, 'chlorine', 35.45),
    'Ar' : (18, 'argon', 39.948),
    'K' : (19, 'potassium', 39.098),
    'Ca' : (20, 'calcium', 40.078),
    'Sc' : (21, 'scandium', 44.956),
    'Ti' : (22, 'titanium', 47.867),
    'V' : (23, 'vanadium', 50.942),
    'Cr' : (24, 'chromium', 51.996),
    'Mn' : (25, 'manganese', 54.938),
    'Fe' : (26, 'iron', 55.845),
    'Co' : (27, 'cobalt', 58.933),
    'Ni' : (28, 'nickel', 58.693),
    'Cu' : (29, 'copper', 63.546),
    'Zn' : (30, 'zinc', 65.38),
    'Ga' : (31, 'gallium', 69.723),
    'Ge' : (32, 'germanium', 72.630),
    'As' : (33, 'arsenic', 74.922),
    'Se' : (34, 'selenium', 78.971),
    'Br' : (35, 'bromine', 79.904),
    'Kr' : (36, 'krypton', 83.798),
    'Rb' : (37, 'rubidium', 85.468),
    'Sr' : (38, 'strontium', 87.62),
    'Y' : (39, 'yttrium', 88.906),
    'Zr' : (40, 'zirconium', 91.224),
    'Nb' : (41, 'niobium', 92.906),
    'Mo' : (42, 'molybdenum', 95.95),
    'Ru' : (44, 'ruthenium', 101.07),
    'Rh' : (45, 'rhodium', 102.91),
    'Pd' : (46, 'palladium', 106.42),
    'Ag' : (47, 'silver', 107.87),
    'Cd' : (48, 'cadmium', 112.41),
    'In' : (49, 'indium', 114.82),
    'Sn' : (50, 'tin', 118.71),
    'Sb' : (51, 'antimony', 121.76),
    'Te' : (52, 'tellurium', 127.60),
    'I' : (53, 'iodine', 126.90),
    'Xe' : (54, 'xenon', 131.29),
    'Cs' : (55, 'caesium', 132.91),
    'Ba' : (56, 'barium', 137.33),
    'La' : (57, 'lanthanum', 138.91),
    'Ce' : (58, 'cerium', 140.12),
    'Pr' : (59, 'praseodymium', 140.91),
    'Nd' : (60, 'neodymium', 144.24),
    'Sm' : (62, 'samarium', 150.36),
    'Eu' : (63, 'europium', 151.96),
    'Gd' : (64, 'gadolinium', 157.25),
    'Tb' : (65, 'terbium', 158.93),
    'Dy' : (66, 'dysprosium', 162.50),
    'Ho' : (67, 'holmium', 164.93),
    'Er' : (68, 'erbium', 167.26),
    'Tm' : (69, 'thulium', 168.93),
    'Yb' : (70, 'ytterbium', 173.05),
    'Lu' : (71, 'lutetium', 174.97),
    'Hf' : (72, 'hafnium', 178.49),
    'Ta' : (73, 'tantalum', 180.95),
    'W' : (74, 'tungsten', 183.84),
    'Re' : (75, 'rhenium', 186.21),
    'Os' : (76, 'osmium', 190.23),
    'Ir' : (77, 'iridium', 192.22),
    'Pt' : (78, 'platinum', 195.08),
    'Au' : (79, 'gold', 196.97),
    'Hg' : (80, 'mercury', 200.59),
    'Tl' : (81, 'thallium', 204.38),
    'Pb' : (82, 'lead', 207.2),
    'Bi' : (83, 'bismuth', 208.98),
    'Th' : (90, 'thorium', 232.04),
    'Pa' : (91, 'protactinium', 231.04),
    'U' : (92, 'uranium', 238.03)}













