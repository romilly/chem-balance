import re
from collections import defaultdict
from chemistry.aw import edict
from maths.reduced_echelon import reduced_row_echelon_form, catenate, identity

eq = '[A-Z][a-z]?[\d]*'
expr_form = re.compile('(\((?:%s)*\)[\d]*|%s)' % (eq,eq))
compound_form = re.compile('([A-Z][a-z]?)([\d]*)')
bracketed_form = re.compile('\(((?:%s)*)\)([\d]*)' % eq)


class Compound():
    def __init__(self, formula):
        self._formula = formula
        self._elements = self.atoms_in(formula)

    def formula(self):
        return self._formula

    def elements(self):
        return self._elements.keys()

    def __getitem__(self, item):
        return self._elements[item]

    def __setitem__(self, key, value):
        self._elements[key] = value

    def atoms_in(self, formula):
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
        if abs(item) > 1E-4:
            return False
    return True


def balance(equation):
    compounds = [Compound(formula) for formula in re.split('\+|=', equation)]
    elements = set()
    for compound in compounds:
        for element in compound.elements():
            elements.add(element)
    m = [[compound[element] for compound in compounds] for element in elements]
    mt = list((zip(*m)))
    me = catenate(mt, identity(len(mt)))
    ref = reduced_row_echelon_form(me)
    zrows = [row for row in ref if all_zero(row[:len(m)])]
    if len(zrows) != 1:
        return
    factors = list(map(lambda x: int(abs(x)), zrows[0]))[len(m):]
    reagents, products = equation.split('=')
    reagents = reagents.split('+')
    products = products.split('+')
    br = '+'.join([(str(factor) if factor > 1 else '')+reagent for (factor, reagent) in zip(factors[:len(reagents)], reagents)])
    pr = '+'.join([(str(factor) if factor > 1 else '')+product for (factor, product) in zip(factors[-len(products):],products)])
    return br+'='+pr











