import re
from collections import defaultdict
from aw import edict

eq = '[A-Z][a-z]?[\d]*'
expr_form = re.compile('(\((?:%s)*\)[\d]*|%s)' % (eq,eq))
compound_form = re.compile('([A-Z][a-z]?)([\d]*)')
bracketed_form = re.compile('\(((?:[A-Z]|[a-z]|[\d])*)\)([\d]*)')

def count(string):
    return int(string) if string else 1

def amu(element):
    return edict[element][2]

def atoms_in(compound):
    dict = defaultdict(int)
    for expr in expr_form.findall(compound):
        if expr.startswith('('):
            m = bracketed_form.match(expr)
            chain = m.group(1)
            factor = count(m.group(2))
            d = atoms_in(chain)
            for key in d.keys():
                dict[key] += factor * d[key]
        else:
            for (element, qty) in compound_form.findall(expr):
                dict[element] += count(qty)
    return dict

def molecular_weight(compound):
    quantities = atoms_in(compound)
    return sum([amu(element) * quantities[element] for element in quantities.keys()])

