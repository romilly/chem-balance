import re
from collections import defaultdict

# becomes = re.compile('->')
expr_form = re.compile('(\((?:[A-Z][a-z]?[\d]*)*\)[\d]*|[A-Z][a-z]?[\d]*)')
compound_form = re.compile('([A-Z][a-z]?)([\d]*)')
bracketed_form = re.compile('\(((?:[A-Z]|[a-z]|[\d])*)\)([\d]*)')

def atoms_in(compound):
    dict = defaultdict(int)
    for expr in expr_form.findall(compound):
        print(expr)
        if expr.startswith('('):
            m = bracketed_form.match(expr)
            chain = m.group(1)
            factor = int(m.group(2)) if m.group(2) else 1
            d = atoms_in(chain)
            for key in d.keys():
                dict[key] += factor * d[key]
        else:
            for (element, qty) in compound_form.findall(expr):
                dict[element] += int(qty) if qty else 1
    return dict

# (lhs, rhs) = becomes.split('Fe2O3+SO2->CuO3+H2O+CO2')
# dod = {}
# for item in items:
#     dod[item] = atoms_in(item)
# print(dod)
# print (lhs, rhs)
