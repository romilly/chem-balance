# chem-balance

A collection of Python Scripts for Chemistry Students.

The code parses molecular formulae, calculates molecular weights
and can balance chemical equations.

The code has been tested with Python version 3.5. 
Changes to the `math` and `fractions` packages
mean that it may not work with earlier versions.

The tests use PyHamcrest which is the only external requirement.

I will document the code but these examples may be enough to get you going:

    from balancer import Compound, balance 
    
    print(Compound('H2SO4').molecular_weight())
    
    print(balance('Fe+Cl2=FeCl3'))
