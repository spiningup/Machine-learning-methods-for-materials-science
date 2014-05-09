import pylada.periodic_table as pt


sanderson = {}
pauling = {}
radius = {}
row = {}
column = {}
Eionization = {}
atomic_number = {}
Eaffinity = {}

for symbol in pt.symbols:
    sanderson[symbol] = getattr(pt, symbol).sanderson
    pauling[symbol] = getattr(pt, symbol).pauling
    row[symbol] = getattr(pt, symbol).row
    column[symbol] = getattr(pt, symbol).column
    atomic_number[symbol] = getattr(pt, symbol).atomic_number
    try:
        radius[symbol] = getattr(pt, symbol).atomic_radius.magnitude
    except:
        pass

    try:
        Eionization[symbol] = float('%6.2f'%(getattr(pt, symbol).ionization_energies[0].magnitude / 1e6))
    except:
        pass
    
    try:     
        Eaffinity[symbol] = float('%6.2f'%(getattr(pt, symbol).electron_affinity.magnitude / 1e3))
    except:
        pass

print sanderson
print pauling
print radius
print row
print column
print Eionization
print atomic_number
print Eaffinity
