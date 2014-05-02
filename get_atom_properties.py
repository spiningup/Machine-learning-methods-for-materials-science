import pylada.periodic_table as pt


sanderson = {}
pauling = {}
radius = {}

for symbol in pt.symbols:
    sanderson[symbol] = getattr(pt, symbol).sanderson
    pauling[symbol] = getattr(pt, symbol).pauling
    try:
        radius[symbol] = getattr(pt, symbol).atomic_radius.magnitude
    except:
        pass

print sanderson
print pauling
print radius
