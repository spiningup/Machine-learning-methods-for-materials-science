import numpy as np
from read_json import read_json
from atomic_constants import pauling, radius, Eatom, Emadelung, charge, mus, Eion

def write_csv(strtype='general'):
    if strtype == "general":
        mset = read_json(filename = "data.json")
        f = open('general.csv', 'w')
        print >> f, "(E_DFT, mean.mass, mean.elecneg, mean.radius, mean.ionization"
        for atoms in mset:
            elecneg = []
            rad = []
            Eionization = []
            for name in atoms.names:
                elecneg.append(pauling[name])
                rad.append(radius[name])
                Eionization.append(Eion[name])
            sum_elecneg = np.sum(elecneg) / len(elecneg)
            std_elecneg = np.std(elecneg) 
            sum_mass = np.sum(atoms.masses) / len(atoms.masses)
            std_mass = np.std(atoms.masses)
            sum_radius = np.sum(rad) / len(rad)
            std_radius = np.std(rad)
            sum_Eionization = np.mean(Eionization)
            print >>f, "%6.2f, %6.2f, %6.2f, %6.2f, %6.2f"%(atoms.Eref, sum_mass, 
                                                            sum_elecneg, 
                                                            sum_radius, sum_Eionization)
                                                                                    
    elif strtype == "RS":
        mset = read_json(filename = "data_RS.json")
        f = open('RS.csv', 'w')
        print >> f, "(name0, name1, el, sum.mass, dmass, sum.elecneg, delecneg, calcvol, sum.radius, draius, dpos, Emadelung, Ecoh"
        for atoms in mset:
            elecneg1 = pauling[atoms.names[0]]
            elecneg2 = pauling[atoms.names[1]]
            volscaled = atoms.calcvol**(1./3.) / (radius[atoms.names[0]] * radius[atoms.names[1]]) * 10**4
            delecneg = np.abs(elecneg1-elecneg2)
            sqrtneg = np.std([elecneg1, elecneg2]) #np.sqrt(elecneg1*elecneg2)
    #        Elatt = Emadelung["%s"%(atoms.icsdno)]
            for el in atoms.names:
                if el in charge.keys():
                    Elatt = Emadelung["%s"%(atoms.icsdno)] / charge[el]**2
                    break
            
            dmass = np.abs(atoms.masses[0] - atoms.masses[1])
            d = atoms.positions[0] - atoms.positions[1]
            dpos = np.sqrt(np.inner(d, d))
    
            print >>f, "%s, %s, %s, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f"%(atoms.names[0], atoms.names[1], el, atoms.masses[0] + atoms.masses[1], np.abs(atoms.masses[0] - atoms.masses[1]), (elecneg1 + elecneg2), delecneg, atoms.calcvol, radius[atoms.names[0]] + radius[atoms.names[1]], np.abs(radius[atoms.names[0]] - radius[atoms.names[1]]), dpos, Elatt, atoms.Eref)


    elif strtype == "exptgap":
        mset = read_json(filename = "exptgap.json")
        f = open('exptgap.csv', 'w')
        print >> f, "(formula, std.mass, sum.mass, diff.mass, std.elecneg, sum.elecneg, diff.elecneg, std.radius, sum.radius, diff.radius, bandgap"
        for atoms in mset:
            elecneg = []
            rad = []
            for name in atoms.names:
                elecneg.append(pauling[name])
                rad.append(radius[name])
            sum_elecneg = np.sum(elecneg) / len(elecneg)
            std_elecneg = np.std(elecneg) 
            diff_elecneg = np.max(elecneg) - np.min(elecneg)
            sum_mass = np.sum(atoms.masses) / len(atoms.masses)
            std_mass = np.std(atoms.masses)
            diff_mass = np.max(atoms.masses) - np.min(atoms.masses)
            sum_radius = np.sum(rad) / len(rad)
            std_radius = np.std(rad)
            diff_radius = np.max(rad) - np.min(rad)
            print >>f, "%s,  %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f"%(atoms.formula, std_mass, sum_mass, diff_mass, 
                                                                                                    std_elecneg, sum_elecneg, diff_elecneg,
                                                                                                    std_radius, sum_radius, diff_radius, atoms.Eref)

if __name__ == "__main__":

    write_csv("general")
#    write_csv("RS")
#    write_csv("exptgap")
