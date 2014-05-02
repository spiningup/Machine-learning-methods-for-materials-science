from pylab import *
from matplotlib.ticker import FuncFormatter
from read_json import attribute_tolist, get_unique_elements

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'


def plot_Eref(mset):
    Eref = attribute_tolist(mset, attr="Eref")
    plot_all(Eref, "Total energy per atom (eV)")


def plot_all(d, xl):
    hist(d, 50, normed=False)
    xlabel(xl, fontsize=14)
#    formatter = FuncFormatter(to_percent)
#    gca().yaxis.set_major_formatter(formatter)
    show()

def plot_natoms(mset):
    natoms = []
    for atom in mset:
        natoms.append(atom.natoms)
    plot_all(natoms, "Number of atoms in simulated cell")

def plot_error_in_volume(mset):
    volerror = []
    for atoms in mset:
        if atoms.exptvol is not None:
            volerror.append((atoms.calcvol - atoms.exptvol) / atoms.exptvol * 100)
            if np.abs(volerror[-1]) > 50:
                print atoms.formula, atoms.icsdno, atoms.exptvol, atoms.calcvol, volerror[-1]

    plot_all(volerror,"Error in volume (%) : (Vcalc - Vexpt) / Vexpt * 100")


def plot_elements(mset):
    a = get_unique_elements(mset)
    X = np.arange(len(a))
    bar(X, a.values(), align="center", width=0.5)
    xticks(X, a.keys())
    ymax = max(a.values()) + 20
    ylim(0, ymax)
    show()


if __name__ == "__main__":
    from read_json import read_json
    mset = read_json("data.json")
    plot_error_in_volume(mset)
    plot_Eref(mset)
    plot_natoms(mset)
    plot_elements(mset)
