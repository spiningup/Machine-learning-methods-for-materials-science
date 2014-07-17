import numpy as np
import sys
from pylada import vasp
from pylada.crystal import read

def find_band_crossing(E, eig):
    nkpt, nbands = eig.shape
    bandidx = []
    for i in range(nbands):
        if min(eig[:,i]) < E and max(eig[:, i]) > E:
            bandidx.append(i)
    return bandidx

def save_bxsf(filename, E, eig, bandidx, nk, Blatt):
    if len(bandidx) == 0:
        print "No band crossing"
        return None

    f = open(filename, 'w')
    f.write("BEGIN_INFO\n")
    f.write(" Fermi Energy: %f\n" %(E))
    f.write("END_INFO\n\n")

    f.write("BEGIN_BLOCK_BANDGRID_3D\n")
    f.write(" Bugfixed_by_mbB_28.06.2013\n")
    f.write(" BEGIN_BANDGRID_3D\n")
    f.write(" %d \n"%(len(bandidx)))
    f.write(" %d %d %d\n"%(nk[0]+1, nk[1]+1, nk[2]+1))
    f.write(" 0 0 0 \n") # origin gamma
    f.write(" %f %f %f \n" %(Blatt[0,0], Blatt[0,1], Blatt[0,2]))
    f.write(" %f %f %f \n" %(Blatt[1,0], Blatt[1,1], Blatt[1,2]))
    f.write(" %f %f %f \n" %(Blatt[2,0], Blatt[2,1], Blatt[2,2]))

    for i in range(len(bandidx)):
        f.write(" BAND: %d \n" %(bandidx[i])) # band index starts from 0
        for iz in range(nk[2]+1):
            for iy in range(nk[1]+1):
                for ix in range(nk[0]+1):
                    f.write(" %f"%(eig[ix%nk[0], iy%nk[1], iz%nk[2], bandidx[i]]))
                f.write("\n")
            f.write("\n")
    f.write(" END_BANDGRID_3D\n")
    f.write("END_BLOCK_BANDGRID_3D\n")
    f.close()
    return 

def read_nk(dirname):
    item = open(dirname+'/KPOINTS', 'r').readlines()
    assert "Gamma" in item[2]
    nkpt = item[3].split()
    nk = [int(k) for k in nkpt]
    return nk
    
def index_transform(nk0, old):
    if nk0 % 2 == 0:
        if old <= nk0 // 2:
            new = old + nk0 // 2 - 1
        else:
            new = nk0 - 1 - old
    else:
        print "To be implement"
        XX
    return new

def sort_eig(nk, eigold):
    eig = np.zeros((nk[0], nk[1], nk[2], eigold.shape[1])) 
    # x varies the fastest in the eig
    ik = 0
    for iz in range(nk[2]):
        for iy in range(nk[1]):
            for ix in range(nk[0]):
                eig[ix, iy, iz, :] = eigold[ik, :]
                ik += 1
    return eig

def get_Nb(idx, bandtype="valence"):
    Nb = 0

    for iband in idx:
        if bandtype == "valence":
            mineig = np.max(eig[:,:,:,iband])
        elif bandtype == "conduction":
            mineig = np.min(eig[:,:,:,iband])
        else: XX
        Nb += len(np.where(np.abs(eig[:,:,:,iband] - mineig) < 1e-4)[0])

    return Nb
    

def get_neighbors(ix, iy, iz, iband, eig, n_neighbors=5):
    nghs = []
    imin = - n_neighbors // 2
    imax = n_neighbors // 2 + 1
    for i in range(imin, imax):
        for j in range(imin, imax):
            for k in range(imin, imax):
                if i == 0 and j == 0 and k == 0: continue
                nghs.append(eig[(ix+i+nk[0])%nk[0], (iy+j+nk[1])%nk[1], (iz+k+nk[2])%nk[2], iband] )

    return nghs

def get_Nb_by_nghs(idx, nk, eig, band_type="valence", n_neighbors=5):                
    Nb = 0
    for iband in idx:
        for ix in range(nk[0]):
            for iy in range(nk[1]):
                for iz in range(nk[2]):
                    nghs = get_neighbors(ix, iy, iz, iband, eig, n_neighbors)
                    if band_type == "valence":
                        if np.all(eig[ix, iy, iz, iband] >= nghs): Nb += 1
                    elif band_type == "conduction":
                        if np.all(eig[ix, iy, iz, iband] <= nghs): Nb += 1
                    else: XX
    return Nb


if __name__ == "__main__":
    dirname = sys.argv[1]
    calc = vasp.Extract(dirname)
    
    eig = calc.eigenvalues.magnitude
    fermi = calc.fermi_energy.magnitude
    nk = read_nk(dirname)
    
    # get vbm and cbm
    vbm=max([x for x in eig.flatten() if float(x)<=fermi])
    cbm=min([x for x in eig.flatten() if float(x) >fermi])
    print "fermi", fermi, vbm, cbm
    
    # find bands cross vbm/cmb + window
    window = 0.1 # eV
    vbm_idx = find_band_crossing(vbm-window, eig)
    cbm_idx = find_band_crossing(cbm+window, eig)
    print "Number of bands crossing (vbm, cbm)", len(vbm_idx), len(cbm_idx)
    
    # sort eig to eig[nkx, nky, nkz, nbands] with the correct kpoints order
    eig = sort_eig(nk, eig)
    
    # get Nb
    Nb_v = get_Nb(vbm_idx, bandtype="valence")
    Nb_c = get_Nb(cbm_idx, bandtype="conduction")
    print "band degeneracy (vbm, cbm)", Nb_v, Nb_c
    
    Nb2_v = get_Nb_by_nghs(vbm_idx, nk, eig, band_type="valence", n_neighbors=7)
    Nb2_c = get_Nb_by_nghs(cbm_idx, nk, eig, band_type="conduction", n_neighbors=7)
    print "band degeneracy by finding peaks", Nb2_v, Nb2_c
    
    # save to bxsf
    atoms = read.poscar(dirname+"/CONTCAR")
    Acell = atoms.cell
    Bcell = np.linalg.inv(Acell.T)
    
    save_bxsf(dirname+'/vbm.bxsf', vbm-window, eig, vbm_idx, nk, Bcell)
    save_bxsf(dirname+'/cbm.bxsf', cbm+window, eig, cbm_idx, nk, Bcell)
