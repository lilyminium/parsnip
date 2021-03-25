import yaml


import os
from rdkit import Chem
import MDAnalysis as mda
from MDAnalysis.exceptions import SelectionError

def write_universes_to_file(universes=[], suffix="pdb",
                            name="psiresp", dirname="rdmols"):
    try:
        os.makedirs(dirname)
    except FileExistsError:
        pass
    names = []
    for i, u in enumerate(universes, 1):
        molname = f"{name}_m{i:03d}"
        names.append(molname)
        filename = os.path.join(dirname, f"{molname}.{suffix}")
        u.atoms.write(filename)
    return names

def get_cap_constraints(universes=[]):
    chrequivs = []
    chrconstrs = []
    for i, u in enumerate(universes):
        try:
            ag = u.select_atoms("altloc + -")
        except SelectionError:
            continue
        else:
            # all cap atoms should sum to 0
            constr = []
            for ix in ag.indices:
                constr.append([i, int(ix)])
            chrconstrs.append([0, constr])

        # cap atoms with + should be equiv
        ag = u.select_atoms("altloc +")
        for at in ag:
            ix = int(at.index)
            other = [int(at.resid), int(at.id) - 1]
            chrequivs.append([[i, ix], other])
    return chrequivs, chrconstrs
        

def make_psiresp_job(rdmols=[], suffix="pdb", name="psiresp",
                     dirname="rdmols", jobfile="job.yml",
                     n_confs=50, minimize=True, basis="6-31g*",
                     method="scf", opt=True,
                     force=False, verbose=True):
    charges = [Chem.GetFormalCharge(mol) for mol in rdmols]
    universes = [mda.Universe(mol, format="RDKIT") for mol in rdmols]
    names = write_universes_to_file(universes=universes, suffix=suffix,
                                    name=name, dirname=dirname)
    chrequivs, chrconstrs = get_cap_constraints(universes)

    moldct = {}
    for molname, charge in zip(names, charges):
        moldct[molname] = {"charge": charge}

    dct = dict(
        n_confs=n_confs,
        minimize=minimize,
        basis=basis,
        method=method,
        opt=opt,
        force=force,
        verbose=verbose,
        coords=None,
        rdmol=os.path.join(dirname, f"{{name}}.{suffix}"),
        out=f"{{name}}_charged.pdbqt",
        molecules=moldct,
    )

    chrequiv_yml = yaml.dump({'inter_chrequiv': chrequivs},
                              default_flow_style=True)[1:-2]
    chrconstr_yml = yaml.dump({'inter_chrconstr': chrconstrs},
                               default_flow_style=True)[1:-2]
    other = yaml.dump(dct, default_flow_style=False)

    with open(jobfile, "w") as f:
        f.write("\n".join([other, chrequiv_yml, chrconstr_yml]))



def write_polymer_script(unit_indices, universes, names, suffix="pdb"):
    monomers = {}

    for uname, u in zip(names, universes):
        mon = {"file": f"{uname}_charged.topcrd",}
    for uix in unit_indices:
        uname = names[uix]
        