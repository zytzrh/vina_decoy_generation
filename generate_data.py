# -----------------------------------------------------------------------------
# Optimized Structure
# -----------------------------------------------------------------------------
import os
import fileinput
import sys
import get_pdbinfo
from run_features import get_input
import env_path
import glob
from tqdm import tqdm


#### Do minimization, structure file should have hydrogen atom!!!####
## generate the refined data list
## generate the box for mol2 file


def get_box(fn, inlig):  # SOS
    if inlig.split(".")[-1] == "mol2":
        inputfile = open("../" + inlig)
        x = []
        y = []
        z = []
        flag = False
        for line in inputfile:
            if line[0:13] == "@<TRIPOS>ATOM":
                flag = True
                continue
            elif line[0:13] == "@<TRIPOS>BOND":
                flag = False
            if flag:
                x.append(float(line[16:26].split()[0]))
                y.append(float(line[27:37].split()[0]))
                if line[37] == "0":
                    z.append(float(line[38:48].split()[0]))
                else:
                    z.append(float(line[37:48].split()[0]))
    elif inlig.split(".")[-1] == "pdb":
        lines = get_pdbinfo.pdbinfo(name=fn, file="../" + inlig).getAtoms()
        x, y, z = [], [], []
        for line in lines:
            coords = get_pdbinfo.pdbinfo(name=fn, lines=[line]).getCoords()
            x.append(float(coords[0][0]))
            y.append(float(coords[0][1]))
            z.append(float(coords[0][2]))

    x_center = (max(x) + min(x)) / 2
    y_center = (max(y) + min(y)) / 2
    z_center = (max(z) + min(z)) / 2
    size_x = max(x) - min(x) + 10
    size_y = max(y) - min(y) + 10
    size_z = max(z) - min(z) + 10
    new_file = open("box.txt", "w")
    new_file.write("center_x = " + str(x_center) + "\n")
    new_file.write("center_y = " + str(y_center) + "\n")
    new_file.write("center_z = " + str(z_center) + "\n")

    new_file.write("size_x = " + str(size_x) + "\n")
    new_file.write("size_y = " + str(size_y) + "\n")
    new_file.write("size_z = " + str(size_z) + "\n")
    new_file.close()


def genpdbqt(fn, ligpdb, propdb):
    propdbqt = fn + "_rec.pdbqt"
    ligpdbqt = fn + "_lig.pdbqt"
    cmd1 = "$MGLPY $MGLUTIL/prepare_receptor4.py -r ../" + propdb + " -o " + propdbqt + " -U '_' > out1.tmp"
    cmd2 = "$MGLPY $MGLUTIL/prepare_ligand4.py -l ../" + ligpdb + " -o " + ligpdbqt + " -U '_' -Z > out2.tmp"
    os.system(cmd1)
    os.system(cmd2)


def runmin(fn):
    cmd = "$VINADIR/vina --receptor " + fn + "_rec.pdbqt --ligand " + fn \
          + "_lig.pdbqt --config box.txt --local_only --out " + fn + "_lig_min.pdbqt >out_min.txt"  # SOS: local_only
    os.system(cmd)


def chanPdb(fn):
    cmd = "obabel -ipdbqt " + fn + "_lig_min.pdbqt  -opdb -O " + fn + "_lig_min.pdb"
    os.system(cmd)


def get_Co(datadir, fn, inlig, st):
    """
    Get Vina optimized structure

    :param datadir:datadir for input
    :param fn: input index
    :param inlig: inlig
    :param st: water type   #TODO

    """
    os.chdir(datadir)
    olddir = os.getcwd()
    os.system("mkdir vinamin_rigid")
    os.chdir("vinamin_rigid")
    if st != "" and "_" not in st:
        st = "_" + st

    inpro = fn + "_protein" + st + ".pdb"
    outlig = fn + "_lig_min" + st + ".pdb"

    genpdbqt(fn, inlig, inpro)
    get_box(fn, inlig)
    runmin(fn)
    chanPdb(fn)
    os.system("cp " + fn + "_lig_min.pdb ../" + outlig)
    os.chdir(olddir)
    # os.system("rm -r vinamin_rigid")

if __name__ == "__main__":
    PDBbind_path = "/NAS2020/Workspaces/DRLGroup/zyt/dataset/PDBbind-2016-benchmark/refined-set-dockrefined"
    target_ids = [f.split('/')[-1] for f in
                glob.glob(os.path.join(PDBbind_path, "*"), recursive=False)]
    for target_id in tqdm(target_ids):
        print(f">>>Process {target_id}")
        try:
            datadir = os.path.join(PDBbind_path, target_id)
            # datadir = "/NAS2020/Workspaces/DRLGroup/zyt/Mila/vina_decoy_generation/example/1e1v"
            # pdbid = "1e1v"
            inlig_rdkit, inlig_pdb, inpro_pro, inpro_water = get_input(datadir, target_id)
            get_Co(datadir=datadir,
                   fn=target_id,
                   inlig=inlig_pdb,
                   st="")
        except:
            continue






