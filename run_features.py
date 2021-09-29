import os
import sys

import rdkit
from rdkit import Chem
import get_pdbinfo



def renumber(fmt, infile, outfile):
    """
    Rename atoms in file based on order
    """
    if fmt == "mol2":
        lines = [line for line in open(infile)]
        atom_index = lines.index("@<TRIPOS>ATOM\n")
        bond_index = lines.index("@<TRIPOS>BOND\n")
        atom_lines = lines[atom_index + 1: bond_index]
        atoms = set([atom.split()[5].split(".")[0] for atom in atom_lines])
        atom_dic = {key: 1 for key in atoms}
        atom_lines_new = []
        for line in atom_lines:
            atom_key = line.split()[5].split(".")[0]
            atom_old = line.split()[1]
            atom_new = atom_key.upper() + str(atom_dic[atom_key])
            if len(atom_old) > len(atom_new):
                atom_new = atom_new + (len(atom_old) - len(atom_new)) * " "
            elif len(atom_old) < len(atom_new):
                atom_old = atom_old + (len(atom_new) - len(atom_old)) * " "

            newline = line.replace(atom_old, atom_new, 1)
            atom_lines_new.append(newline)
            atom_dic[atom_key] += 1

        outfile = open(outfile, "w")
        outfile.write("".join(lines[0:atom_index + 1]))
        outfile.write("".join(atom_lines_new))
        outfile.write("".join(lines[bond_index:]))
        outfile.close()
    elif fmt == "pdb":
        lines = [line for line in open(infile)]
        atom_lines = get_pdbinfo.pdbinfo(file=infile).getAtoms()
        atom_index = lines.index(atom_lines[0])
        bond_index = lines.index(atom_lines[-1]) + 1
        mol = Chem.MolFromPDBFile(infile, removeHs=False)
        atom_list = [atom.GetSymbol() for atom in mol.GetAtoms()]
        atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
        atom_dic = {key: 1 for key in atoms}
        atom_lines_new = []
        for idx, line in enumerate(atom_lines):
            atom_key = atom_list[idx]
            atom_old = get_pdbinfo.atmn(line).strip()
            atom_new = atom_key.upper() + str(atom_dic[atom_key])
            if len(atom_old) > len(atom_new):
                atom_new = atom_new + (len(atom_old) - len(atom_new)) * " "
            elif len(atom_old) < len(atom_new):
                atom_old = atom_old + (len(atom_new) - len(atom_old)) * " "
            newline = line.replace(atom_old, atom_new, 1)
            atom_lines_new.append(newline)
            atom_dic[atom_key] += 1

        outfile = open(outfile, "w")
        outfile.write("".join(lines[0:atom_index]))
        outfile.write("".join(atom_lines_new))
        outfile.write("".join(lines[bond_index:]))
        outfile.close()


def get_input(datadir, fn):
    """
    Get input files based on pdbid
    return: inlig_rdkit --> inlig mol2 or sdf for RDkit using
            inlig3 --> inlig pdb file
            inpro1 --> inpro with only protein
            inpro2 --> inpro with both proteina and water
    """
    olddir = os.getcwd()
    os.chdir(datadir)
    ### get ligand ###
    ### input should be provided as either of sdf file (best choice) or mol2 file ###
    inlig1 = fn + "_ligand.mol2"
    inlig2 = fn + "_ligand.sdf"
    inlig3 = fn + "_ligand.pdb"
    ### check ligand input file ###
    inlig_rdkit = None
    if inlig1 in os.listdir("."):
        inlig = inlig1
        try:
            mol = Chem.MolFromMol2File(inlig, removeHs=False)
            if mol == None:
                if inlig2 in os.listdir("."):
                    inlig = inlig2
                    mol = Chem.SDMolSupplier(inlig, removeHs=False)[0]
                    if mol != None:
                        inlig_rdkit = inlig2
            else:
                inlig_rdkit = inlig1
        except:
            pass
    else:
        inlig = inlig2
        try:
            mol = Chem.SDMolSupplier(inlig, removeHs=False)[0]
            if mol != None:
                inlig_rdkit = inlig2
        except:
            pass
    ### correcting atom name for sasa calculation ###
    ### if inlig sdf or mol2 file might have problems, you should also provide a pdb file that can be used to conduct other calculation ###
    try:
        if inlig3 not in os.listdir(".") and inlig_rdkit == None:
            ### if sdf and mol2 files can't be processed by RDKit, we need to generate pdb file and continue other calculations ###
            ### this is not good, since if that molecule has large problem, the conversion process might be wrong; but this can work when the problem is caused by RDKit ###
            print("sdf and mol2 can't be processed")
            infmt = inlig.split(".")[-1]
            outlig = inlig3
            cmd = "obabel -i" + infmt + " " + outlig + " -opdb -O " + outlig    #CHANGED
            # obabel -i sdf 1vso_ligand.pdb -opdb -O 1vso_ligand.pdb
            os.system(cmd)
    except:
        print(f"<<<<<<<<{inlig} failed!")
        raise RuntimeError()

    if inlig3 not in os.listdir("."):
        inlig = inlig_rdkit
        infmt = inlig.split(".")[-1]
        if infmt == "mol2":
            outlig_num = inlig.split(".")[0] + "_rename.mol2"
            renumber(infmt, inlig, outlig_num)
            outlig = inlig3.split(".")[0] + "_rename.pdb"
            cmd = "obabel -i" + infmt + " " + outlig_num + " -opdb -O " + outlig
            os.system(cmd)
        elif infmt == "sdf":
            outlig = inlig3
            cmd = "obabel -i" + infmt + " " + inlig + " -opdb -O " + outlig
            os.system(cmd)
            outlig_num = outlig.split(".")[0] + "_rename.pdb"
            renumber('pdb', outlig, outlig_num)
    else:
        infmt = inlig3.split(".")[-1]
        outlig_num = inlig.split(".")[0] + "_rename.pdb"
        renumber(infmt, inlig3, outlig_num)

    inlig3 = fn + "_ligand_rename.pdb"

    ### get protein ###
    ### at least one protein structure should be provided with all waters ###
    inpro1 = fn + "_protein.pdb"
    inpro2 = fn + "_protein_all.pdb"    #CHANGED
    os.system(f"cp {inpro1} {inpro2}")  #CHANGED
    if inpro1 not in os.listdir("."):
        inpro = inpro2
        outpro = open(inpro1, "w")
        protein_lines = get_pdbinfo.pdbinfo(fn, file=inpro).getProteinWaters()[0]
        outpro.write("".join(protein_lines))
        outpro.close()
    ### check input structures ###
    if inlig_rdkit != None and os.path.isfile(inlig_rdkit) and os.stat(inlig_rdkit).st_size != 0:
        print("Ligand for conformation stability:" + inlig_rdkit)
    else:
        print(
            "Warning:input ligand should be checked, skip ligand stability calculation, use default(dE:-300, RMSD:300)")

    if os.path.isfile(inlig3) and os.stat(inlig3).st_size != 0:
        print("Ligand for Vina, SASA, BA, ION:" + inlig3)
    else:
        sys.exit("Error: ligand input (pdb)")
    if os.path.isfile(inpro1) and os.stat(inpro1).st_size != 0:
        print("Protein without water molecules:" + inpro1)
    else:
        sys.exit("Error: protein input without water")
    if os.path.isfile(inpro2) and os.stat(inpro2).st_size != 0:
        print("Protein with water molecules:" + inpro2)
    else:
        sys.exit("Error: protein input with water")
    os.chdir(olddir)

    return inlig_rdkit, inlig3, inpro1, inpro2