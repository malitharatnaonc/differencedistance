# MALITHA RATNAWEERA
# VERSION 1.0 2021-08-09

# imports
# requires base env and following packages: pandas, numpy
import pandas as pd 
import argparse
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import os

"""------------------------------------------------------------------------------------------------"""

# Definitions
aminos = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E",
    "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I",
    "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
    "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "SM"
}


"""------------------------------------------------------------------------------------------------"""

## Important calculation functions ##
"""-------------------------------"""


def str2bool(v):
    """Conversion from an input string to a boolean
    Adapted from: https://github.com/symonsoft/str2bool
    """

    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def diff(x):
    """Generates all possible differences between values in a list.
    
    Parameters
    ----------
    x : list
        List of numbers

    Returns
    -------
    diff_x : list
        List of differences between all values in x 
    
    Examples
    --------
    >>> diff([1, 2, 3])
    [0, -1, -2, 1, 0, -1, 2, 1, 0]

    """

    
    diff_x = []
    for a1 in x:
        for b1 in x:
            diff_x.append(a1 - b1)

    return diff_x



def absolute(a, b, c):
    """Uses the list of differences to construct an absolute difference"""
    xdiff = diff(a)
    ydiff = diff(b)
    zdiff = diff(c)

    absolute = []
    for l in range(len(xdiff)):
        s = (xdiff[l] ** 2) + (ydiff[l] ** 2) + (zdiff[l] ** 2)
        root = math.sqrt(s)
        absolute.append(root)

    return absolute



def checker(df):
    """
    In a filetype, searched for any duplicate residues and keeps
    one with higher occupancy or lower bfactor

    Parameters
    ----------
    df : pd.DataFrame object
        Input dataframe

    Returns
    -------
    dfout : pd.DataFrame object
        Dataframe filtered on best residue rotamer (generally should not matter for CA)
    """
    
    newdf = df.sort_values(['Occupancy', 'Temperature Factor'], ascending=(False, True)).drop_duplicates(subset = ['Residue Number', 'Chain Identifier']).sort_index()
    reje = df.loc[~df.index.isin(newdf.index.tolist())]
    
    print("Duplicates exist for following (residues deleted shown):")
    print(" - " + "AA-Res" + "\t" + "Indicator" + "\t" + "Chain")
    print("-" * 40)
    log = [" - " + str(c) + "-" + str(b) + "\t" + str(a) + "\t" + str(d) for a, b, c, d in zip(reje["Location Indicator"].tolist(),
                                                    reje["Residue Number"].tolist(),
                                                    reje["Residue Name"].tolist(),
                                                    reje["Chain Identifier"].tolist())]
    print("\n".join(log))
    
    return newdf




def residue_corrector(af, bf, aChains=None, bChains=None):
    """For the two peptides, residue numbers are compared to ensure they match.
    
    Parameters
    ----------
    af, bf : pd.DataFrame
        Columns of the residue numbers for two proteins
    aChains, bChains : list
        List containing chains to compare


    Returns
    -------
    f1 : pd.DataFrame
        filtered a dataframe
    f2 : pd.DataFrame
        filterd b dataframe

    """

    a, b = af.filt, bf.filt
    assert len(aChains) == len(bChains), "Lists of chains provided not of equal length."

    # checking comparative chains
    f1 = pd.DataFrame([], columns = a.columns)
    f2 = pd.DataFrame([], columns = a.columns)
    for aC, bC in zip(aChains, bChains):
        a1, b1 = a[a["Chain Identifier"] == aC], b[b["Chain Identifier"] == bC]
        for ind, row in a1.iterrows():
            if row["Residue Number"] in b1["Residue Number"].tolist():
                if b1[b1["Residue Number"] == row["Residue Number"]]["Residue Name"].values  == row["Residue Name"]:
                    f1 = f1.append(row)
                    f2 = f2.append(b1[b1["Residue Number"] == row["Residue Number"]])

    info_a = []
    for x, y, z, zb in zip(f1["Residue Number"].tolist(), f1["Chain Identifier"].tolist(), f1["Residue Name"].tolist(), f2["Residue Name"].tolist()):
        if z == zb:
            if z in aminos.keys():
                aa = aminos[z]
            else:
                aa = "X"
            info_a.append(str(y)+"-"+str(x)+"-"+str(aa))
        else:
            info_a.append(str(y)+"-"+str(x)+"-X")
    
    noncomm1 = [str(y) + "\t" + str(z) + "-" + str(x) for x, y, z in zip(a["Residue Number"].tolist(), f1["Chain Identifier"].tolist(), f1["Residue Name"]) if not (x in f1[f1["Chain Identifier"].values == y]["Residue Number"].tolist())]
    noncomm2 = [str(y) + "\t" + str(z) + "-" + str(x) for x, y, z in zip(b["Residue Number"].tolist(), f2["Chain Identifier"].tolist(), f1["Residue Name"]) if not (x in f2[f2["Chain Identifier"].values == y]["Residue Number"].tolist())]

    print("Following residue(s) omitted from first pdb {0}".format(af.pdbfile))
    print("\n".join(noncomm1))
    
    print("\n")
    
    print("Following residue(s) omitted from second pdb {0}".format(bf.pdbfile))
    print("\n".join(noncomm2))
    
    return f1, f2, info_a


"""------------------------------------------------------------------------------------------------"""

## Creating a class for individual DDM files ##
"""-----------------------------------------"""


class DDM():

    def __init__(s, pdbfilename):

        # import pdb file and store as a pandas dataframe
        s.pdbfile = str(pdbfilename)
        if s.pdbfile.split('.')[-1] == 'pdb':
            try:
                with open(s.pdbfile) as file:
                    filelines = [f for f in file.readlines() if f.startswith("ATOM")]
            except:
                raise ValueError("Cannot find the file {0}".format(s.pdbfile))
        else:
            raise ValueError('Please use the format: pdbID.pdb')

        columns = ["Name", "Atom Name", "Location Indicator",
                "Residue Name", "Chain Identifier", "Residue Number",
                "Code for Insertion of Residues", "X", "Y", "Z", "Occupancy",
                "Temperature Factor", "Element Symbol", "Charge on Atom"]
    
        results = {}
        if len(filelines[-1]) == 81:
            for line in filelines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    results[int(line[6:11])] = (line[0:6].replace(" ",""), line[12:16].replace(" ",""), line[16], 
                                        line[17:21].replace(" ",""), line[21], int(line[22:26]),
                                        line[26], float(line[30:38]), float(line[38:46]),
                                        float(line[46:54]), float(line[54:60]), float(line[60:66]),
                                        line[76:78], line[78:80])

            # collect pdb data into a pd Dataframe
            s.data = pd.DataFrame.from_dict(results, columns=columns, orient='index')
            s.data.index.name = "Atom Number"

        elif len(filelines[-1]) != 81:
            raise ValueError("The pdb file you've entered isn't complete.")


    def __sub__(s, other):

        print("\n\n\n")
        print("Starting delta difference matrix calculations.")
        print("-" * 60)
        print("\nComparing dataframes to ensure residues are matched.\n")
        dfs = residue_corrector(s, other, aChains = s.chains, bChains = other.chains)

        first, second, info = dfs[0], dfs[1], dfs[2]
    
        ag_x = np.array(first['X'])
        ag_y = np.array(first['Y'])
        ag_z = np.array(first['Z'])
    
        bg_x = np.array(second['X'])
        bg_y = np.array(second['Y'])
        bg_z = np.array(second['Z'])

        ag = absolute(ag_x, ag_y, ag_z)
        bg = absolute(bg_x, bg_y, bg_z)

        ag = np.reshape(ag, (len(ag_x), len(ag_x)))
        bg = np.reshape(bg, (len(bg_x), len(bg_x)))
    
        delta = ag - bg
        mat = np.matrix(delta)
        ddm = pd.DataFrame(data = mat, columns = info, index = info)
        a_dist = pd.DataFrame(data = ag, columns = info, index = info)
        b_dist = pd.DataFrame(data = bg, columns = info, index = info)


        pathWind = os.getcwd()
        if s.outputs:
            a_dist.to_csv(os.path.join(pathWind, s.pdbfile.split(".")[0] + "_distances.csv"), sep=",", float_format="%.5f")
            b_dist.to_csv(os.path.join(pathWind, other.pdbfile.split(".")[0] + "_distances.csv"), sep=",", float_format="%.5f")
        
        filen = other.pdbfile.split(".")[0] + "--" +  s.pdbfile.split(".")[0] + ".csv"
        
        ddm.to_csv(os.path.join(pathWind, filen), sep=",", float_format="%.5f")
        print("\n" + "-" * 100)
        print("Delta difference matrix of conformation {0} >> {1} saved as '{2}'.".format(other.pdbfile, s.pdbfile, filen))
        print("-" * 100)


        # the difference distance matrix alongside a key reference is outputted as a file.

        return delta



    def pdb_filter(s, chains, atomtype="CA", res_range=None, output=True):
        

        # chain should be inputted as a list to begin with
        s.chains, s.atomtype, s.resrange = list(chains.split(",")), atomtype, res_range
        s.outputs = output
        f = (len(s.pdbfile) + 16)
        b = f - 1
        print("\n\n\n")
        print("-" * (f * 3))
        print("|" + " " * b + "Correcting file {0}".format(s.pdbfile) + " " * b + "|")
        print("-" * (f * 3))

        
        if s.atomtype in ["CA", "CO", "CB"]:
            s.atomtype = s.atomtype
        else:
            raise ValueError("Please provide an atomtype of the following: [CA, CO, CB].")

        data = s.data[s.data["Atom Name"] == atomtype]
        if s.resrange:
            print("Filtering pdb file ({0}) with chains {1}, atom type '{2}', and residue range between {3}-{4}.".format(
                                                                            s.pdbfile, s.chains, s.atomtype,
                                                                            s.resrange[0], s.resrange[1]))
            if type(s.resrange) == tuple and len(s.resrange) == 2:
                s.filt = data[np.logical_and(data["Chain Identifier"].isin(s.chains),
                                        data["Residue Number"] >= s.resrange[0], 
                                        data["Residue Number"] <= s.resrange[1])]
            else:
                raise ValueError("Residue range specified should be in the format XXX,YYY.")
        else:
            s.filt = data[data["Chain Identifier"].isin(s.chains)]
            print("\n\nFiltering pdb file ({0}) with chains {1} and atom type '{2}'.".format(s.pdbfile, 
                                                                                    s.chains, 
                                                                                    s.atomtype))
        
        # correcting any duplicate residue events for multiple occurrences
        print("\n\nChecking for duplicate residues...\n" + "-" * 60 + "\n")
        s.filt = checker(s.filt)

        # output the filtered dataset as a file
        if output:
            outname = s.pdbfile.split(".")[0] + "_" + "-".join(s.chains) + ".csv"
            s.filt.to_csv(outname)

        return s



"""------------------------------------------------------------------------------------------------"""

def run(args):
    """Runs the sript to produce a difference distance matrix.

    """
    if args.a1.split(".")[-1] != "pdb" and args.b1.split(".")[-1] != "pdb":
        print("File types are not both of the format 'ABC.pdb'.")
    else:
        a1 = DDM(args.a1).pdb_filter(chains=args.aC, output=args.nv)
        b1 = DDM(args.b1).pdb_filter(chains=args.bC, output=args.nv)
        a1 - b1




def main():
    """Function to collect user inputs."""

    parser = argparse.ArgumentParser(description="Creates data for the difference distance matrix.")
    parser.add_argument("-a", help="PDB file a (the final conformation)", dest="a1", type=str, required=True)
    parser.add_argument("-b", help="PDB file b (the starting conformation)", dest="b1", type=str, required=True)
    parser.add_argument("-aC", help="Chains to use for pdb file a. (separated by a comma)", dest="aC", type=str, default=None)
    parser.add_argument("-bC", help="Chains to use for pdb file b. (separated by a comma)", dest="bC", type=str, default=None)
    parser.add_argument("--verbose", help="Outputs all intermediate information files.", 
                        dest="nv", type=str2bool, nargs='?', const=True, default=False)
    #parser.add_argument("-t", help="Title of the plot. Start with 'CONTXT:' if you wish to add the change.",
    #                    dest="title", type=str, default='None')
    #parser.add_argument("-o", help="Output Name", dest="output", type=str, default='None')
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()




