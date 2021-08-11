# Generate Difference Distance Matrix (DDM)

Use this script to generate a difference distance matrix file from two structure files (.pdb) of the same protein. It does not yet create the plots but instead outputs the values for the matrix alongside residue numbers and chain information for easy manipulation in R.

Usage:
```
ddm.py -a [end pdb file] -b [start pdb file] -aC [chains for end pdb] -bC [chains for start pdb] [options]
```

## Getting Started
* Make sure you have Python 3 installed and have the **pandas** and **numpy** packages installed.
* Ensure your starting conformation and end conformation pdb files are in the current directory.

## Required Fields
* ```-a```	PDB file a (the final conformation)
* ```-b```	PDB file b (the starting conformation)
* ```-aC```	Chains to use for pdb file a. (separated by a comma)
* ```-bC```	Chains to use for pdb file b. (separated by a comma) [should be equivalent to those of pdb file a]

## Optional Fields
* ```--verbose```	Choice of outputting all intermediate information files in csv format.

## Running Example
Say we wish to obtain the ddm for the conformational shifts from 7af1.pdb to 7apv.pdb, only looking at chain A for both, we would then run:
```
ddm.py -a 7apv.pdb -aC A -b 7af1.pdb -bC A
```

If instead we realised that chain A of 7af1.pdb is equivalent to chain B of 7apv.pdb:
```
ddm.py -a 7apv.pdb -aC B -b 7af1.pdb -bC A
```

To overlay multiple equivalent chains (A=A, B=B in this case):
```
ddm.py -a 7apv.pdb -aC A,B -b 7af1.pdb -bC A,B
```

## Functionality
**COMPLETE FIELD**

## Updates
* 2021-03-29 - Version 1.0

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.