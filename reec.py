# Receptor Energetic Efficacy Calculator

import sys
import re
import os
import os.path
import json
import subprocess
import pprint

if len(sys.argv) < 3:
    print("Both a protein and a ligand are required.")
    exit

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.getcwd())
os.chdir("data")

import data.globals
import data.protutils
import data.odorutils
import data.dyncenter

data.protutils.load_prots()
data.odorutils.load_odors()

protid = sys.argv[1]
if not protid in data.protutils.prots:
    print("Bad receptor ID " + protid)
    exit()
ligid = sys.argv[2]
active = True
if (len(sys.argv) > 3):
    c = sys.argv[3][0].upper()
    active = (c == 'A')

if ligid in data.odorutils.odors:
    odor = data.odorutils.odors[ligid]
else:
    odor = data.odorutils.find_odorant(ligid)

if not odor:
    print("Bad odorant " + protid)
    exit()

fam = data.protutils.family_from_protid(protid)
uname = odor["full_name"].replace(' ', '_')
acvi = "active"
if not active: acvi = f"in{acvi}"
dockname = f"out/{fam}/{protid}/{protid}~{uname}.{acvi}.dock"
pdbname = f"out/{fam}/{protid}/{protid}~{uname}.{acvi}.model1.pdb"

cmd = ["test/occlusion_test", dockname]
print(" ".join(cmd), "\n")
proc = subprocess.run(cmd, stdout=subprocess.PIPE)
for ln in proc.stdout.decode().split('\n'):
    print(ln)

# TODO: Find mutual nearest atoms between ligand and Y/F6.48.
# This residue is a toggle switch for Class II ORs and the docks of OR8B8,
# a putative mintiness receptor, have l-carvone right up to the aromatic ring,
# while d-carvone lies some distance away. l-Carvone is also a better fit for
# V108(3.36) and V203(5.43), so if you can quantify these interactions then
# that may be another basis for predicting agonism.
