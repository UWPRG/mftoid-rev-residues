import os
import re

path = os.environ.get("GMXLIB")
rdf = os.listdir(path)
charmms = []
latest_version = 0
latest_year = 0
latest_month = 0
latest_dir = ""
months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
for direc in rdf:
    if "charmm36-feb2021" in direc:
        latest_dir = os.path.join(path, direc)
        break
    else:
        match1 = re.match(r"charmm(\d+)-([a-z]+)(\d+).ff", direc)
        match2 = re.match(r"charmm(\d+)", direc)
        if match1:
            version = match1.group(1)
            month = match1.group(2)
            year = match1.group(3)
            if int(version) > latest_version:
                latest_version = int(version)
                latest_year = int(year)
                latest_month = months.index(month)
                latest_dir = os.path.join(path, direc)
            if int(version) == latest_version and latest_year < int(year):
                latest_version = int(version)
                latest_year = int(year)
                latest_month = months.index(month)
                latest_dir = os.path.join(path, direc)
            if int(version) == latest_version and latest_year == int(year) and months.index(month) > latest_month:
                latest_version = int(version)
                latest_year = int(year)
                latest_month = months.index(month)
                latest_dir = os.path.join(path, direc)
        if match2:
            version = match2.group(1)
    #         print("Version", version)
            if int(version) > latest_version:
                latest_version = int(version)
                latest_year = 0
                latest_month = 0
                latest_dir = os.path.join(path, direc)

cd = os.getcwd()
with open(os.path.join(cd, "merged.rtp"), "r") as f1:
    with open(os.path.join(latest_dir, "merged.rtp"), "a+") as f2:
        f2.write(f1.read())
with open(os.path.join(cd, "ffnonbonded.itp"), "r") as f1:
    with open(os.path.join(latest_dir, "ffnonbonded.itp"), "a+") as f2:
        f2.write(f1.read())
with open(os.path.join(cd, "ffbonded/bondtypes.itp"), "r") as f1:
    with open(os.path.join(latest_dir, "ffbonded.itp"), "a+") as f2:
        f2.write(f1.read())
with open(os.path.join(cd, "ffbonded/angletypes.itp"), "r") as f1:
    with open(os.path.join(latest_dir, "ffbonded.itp"), "a+") as f2:
        f2.write(f1.read())
with open(os.path.join(cd, "ffbonded/dihedraltypes.itp"), "r") as f1:
    with open(os.path.join(latest_dir, "ffbonded.itp"), "a+") as f2:
        f2.write(f1.read())
with open(os.path.join(cd, "ffbonded/improper.itp"), "r") as f1:
    with open(os.path.join(latest_dir, "ffbonded.itp"), "a+") as f2:
        f2.write(f1.read())
with open(os.path.join(cd, "atomtypes.atp"), "r") as f1:
    lines = f1.readlines()  # Read all lines into a list
    with open(os.path.join(latest_dir, "atomtypes.atp"), "a+") as f2:
        # Write each line except the last one with a newline
        f2.writelines(lines[:-1]) 
        # Write the last line without adding an extra newline at the end
        f2.write(lines[-1].rstrip('\n'))

