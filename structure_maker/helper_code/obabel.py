import os
import subprocess
directory = "minima_pdb/"
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    strf = str(f)
    command = "obabel " + strf + " -h -O " + strf[:-4] + "_h.pdb" 
    subprocess.run(command.split(), stdout=subprocess.PIPE)
