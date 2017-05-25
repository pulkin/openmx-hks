from matplotlib import pyplot
import os, subprocess, json

fname = "../01_hks2json/lead-chain.json"

def pprint(x):
    print("[python] {}".format(x))
    
pprint("Checking if file {} exists".format(fname))
if not os.path.isfile(fname):
    pprint("File does not exist, running previous example")
    subprocess.Popen(["./run"], cwd = "../01_hks2json").wait()
    
pprint("Opening {}".format(fname))
with open(fname,'r') as f:
    data = json.load(f)

pprint("Retrieving JSON data")
vectors = data["vectors"]
H = data["H"]
S = data["S"]

pprint("Finding diagonal block")
index = vectors.index([0,0,0])
H_r = H["r"][index]
H_i = H["i"][index]
S_r = S["r"][index]
S_i = S["i"][index]

pprint("Plotting")
for x, (m, name) in enumerate((
    (H_r,"Hamiltonian real"),
    (H_i,"Hamiltonian imaginary"),
    (S_r,"Overlap real"),
    (S_i,"Overlap imaginary"),
)):
    pyplot.subplot(221+x)
    pyplot.imshow(m, interpolation = 'none')
    pyplot.title(name)
    pyplot.colorbar()

pyplot.show()
