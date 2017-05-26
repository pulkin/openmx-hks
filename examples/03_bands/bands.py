import numpy
from matplotlib import pyplot, gridspec
import os, subprocess, json

import sys
sys.path.append("../common")
import tb

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
vectors = list(tuple(i) for i in data["vectors"])
H = numpy.array(data["H"]["r"]) + 1j*numpy.array(data["H"]["i"])
S = numpy.array(data["S"]["r"]) + 1j*numpy.array(data["S"]["i"])

# Subtract Fermi level and convert to eV
H -= data["fermi"]*S
H *= 27.2113845

pprint("Converting to tight binding objects")
H = tb.TightBinding(H, vectors = vectors).squeezed()
S = tb.TightBinding(S, vectors = vectors).squeezed()

pprint("Calculating band structure")
k = numpy.linspace(-.5,.5,50)
values = H.eig_path(k[:,numpy.newaxis], b = S)

pprint("Plotting band structure")
gs = gridspec.GridSpec(1, 4)
pyplot.subplot(gs[0, :3])
for i in range(values.shape[1]):
    pyplot.plot(k,values[:,i], color = "black")
pyplot.xlim(-.5,.5)
pyplot.xlabel("k")
pyplot.ylabel("Energy (eV)")
pyplot.title("Band structure")

pprint("Calculating channels count")
e = numpy.linspace(-10,10,50)
e = e + (e[1]-e[0])*0.01j
device_h = H.periodic_device()
device_s = S.periodic_device()
transmission = []
for energy in e:
    transmission.append(tb.MTDCalculator(device_s*energy - device_h).transmission(0,1).real)

pprint("Plotting channels count")
pyplot.subplot(gs[0, 3], sharey = pyplot.gca())
pyplot.plot(transmission,e.real)
pyplot.ylim(-10,10)
pyplot.xlabel("T")
pyplot.xticks([0,1,2,3])
pyplot.grid()
pyplot.title("Channels")

pyplot.suptitle("Carbon chain, compare to http://www.openmx-square.org/openmx_man3.8/node114.html")

pyplot.show()
