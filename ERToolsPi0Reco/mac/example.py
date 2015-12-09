import sys
from ROOT import gSystem
gSystem.Load("libpi0_ERToolsPi0Reco")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing ERToolsPi0Reco..."

sys.exit(0)

