import sys
from ROOT import gSystem
gSystem.Load("libpi0_ERToolsPi0Train")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing ERToolsPi0Train..."

sys.exit(0)

