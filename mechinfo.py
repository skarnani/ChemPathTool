#
# Functions to supply reaction mechanism data for use by ChemPathTool.
# These functions require the Cantera package, which is available from
# www.cantera.org.
#
# D. G. Goodwin, 8/7/02
#

#Commented CanteraError out because it doesn't appear in new versions
#from cantera import CanteraError 
from cantera import Solution
#Changed from cantera.gases import IdealGasMix

# Store the object representing the mechanism here
class __data:
    g = None
    nsp = 0

def readSolution(solution):
    """Read in a solution file."""
    g = solution
    __data.g = g
    __data.nsp = g.n_species
    
def readMechanism(infile, thermo=""):
    """Read in a mechanism file."""
    g = Solution(infile, thermo)
    __data.g = g
    __data.nsp = g.n_species

def setState(T, P, X):
    """Set the gas temperature [K], pressure [Pa],
    and mole fractions."""
    __data.g.TPX = T, P, X

def elementNames():
    return __data.g.element_names

def indexElt(elt):
    return __data.g.element_index(elt)

def indexSpec(elt):
    return __data.g.species_index(elt)

def elementAtomicWt():
    return __data.g.atomic_weights

def speciesNames():
    return __data.g.species_names

def numberOfElementXinSpeciesY(X,Y):
    _m = __data.g.element_index(X)
    _k = __data.g.species_index(Y)
    return __data.g.n_atoms(_k, _m)

def numSpecies():
    return __data.g.n_species

def numReactions():
    return __data.g.n_reactions

def specCoeffsInReaction(r):
    """Return a list of pairs of (species name, coefficient) for reaction
    number r. Only those species with non-zero stoichiometric coefficient
    are included."""
    c = []
    for k in range(__data.nsp):
        nu = (__data.g.product_stoich_coeff(k,r) -
              __data.g.reactant_stoich_coeff(k,r))
        if (nu <> 0):
            c.append((__data.g.species_name(k),nu))
    return c

def reactionString(r):
    return __data.g.reaction_type(r)

def fwdRatesOfProgress():
    return __data.g.forward_rates_of_progress

def revRatesOfProgress():
    return __data.g.reverse_rates_of_progress
            

# test
if __name__ == "__main__":

    from cantera import one_atm
    readMechanism('gri30.xml')
    setState(2000.0, one_atm, [1.0]*numSpecies())
    print __data.g

    print
    print elementNames()
    print indexElt('N')
    print elementAtomicWt()

    print
    print speciesNames()
    for s in speciesNames():
        print s,':'
        for e in elementNames():
            print numberOfElementXinSpeciesY(e, s),e,
        print

    for i in range(numReactions()):
        print
        print reactionString(i)
        print specCoeffsInReaction(i)

    qf = fwdRatesOfProgress()
    qr = revRatesOfProgress()
    qnet = qf - qr
    print qf
    print qr
    print qnet
    
    

    

    
