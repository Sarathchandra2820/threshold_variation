import numpy as np
import math
#import scm.plams
#from scm import plams
from scm.plams import *
import sys

if len(sys.argv) != 4:
    print("Usage: python gwrose_script.py <molecule> <distance> <method>")
    sys.exit(1)

molecule = sys.argv[1]  # Molecule type (pyrene, pyridine, ethylene)
add_virt_threshold = float(sys.argv[2])
frag_th = float(sys.argv[3])  # Dimer separation distance
#method = "qsGW"
theta = 0
method = "qsGW"  # Calculation method (G0W0, evGW, TDDFT)
distance = 3
num_exc = 5

init()


@add_to_class(AMSResults)
def get_excitations(results,symmetry=''):
    """Returns excitation energies (in eV) and oscillator strenghts (in Debye)."""
    if results.job.ok():
        exci_energies_au = results.readrkf("Excitations SS "+symmetry, "excenergies", file="engine")
        oscillator_str_au = results.readrkf("Excitations SS "+symmetry, "oscillator strengths", file="engine")
        # The results are stored in atomic units. Convert them to more convenient units:
        exci_energies = Units.convert(exci_energies_au, "au", "eV")
        oscillator_str = Units.convert(oscillator_str_au, "au", "Debye")
        return exci_energies, oscillator_str
    else:
        return [], []

@add_to_class(AMSResults)
def get_coupling_matrix(results,symmetry='A -1 -1'):
    """Returns coupling matrix (in eV) between localized and CT state. Diagonal has the uncoupled excitation energies."""
    if results.job.ok():
        n_states = results.readrkf("Excitations SS "+symmetry, "nr of excenergies", file="engine")
        coupling_matrix_au = results.readrkf("Excitations SS "+symmetry, "electronic couplings", file="engine")
        # The results are stored in atomic units. Convert them to eV and reshape into square matrix:
        coupling_matrix = np.array(Units.convert(coupling_matrix_au, "au", "eV")).reshape((n_states,n_states))
        return coupling_matrix
    else:
        return []

@add_to_class(AMSResults)
def get_coupled_excitations(results,symmetry='A -1 -1'):
    """Returns coupled excitation energies (in eV)."""
    if results.job.ok():
        n_states = results.readrkf("Excitations SS "+symmetry, "nr of excenergies", file="engine")
        coupling_matrix_au = results.readrkf("Excitations SS "+symmetry, "electronic couplings", file="engine")
        # The results are stored in atomic units. Convert them to eV and reshape into square matrix:
        coupling_matrix = np.array(Units.convert(coupling_matrix_au, "au", "eV")).reshape((n_states,n_states))
        eigenvalues, eigenvectors = np.linalg.eigh(coupling_matrix)
        return eigenvalues, eigenvectors
    else:
        return []
    

def rotation_matrix(axis, angle_degrees):
    """
    Generate a 3x3 rotation matrix for a given axis and angle (in degrees).
    """
    angle = np.radians(angle_degrees)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    cos = np.cos(angle)
    sin = np.sin(angle)
    return np.array([
        [cos + x*x*(1 - cos),     x*y*(1 - cos) - z*sin, x*z*(1 - cos) + y*sin],
        [y*x*(1 - cos) + z*sin,   cos + y*y*(1 - cos),   y*z*(1 - cos) - x*sin],
        [z*x*(1 - cos) - y*sin,   z*y*(1 - cos) + x*sin, cos + z*z*(1 - cos)]
    ])
    

# Generic settings for all calculations, used like this for the moving fragment 2.

settings = Settings()
#settings.input.ams.system.symmetrize = ''
settings.input.ams.task = "SinglePoint"
settings.input.ADF.basis.type = 'TZP'
settings.input.adf.symmetry = "NOSYM"
settings.input.adf.basis.core = "None"
settings.input.adf.numericalQuality = 'VeryGood'
#settings.input.adf.xc.LIBXC = "HYB_GGA_XC_CAMY_B3LYP"
settings.input.adf.xc.LIBXC = "LRC_WPBEH"
#settings.input.adf.xc.HYBRID = "CAMY-B3LYP"
#settings.input.adf.xc.XCFUN = ""
settings.input.adf.noprint = "Logfile"
settings.input.ADF.dependency.enabled = 'Yes'
settings.input.ADF.dependency.bas = 0.005
settings.input.ADF.RIHartreeFock.DependencyThreshold = 0.005
#settings.input.ADF.GW.enabled = 'Yes'
settings.input.adf.save = "TAPE15"

# Define GW settings dynamically based on method
if method == "G0W0":
    settings.input.ADF.GW.SelfConsistency = 'G0W0'
    settings.input.ADF.GW.nstates = "-1"
    settings.input.ADF.GW.DIIS = "3"
    settings.input.ADF.GW.scissorshift = 'True'
    settings.input.ADF.GW.nIterations = '15'
    if molecule!="pyrene":
        settings.input.adf.MBPT.formalism = 'RI'
        settings.input.adf.MBPT.frequencyGridType = 'GAUSSLEGENDRE'
        settings.input.adf.MBPT.nFrequency = '16'
        settings.input.adf.MBPT.useGreenXGrids = 'True'
        settings.input.adf.MBPT.nTime = '32'
    


elif method == "evGW":
    settings.input.ADF.GW.SelfConsistency = 'evGW'
    settings.input.ADF.GW.nstates = "-1"
    settings.input.ADF.GW.DIIS = "3"
    settings.input.ADF.GW.scissorshift = 'True'
    settings.input.ADF.GW.nIterations = '30'
    if molecule!="pyrene":
        settings.input.adf.MBPT.formalism = 'RI'
        settings.input.adf.MBPT.frequencyGridType = 'GAUSSLEGENDRE'
        settings.input.adf.MBPT.nFrequency = '16'
        settings.input.adf.MBPT.useGreenXGrids = 'True'
        settings.input.adf.MBPT.nTime = '32'
    

elif method == "qsGW":
    settings.input.ADF.GW.SelfConsistency = 'qsGW'
    settings.input.ADF.GW.nstates = "-1"
    settings.input.ADF.GW.DIIS = "3"
    settings.input.ADF.GW.scissorshift = 'True'
    settings.input.ADF.GW.nIterations = '30'

# Fragment 1 settings (run only once for isolated monomer excitation energies)
frag1_settings = settings.copy()
frag1_settings.input.ADF.TDA = ''
if method != "TDDFT":  # Ensure BSE is only included for GW methods
    frag1_settings.input.ADF.excitations.BSE = ''
frag1_settings.input.ADF.excitations.OnlySing = ''
frag1_settings.input.ADF.excitations.Lowest = num_exc

# Reference supermolecule (with symmetry) settings
refsm_settings = settings.copy()
refsm_settings.input.ADF.TDA = ''
if method != "TDDFT":  # Ensure BSE is only included for GW methods
    refsm_settings.input.ADF.excitations.BSE = ''
refsm_settings.input.ADF.excitations.OnlySing = ''
refsm_settings.input.ADF.excitations.descriptors = ''
refsm_settings.input.ADF.excitations.Lowest = num_exc
refsm_settings.input.ADF.MODIFYEXCITATION.UseOccVirtRange = '0 1'
#refsm_settings.input.ADF.excitations.Exact = {}  
#refsm_settings.input.ADF.excitations.Exact.End = ''

# Supermolecule (without symmetry) settings with intrinsic fragment orbitals
sm_settings = refsm_settings.copy()
sm_settings.input.ADF.Rose.nfragments = 2
sm_settings.input.ADF.Rose.Additional_virtuals_cutoff = add_virt_threshold
sm_settings.input.ADF.Rose.Frag_threshold = frag_th


if method in ["G0W0", "evGW"]:  # Correct condition check
    sm_settings.input.ADF.Rose.GWenergies = ''

if method == "qsGW":  # Fixed missing colon
    sm_settings.input.ADF.Rose.GWorbitals = ''

sm_settings.input.ADF.ExcitonTransfer.useRose = ''
sm_settings.input.ADF.ExcitonTransfer.Localize = 'OccupiedAndVirtual'
sm_settings.input.ADF.ExcitonTransfer.FullRun = ''
sm_settings.input.ADF.ExcitonTransfer.SecondOrder = 'True'
sm_settings.input.ADF.ExcitonTransfer.Output = 'AllAndFilteredCouplings'



# Load molecule
geometry_file = f"geometry/{molecule}.xyz"
frag1 = Molecule(geometry_file)
#frag1.translate([distance/2,0.0,0.0])

#try: plot_molecule(frag1) # plot molecule in a Jupyter Notebook in AMS2023+
#except NameError: pass
print(frag1)

frag1job = AMSJob(molecule=frag1, settings=frag1_settings, name='fragment_1')

# set to True if you want to see the input that is used
if False:
    print("-- input to the job --")
    print(frag1job.get_input())
    print("-- end of input --")

frag1job.run()

#exci_energies, exci_vectors = frag1job.results.get_excitations()
#print(exci_energies)

#distances = [3.0,4.0,5.0,6.0,8.0,10.0]


#sm_excitation_energies = []
#sm_excitation_vectors = []
#sm_excitation_couplings = []

axis = ([1,0,0])

frag2 = frag1.copy()
frag2.translate([distance,0.0,0.0])

R = rotation_matrix(axis, theta)

# Apply rotation
frag2.rotate(R)  # This now works

frag2job = AMSJob(molecule=frag2, settings=frag1_settings, name=f"frag2_{distance:.2f}")
#frag2job = AMSJob(molecule=frag2, settings=settings, name=f"frag2_{distance}")
frag2job.run()

# symmetries = ['A.g','B1.g','B2.g','B3.g','A.u','B1.u','B2.u','B3.u']
# refsm_excitation_energies  = {}
# refsm_oscillator_strengths = {}
# for symm in symmetries:
#     refsm_excitation_energies[symm] = []
#     refsm_oscillator_strengths[symm] =[]


# refsm = frag1 + frag2
# refsmjob = AMSJob(molecule=refsm, settings=refsm_settings, name=f"supermolecule_{distance:.2f}")
# refsmjob.run()
# for symm in symmetries:
#         try:
#             exci_energies, oscillator_str = refsmjob.results.get_excitations(symm)
#             if not isinstance(exci_energies, list):
#                 exci_energies, oscillator_str = [exci_energies], [oscillator_str]
#         except:
#             exci_energies, oscillator_str = [],[]
#         refsm_excitation_energies[symm].append(exci_energies)
#         refsm_oscillator_strengths[symm].append(oscillator_str)
#exci_energies, oscillator_str = refsmjob.results.get_excitations()
# with open(f'output_{distance}.txt', 'w') as f:
#     for key, value in refsm_excitation_energies.items():
#         f.write(f'{key}: {value}\n')


'''
refsm = frag1 + frag2
refsmjob = AMSJob(molecule=refsm, settings=refsm_settings, name=f"supermolecule_{distance:.2f}")
refsmjob.run()
exci_energies, oscillator_str = refsmjob.results.get_excitations()
print(exci_energies)
np.savetxt(f"output_{distance}.txt", exci_energies)
'''
        

sm = frag1 + frag2
print(sm)




#smjob = AMSJob(molecule=sm, settings=sm_settings, name=f"sm_fragments_{distance:.2f}")
smjob = AMSJob(molecule=sm, settings=sm_settings, name=f"sm_fragments")
@add_to_instance(smjob)
def prerun(self):
    import shutil
    from os.path import join
    shutil.copy(frag1job.results.rkfpath(file='engine'),join(self.path,'frag1.rkf'))
    shutil.copy(frag2job.results.rkfpath(file='engine'),join(self.path,'frag2.rkf'))
smjob.settings.link_files = False # to prevent an error when rerunning calculations we explicitly copy files instead of linking
sm_settings.input.ADF.Rose.FragmentFile = "frag1.rkf"
sm_settings.input.ADF.Rose.FragmentFile = "frag2.rkf"
smjob.run()
exci_couplings = smjob.results.get_coupling_matrix()

#sm_excitation_energies.append(exci_energies)
#sm_excitation_vectors.append(exci_vectors)
exci_couplings=exci_couplings.reshape((4*num_exc,4*num_exc))


np.savetxt(f"output_{add_virt_threshold : .2f}_{frag_th : .2f}.txt", exci_couplings)

finish()

