# Bounded One-Dimensional Plasma-Simulation-Code


This is a bounded one-dimensional kinetic particle-in-cell code with Monte-Carlo collisions and secondary electrons that is capable of accounting for any type of species, as well as any number of different types. Users supply functions or constants that return the collision cross-sections and secondary electron yields. Sample Argon and atomic Hydrogen functions are supplied. Code features include:

- Atomic Collisions
  - Elastic electron-atom scattering
  - Electron-impact ionization
  - Elastic ion-atom scattering (temporarily implement binary scattering)
  - Ion-atom charge-exchange collisions
- Plasma-Surface Interactions
  - Electron-surface reflection
  - Electron-surface back-scattering
  - Electron-impact secondary electrons
  - Ion-impact secondary electrons
- Two bounded walls with user defined constant or time-dependent voltages
- Optional steady state boundary where an electron-ion pair is put into the simulation domain everytime an electron arrives at a surface. (Must use injection source with this option)
- Optional injection source with user defined injection options

- Optional oblique magnetic field
- When upper particle limit is reached the population is reduced and particles are re-weighted. The user has the following options:
  - Select the number of velocity groups
  - Select the maximum number of particles allowed per group
  - Select the amount to decrease the population by
- Option to pick up a code run from where the previous one left off
- User defines system and gas parameters via input files
- Built in diagnostics calculates graph data by taking averages over a period of time that is determined by the user(user determines start time and duration). This graph data is printed to text files for use by data visualization applications.Current diagnostics include:
  - Rho, (net charge density)
  - Phi, (electrostatic potential)
  - Electric field
  - Averaged electron velocity through domain
  - Averaged ion velocity through domain
  - Averaged particle velocity distributions at locations chosen by the user
  - Instantaneous quantities in time:
    - Total particle number of each species
    - Electrostatic potential at three locations chosen by the user
    - Charge density for each species at two locations specified by the user



This code is still a work in progress. Future developements include:

- Have the user define everything in a single input file, via third party software. Currently several input files are used to define parameters and optional user-defined functions.
- The code which decides whether or not an electron collision occurs is still a work in progress due to the need to limit the number of ionization events. Currently an elastic collision occurs if a random number between zero and one is less than the elastic probability and likewise for ionization. In the case of ionization a further condition that the number of ionization events is less than a user defined limit is also an option as there are current issues with runaway ionization.
- Clean up code: some variables are re-declared in different classes.
- The current method that returns a particle velocity randomly taken off a Maxwellian distribution is a temporary place holder for a more accurate technique.
- Option to watch real time graphics data.
- User defined graph quantities.
- More collision types.
- The ability to handle relativistic velocities

REFERENCES:
S J Araki and R E Wirz, IEEE Trans. Plasma Sci., vol 41 (2013) 470
CK Birdsall and AB Langdon, Plasma Physics Via Computer Simulation
CK Birdsall, IEEE Trans. Plasma Sci., vol 19 (1991) 65
Peter Burger, Phys. Fluids, vol 10, (1967) 658
P Diomede et al, American Institute of Chemical Engineers J., vol 59, (2013) 3214
P Diomede et al, Plasma Sources Sci. Technol. , vol 14 (2005) 459
Z Dunko et al, Plasma Phys. Control Fusion, vol 54 (2012) 124003
Feng Shi et al, Physics of Plasmas, vol 15 (2008) 063503
Ahok Jain et al, Physica Scripta, vol 41 (1990) 321
R K Janev and J J Smith, Atomic and Plasma-Material Interaction Data for Fusion, vol 4 (1993) 
E Kawamura et al, Plasma Sources Sci. Technol., vol 8 (1990) R45
Ralf Krinke and H M Urbassek,  J. Phys. D: Appl. Phys., vol 29 (1996) 378
Particle Weighting:  E E Kunhardt and Y Tzeng, J. Comp. Phys., vol 67 (1986) 279
EE Kunhardt et al, Journal of Computational Physics, vol 67 (1986) 279
Kenichi Nanbu, IEEE Trans. Plasma Sci., vol 28 (2000) 971
K Nanbu and Y Kitatani, Vacuum, Elsevier Science, vol 47 (1996) 1023
K Nanbu and S Yonemura, J. of Computational Physics, vol 145 (1998) 639
A V Phelps, J. Phys. Chem. Ref. Data, vol 20 (1991) 557
Electron-Argon cross-sections: A V Phelps https://jila.colorado.edu/~avp/collision_data/
Richard J Procassini and Charles K Birdsall, Phys. Fluids B: Plasma Physics, vol 3, (1991) 1876
M Radmilovic-Radjenovic et al, Journal of Physics: Conference Series, vol 71 (2007) 012007
M Radmilovic-Radjenovic and Z L Petrovic, Eur. Phys. J. D, vol 54 (2009) 445
J Rodney and M Vaughan, IEEE Transactions on Electron Devices, vol 36 (1989) 1963
Tomonori Takizuka and Hirotada Abe, Journal of Computational Physics, vol 25 (1977) 205
V Vahedi et al, Plasma Sci. Technol., vol 2 (1993) 261
V Vahedi  and M Surendra, Computer Physics Communications, vol 87 (1995) 179
J P Verboncoeur et al, J. Comp. Phys., vol 104, (1993) 321
Argon ion- atom and charge-exchange energy dependent collision cross-sections: pdp1(www.eecs.berkeley.edu)
Electron-impact secondary yield:  Xoopic (www.eecs.berkeley.edu)
Electron-H atom and Ion-H atom elastic energy dependent cross-sections are Excel data fits to data from Aladin database 
