# Analytical Self-Consistent Field

Python package to calculate non-charged polymer brush parameters.
Such as polymer brush thickness, volume fraction and osmotic pressure
at given parameters e.g. distance from grafting surface.

## Abstract

We consider a polymer brush with degree of polymerization N and grafted to an impermeable for the polymer planar surface. Polymer chains are densely grafted to the surface by one of the ends. Number of polymer chains per unit area (expressed in the polymer segment length) σ is a measure of grafting density.

As it is shown in ref[] analytical strong-stretching self-consistent field (SS-SCF) approximation segment results in parabolic segment potential dependency on the distance z from grafting surface.

Mean field Flory-Huggins approximation expresses segment potential in terms of the polymer density. Together with normalization condition ∫<sub>0</sub><sup>H</sup>ϕ(z)dz = N⋅σ  and vanishing osmotic pressure condition Π(z=D)=0 it defines the polymer density ϕ as a function of the distance z from grafting surface.


- log(1 - φ) - 2⋅χ⋅φ - log(1 - φ<sub>H</sub>) - 2⋅χ⋅φ<sub>H</sub> = (3π<sup>2</sup>/8N<sup>2</sup>)(H<sup>2</sup>-z<sup>2</sup>)

On the l.h.s. of the equation is difference in segment potential at a given distance z and segment potential at the edge of the brush, where H is the brush thickness,  ϕ<sub>H</sub> = ϕ(z=H) and  χ<sub>PS</sub> - polymer-solvent interaction parameter.
We define insertion free energy penalty ΔF as a change in free energy when a particle is moved from an infinite distance from brush to a given distance from the grafting surface. Then the free energy penalty can be approximated as


## Usage

### CLI

...

### Python

After the import you can get a function to evaluate
volume fraction, osmotic pressure, insertion free energy penalty, diffusion coefficient, partition coefficient for given parameters.
The next example shows volume concentration calculations

```python
import scf_pb
#calculates local polymer density for sigma 0.02 at grafting surface
phi = scf_pb.phi(N=1000, chi=0.3, z=0, sigma = 0.02)
#0.16132539706803806

#grafting density we want to make calculations with
sigma = [0.01, 0.02, 0.03, 0.04]
#vectorized version of phi function
phi_arr = scf_pb.phi_v(N=1000, chi=0.3, z=0, sigma = sigma)
#array([0.10497596658134178, 0.16132539706803806, 0.20601691513583287,
#       0.24411579950858975], dtype=object)


chi = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
#iterates arguments as cartesian product and calls function, respectively
#(see itertools.product)
phi = scf_pb.phi_v(N=1000, chi=chi, z=0, sigma = sigma)
#array([[0.08060265992912352, 0.12645725063302338, 0.16405649995483096,
#        0.19694738362842357],
#       [0.09388632716508993, 0.14589613297822057, 0.1878441112564706,
#        0.22405062758851813],
#       [0.12342079608503442, 0.18499709443342405, 0.23243025723917354,
#        0.2721155247109588],
#       [0.2691000006100053, 0.30850381016135153, 0.34472994328065776,
#        0.37683369597899286],
#       [0.5305800047307337, 0.5353703249555996, 0.5427356789442022,
#        0.5520003019555626],
#       [0.6842552355869056, 0.685597935046323, 0.6877871390726218,
#        0.6907558236981812]], dtype=object)

len(sigma)
#4
len(chi)
#6
np.shape(phi)
#(6,4)

```

In this example we calculate phi-profile for the all integer values of distance from zero
to the brush edge.

```python
import scf_pb
#calculates brush thickness
D = scf_pb.D(N=1000, chi=0.6, sigma = 0.02)
#69.00098952395685
z = np.arange(D)
#array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12.,
#       13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25.,
#       26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38.,
#       39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51.,
#       52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64.,
#       65., 66., 67., 68., 69.])
phi = scf_pb.phi_v(N=1000, chi=0.6, sigma = 0.02, z=z)
#array([0.30850381016135153, 0.3084937853359111, 0.30846370573486126,
#       0.30841355597239717, 0.3083433103667146, 0.3082529328816076,
#       0.3081423770442657, 0.3080115858388366, 0.307860491575189,
#       0.30768901573216123, 0.30749706877443317, 0.30728454994204335,
#       0.30705134701134706, 0.3067973360260927, 0.3065223809970705,
#       0.3062263335686066, 0.3059090326499181, 0.30557030400913965,
#       0.30520995982754273, 0.30482779821114725, 0.3044236026566486,
#       0.3039971414681756, 0.30354816712097576, 0.3030764155677217,
#       0.302581605482545, 0.30206343743738007, 0.30152159300451575,
#       0.30095573377852713, 0.30036550030990206, 0.29975051094173955,
#       0.2991103605397848, 0.2984446191048372, 0.29775283025510235,
#       0.29703450956445254, 0.2962891427406058, 0.2955161836250567,
#       0.29471505199405323, 0.29388513113686143, 0.2930257651841517,
#       0.29213625615521976, 0.29121586068792504, 0.2902637864095957,
#       0.2892791879003753, 0.28826116219249265, 0.28720874373932315,
#       0.2861208987766324, 0.2849965189844804, 0.28383441434142276,
#       0.2826333050421497, 0.2813918123245418, 0.2801084480211584,
#       0.2787816026117551, 0.27740953150552694, 0.2759903392215435,
#       0.2745219610597288, 0.27300214175752624, 0.2714284105025936,
#       0.2697980515098032, 0.26810806915877206, 0.2663551464075118,
#       0.26453559482249844, 0.26264529405752823, 0.260679617917308,
#       0.25863334317328157, 0.2565005359317917, 0.2542744083906787,
#       0.2519471359464657, 0.2495096203205286, 0.24695117780681652,
#       0.24425912143820783], dtype=object)
```

## Installation
To install the package run in your python environment
```
pip install [path_to_the_package]
```

## Notes
The package is written as a part of PhD thesis work of the author.