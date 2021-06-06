## Radiofrequency AC Zeeman Trapping for Neutral Atoms ## 
I wanted to immortalize some of the good code I've written and neat videos I've generated, so that other atom enthusiasts can benefit from them. 

Files are separated by topic:

### Animations
One animation accompanies Chapter 8, showing a complex representation in 3-D of the current evolution in a lateral wire. 
We also show four Rabi map animations:
1) Matching Appendix C, as it's meant to be
2) time evolution of a Z-wire pulse
3) time evolution of a U-wire pulse
4) F=2 vs. F=1 for the same phase, displaying an anomaly

### Trap position
A mathematica notebook that gave the analytic solutions to trap position, in Chapter 6. 

### Old GUI
This file retains x-y B<sub>DC</sub> orientation, used for the pyramid trap. Only microwave solutions are exact in this GUI, but the RF doesn't consider the whole manifold, only one RF transition at a time. 
Users should run the .m file, rather than opening the figure or editing it in GUIDE.

### New GUI
This GUI locks B<sub>DC</sub> to the z direction, but includes an appropriate Hamiltonian solver for RF AC Zeeman energies. 

### Chip VNA Data
In case it's useful in the future, I wanted to include full S11 spectrum up to 20 GHz for the U-Z-U wires, as measured. 

### Examples
These are a handful of MATLAB codes that went into figures, or are central to the results in some way. 
They include:

1) The Silvester AC skin simulation for arbitrary rectangular conductor, plus the B<sub>x</sub>-field in a line above it. I elect to animate the results somewhat, by playing the data over a few cycles of phase. 
2) Breit-Rabi and AC Zeeman energy curves and Force plots
3) Just AC Zeeman energy
4) Just DC Zeeman energy and Clebsch-Gordan coefficients, over a range of B
5) A supplementary function, which gives C.G. coefficients and resonances, for given B


