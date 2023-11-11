# PAP
This is a Program for hydrodynamics Analysis of Potential flow (PAP) written in Fortran.

The program exploits a state-of-the-art method to conduct hydrodynamics analysis for floating body which involves the calculation of fluid potential (including incident wave potential, diffraction wave potential and radiation wave potential), based on which pressure and forces can be calculated. 

## Assumptions:
1. ideal fluid (no viscosity)
2. no vorticity
3. regular incident wave with linearized free surface boundary condition
4. infinite water depth
5. sigle floating body
6. no advance (the floating body oscillating around its mean position)

The process of hydrodynamics analysis entails the calculation of three-dimensional infinite water depth free-surface Green function induced by a pulsating source, which is the most complicated part of the procedure. Many scholars proposed different methods to calculate Green function and were devoted to enhance the accuracy, stability and speed. The method of this program is based on the work of Delhommeau and his open source software NEMOH. References are listed below.

## References
1. Xie et al, Comparison of existing methods for the calculation of the infinite water depth free-surface Green function for the wave-structure interaction problem, Applied Ocean Research 81 (2018) 150–163.
2. Aurélien Babarit and Gérard Delhommeau. Theoretical and numerical aspects of the open source BEM solver NEMOH. In Proceedings of the 11th European Wave and Tidal Energy Conference. Sept 2015.
3. Ruddy Kurnia and Guillaume Ducrozet. NEMOH: Open-source boundary element solver for computation of first- and second-order hydrodynamics loads in frequency domain. Computer Physics Communications, 292:108885, 2023. doi:10.1016/j.cpc.2023.108885.
4. J. N. Newman, The theory of ship motions, Advances in Applied Mechanics 18 (1978) 221-283.
5. J. N. Newman, Approximations for the Bessel and Struve functions, Mathematics of Computation 43 (1984) 551-556.
6. J. N. Newman, Algorithms for free-surface Green function, Journal of Engineering Mathematics 19 (1985) 57-67.

## How to use
make sure control tool Make and a fortran compiler (such as gfortran) are already installed. Then run the following command to compile

`make`

or

`make run`

Then you get the target *PAP*, type

`./PAP`

to run the program

if you want to clean the modules and object files, run

`make clean`

or if you also want to clear the target, run

`make cleanall`
