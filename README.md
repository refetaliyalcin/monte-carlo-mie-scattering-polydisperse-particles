code works by executing main.m

variables should be altered inside main.m file, they are described with comments inside the code

works on MATLAB

can work in parallel or on a single core

solution of radiative transfer equation in; one layer, pigmented, plane parallel medium using Lorenz Mie theory

ray is incident from air to coating. coating is coated on a substrate

substrate could be air or other material such as silver, glass etc.

the code estimates spectral hemispherical reflectance, transmittance and absorptance

can handle; independent scattering, boundary reflections, absorption in medium.

can't handle; coherent backscattering, dependent scattering and polarized ray tracing.

while calculating the scattering direction the code uses henyey greenstein function approximation

a similar version works on GPU is here: https://github.com/refetaliyalcin/ColorRCMC

if you encounter any problem you can contact me from refetali@gmail.com

in case of use, please cite the article as:, Colored Radiative Cooling Coatings with Nanoparticles, ACS Photonics, Refet Ali Yalçın, Etienne Blandre, Karl Joulain, Jeremie Drevillon, https://doi.org/10.1021/acsphotonics.0c00513
