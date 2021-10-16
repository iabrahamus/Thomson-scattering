# Thomson-scattering
## Relativistic high-intensity laser electron scattering

This project simulates the scattering of high energy photons from non-linear Thomson scattering. The interaction of a free and highly relativsitic electron beam with a realistic 3-dimensional laser field is simulated. In the code the user is provided with alternatives of using the paraxial approximation for calculating fields in the laser beam or calculate laser field parameters with accuracy as high as the 5th order non-paraxial correction term. 

### Project content instructions
[Detailed instructions on how to use the Thomson simulation code](How%20to%20use%20the%20Thomson%20Code.pdf)<br>
[laser parameters and the interaction geometry is adjusted by user through 'parameters.dat' file](parameters.dat)<br>
[The six dimensional phase space of the electron beam is initialized through (x,y,z,px,py,pz) 'ephase_space.dat' file](ephase_space.dat)<br>
[Use MakeFile to create required executables in workding directory](Makefile)<br>
[Python script to create Pegasus Workflow to distribute jobs to required cluster nodes](setup_job.py)<br>
[Main laser-electron interaction code is](Spectra.f)<br>

### Some publications as a result of the Thomson project 

1. https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.16.030705 (Isaac Ghebregziabher et al)
2. https://www.nature.com/articles/nphoton.2013.314 (Nathan Powers, Isaac Ghebregziabher, Gregory Goloving et al)
