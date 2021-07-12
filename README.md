Class HrsTrkCorr  for finding corrections to tracks to account
for residuals in track angles for the reconstruction at the target.

R. Michaels, July 2021

HrsTrkCorr.C  -- implementaion
HrsTrkCorr.h  -- header
Makefile  -- use to complie the library, which can be used for
             codes like Monte Carlos

Necessary data input to HrsTrkCorr : the "holefiles" directory:
/holefiles  -- contains "hole files" and "residual files".
A hole file is the location in tg_th vs tg_ph space where the sieve
holes are located.
A residual file is the residual data -- deviation in mrad.
These exist for the Left HRS and Right HRS, for horizontal beam
positions of 0, -3 and +3 mm (approximately).

/monte  -- directory
A simple example of how to use the HrsTrkCorr library:
See the README file in /monte for further instructions.

There are also git repositories for

Finding the residuals
   https://github.com/rwmichaels/TrkHoles

Extensive testing of HrsTrkHoles
   https://github.com/rwmichaels/TrkCorrTest