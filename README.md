# MoDELib2 fork from the original upstream

## Current Works

- Implement the stacking fault noise generation based on the MD correlation data.

- Create a quickstart guide
- Perform Normalized (devide NR) after the IFFT. Change the Value of NX and NY, because for FFTW, after FFT the size will be changed ( NX*NY -> NX*(NY/2+1))
  now the sampled AutoCorrelation are the same as Input one.
- Solving for SFNoise and SSNoise separately，Add a selectValues function in "UniformPeriodicGrid.h",
  change the value of gridSize and gridSpacing with label "SF".
  Made some modification of 'GlidePlaneNoise.cpp'.
