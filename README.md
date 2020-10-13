# Lifting-Line-Theory
Implementation of Prandtl's Lifting Line Theory. To setup the problem, the following steps are required.

1. Initialize the class. No inputs are required.
2. Load the geometry parameters: Aspect ratio, taper ratio, surface, the chord distribution (as a function of the station angle) and the angle of attack of zero lift.
3. Solve the problem by calling the solve method of this class, which requires as inputs the angle of attack (this term should include any effect regarding the geometric twist), the station angles and the terms of the Fourier Series to be considered for the computation of the terms. For example, if only terms 1 and 3 are to be considered, this input should read [1, 3] and coefficients A1 and A3 would be computed.
