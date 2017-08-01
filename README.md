# ImageReconstruction
A repository of MATLAB code for understanding simple phase retrieval and image reconstruction 
algorithms as well as extracting the point spread functions of real data.

PhaseRetrieval.m is an algorithm taken from Gonsalves, 1982, "Phase Retrieval and diversity in Adaptive Optics." It is intended to extract the phase of an observed point spread function. It creates a point spread function that is aberrated by some amount using Zernike polynomials, and then runs a loop to minimize the error between the observed point spread function and an estimated point spread function.

JointEstimate.m is also taken from Gonsalves' paper. It utilizes the same concept at PhaseRetrieval, but is slightly more advanced and is aimed at reconstructing an image that is aberrated. It first creates an object, and then "images" it with some amount of aberration; we mostly looked at defocus aberration in our simulation. It then creates another image that has an a priori amount of additionaly aberration (we knew how much more defocused the second image was). Using this, it generates a point spread function to match each image, and then reconstructs the image to one without aberration.

These two algorithms use a slightly modified steepest descent algorithm to minimize an error metric between the observed and estimated quantities. The error metric is dependent upon the coefficients of the zernike polynomials, which can end up being a large parameter space.

zernikeps.m is just a function that calculates the zernike polynomials as a matrix.

dotadv.m is a proof of concept script. We wanted to know if the object that was being imaged had an idealized point, could we use the image of that point (i.e. the point spread function) to reconstruct and defocus the image. At least on the computer, this turned out to work perfectly.
