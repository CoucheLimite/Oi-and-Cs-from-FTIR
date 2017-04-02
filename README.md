# Oi-and-Cs-from-FTIR
A small calculator of interstitial oxygen and substitution carbon concentration from FTIR meaasurement. This little GUI requires numpy, matplotlib, scipy and PyQT5. It has been tested only with Python 3.6, but I guess it should also work with other versions of Python.

The calculation alogorithim is based on two standards:

1. The American Society for Testing and Materials (ASTM) (1993), *Standard test method for Interstitial Atomic Oxygen Content of Silicon by Infrared Absorption*, Designation: F 1188- 93a
2. The American Society for Testing and Materials (ASTM) (2000), *Standard Test Method for Substitutional Atomic Carbon Content of Silicon by Infrared Absorption*, Designation: F1391 - 93 (Reapproved 2000)

However, these standards require the tested sample to be 2 mm, which is the same thcikness of the reference sample. In reality, espeacially for the photovoltaic application, this is not practical, Therefore, a modified method has been used to remove the impact of thickness difference between tested sample and reference sample by a factor, which is the thickness ratio of those two (more accurately, the ratio of optical path). By using this factor, the tested sample and reference sample can have different thickness. The referecen sample doesn't have to be 2 mm. However, a thicker referenec sample can ensure a stronger signal thus a more accurate reference spectrum. The optimal value of that factor can be obtained directly from the measured spectra by comparing the heights of the silicon phonon peaks. A more detailed description of this modified method can be found in: A. K. Søiland, *Silicon for solar cells*, Norwegian University of Science and Technology, 2004.

One issue of FTIR measurement is the free carrier absorption (FCA). For the reference sample, it is recommended to use high resistivity samples (above 10 Ωcm) to make sure the phonon peaks are not influenced by FCA. For the tested sample, if the resistivity is below 3 Ωcm for p-type (or below 1 Ωcm for n-type), then the sample spectrum can be affected by FCA. This little program also implemented a small function to fit and remove the FCA in the sample. It assumes a quadratic relation between FCA and the wavelength when doing the fitting.

Another issue is the surface finish of the sample. For the reference sample, it is recommended to have double side polished (DSP) surface to ensure a good signal. The same for the tested sample. However, that's not always the case, then the sample FTIR spectrum can be affected by the wavelenght depdendent refelectance (In some case, the spectrum can be very tilted). This little software uses a simple technique to remove this effect by a subtraction of straight line near the carbon or oxygen peaks. This simple but efficient method is proposed by Energy Research Centre of the Netherlands.

A short summary, if the tested sample is DSP as the reference sample but affected by FCA, then the FCA correction is recommended. Howevre, if the reflectance is impacting the spectrum, the straight line correction is recommened. In general, the stright line correction can be accurate enough for almost all cases.

If you have any trouble or question in using this program, please don't hesitate to contact me by Github or email: yan.zhu@student.unsw.edu.au

Scheduled updated in next version:
1. Add uncertainty to the calculation;
2. Add warning if the concentration is below the detection limit;
3. Add alogorithim to remove the ripples in the spectrum due to thin thickness of the sample;
4. Add a small help function.
