# SEDEX
SEDEX is a public pipeline written in Python to carry out SED fitting using MARCS and BOSZ spectral libraries. Based on its current configuration, SEDEX takes Gaia DR3 source ID,  $T_{\rm{eff}}$, $\sigma_{T_{\rm{eff}}}$, $\textrm{log}g$, and $\sigma_{\textrm{log}g}$ as input, where $T_{\rm{eff}}$ is required to lift the $T_{\rm{eff}}-A_V$ degeneracy. SEDEX automatically downloads broadband photometry via Gaia archive. The current supported photometric systems are Gaia ( $G$, $G_{\rm{BP}}, G_{\rm{RP}}$), 2MASS ( $J, H, K_S$), ALLWISE ( $W1, W2, W3, W4$), Hipparcos ( $H_{P}$), Tycho2 ( $B_T, V_T$), SDSS ( $u, g, r, i, z$), Pan-STARRS1 ( $g, r, i, z, y$), APASS ( $V, B, g, r, i$), and SkyMapper ( $u, v, g, r, i, z$), consisting of 34 bandpasses.

The output of SEDEX includes bolometric flux, angular radius, extinction, luminosity, and radius. Extinction is a direct output, bolometric flux
is calculated by integrating the best-fitting spectrum, luminosity is derived by combining bolometric flux and distance, radius is computed from luminosity and input $T_{\rm{eff}}$, and angular radius is inferred from radius and distance. More details can be found in [Yu et al. 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv220600046Y/abstract). Please cite this paper if you use it in your research.

## Instructions on how to run SEDEX.
- Create a Gaia account and input your user name and password in the file $\texttt{inlist.py}$. 
- Download [MARCS](https://marcs.astro.uu.se/) (standard composition) and [BOSZ](https://archive.stsci.edu/hlsp/bosz/search.php) (solar [C/M], resolution $R=20000$, ~87 GB)  spectral libraries.
- Prepare a list of targets in the file, $\texttt{Data/Input\\_Fits/TOIs/UserInputData.csv}$, with these columns: starID, teff, tefferr, logg, loggerr, feh, feherr, Av, and Averr. Note that starID has to be Gaia DR3 source_id. Av and Averr can be NaN. The sample name, TOIs, can be optionally changed in the file $\texttt{inlist.py}$.
- Configure the inlist file, $\texttt{inlist.py}$. Suggest to use default values.
- Run $\texttt{submit.sh}$, if the single star mode is set, or run $\texttt{submit\\_slurm.sh}$ for HPC, if the batch mode is set.
- Derived parameters are damped into this file:
	$\texttt{Data/Output\\_Fits/TOIs/SEDFits/MARCS/Output\\_SED\\_Fits\\_Final.csv}$.
 The figures of the Blackbody fitting and SED fitting are stored at $\texttt{Figures/TOIs/}$.
