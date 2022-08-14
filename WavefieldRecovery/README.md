# WavefieldRecovery

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> WavefieldRecovery

It is authored by Yijun.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

Or, you could setup the environment with following dependencies:

## Dependencies

The minimum requirements for theis software, and tested version, are `Python 3.x` and `Julia 1.6.0`.
This software requires the following dependencies to be installed:

- [MAT](https://github.com/JuliaIO/MAT.jl). This library can read MATLAB .mat files, both in the older v5/v6/v7 format, as well as the newer v7.3 format.
- [JOLI](https://github.com/slimgroup/JOLI.jl),Julia framework for constructing matrix-free linear operators with explicit domain/range type control and applying them in basic algebraic matrix-vector operations.
- [GenSPGL](https://github.com/slimgroup/GenSPGL.jl). A Julia solver for large scale minimization problems using any provided norm.
- [SeisJOLI](https://github.com/slimgroup/SeisJOLI.jl). Collection of SLIM in-house operators based on JOLI package.
- [Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl). Julia wrapper for the arpack library designed to solve large scale eigenvalue problems.
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/). Julia provides native implementations of many common and useful linear algebra operations which can be loaded with using LinearAlgebra. 

First, you need install the following packages from the stable master branch:
```
using Pkg; Pkg.add(PackageSpec(url="https://github.com/slimgroup/JOLI.jl"))
using Pkg; Pkg.add(PackageSpec(url="https://github.com/slimgroup/SeisJOLI.jl"))
using Pkg; Pkg.add(PackageSpec(url="https://github.com/slimgroup/GenSPGL.jl"))
```

then install dependencies:
```
Pkg.add("MAT")
Pkg.add("Arpack")
Pkg.add("LinearAlgebra")
```

## Software
This software is divided as follows:

 The ground truth data collected from the Gulf of Suez [Full.mat](https://slim.gatech.edu/PublicationsData/zhang2020SEGwrw/).

*data/*:
 
 This directory contains Jittered subsampling indexes[ind.mat].
 
*scripts/*: 

 This directory contains codes to run the corresponding experiments.You can run the 'RecursiveLR.jl' to reproduce the experiments, also you can modify the files to change the settings and design your own experiment.
 
 ```
 Test.jl #The main function of our experiments 
 SubFunctions/ConventionalLR.jl #The subfunction of conventional LR method for each frequency slice 
 SubFunctions/WeightedLR.jl #The subfunction of weighted method for each frequency slice 
 SubFunctions/NLfunForward_weight.jl #The weighted forward function to implement the data misfit constraint
 ```
 
## Citation

If you find this software useful in your research, we would appreciate it if you cite:

```bibtex1
@conference {zhang2020SEGwrw,
	title = {Wavefield recovery with limited-subspace weighted matrix factorizations},
	year = {2020},
	note = {Accepted in SEG},
	month = {4},
	abstract = {Modern-day seismic imaging and monitoring technology increasingly rely on dense full-azimuth sampling. Unfortunately, the costs of acquiring densely sampled data rapidly become prohibitive and we need to look for ways to sparsely collect data, e.g. from sparsely distributed ocean bottom nodes, from which we then derive densely sampled surveys through the method of wavefield reconstruction. Because of their relatively cheap and simple calculations, wavefield reconstruction via matrix factorizations has proven to be a viable and scalable alternative to the more generally used transform-based methods. While this method is capable of processing all full azimuth data frequency by frequency slice, its performance degrades at higher frequencies because monochromatic data at these frequencies is not as well approximated by low-rank factorizations. We address this problem by proposing a recursive recovery technique, which involves weighted matrix factorizations where recovered wavefields at the lower frequencies serve as prior information for the recovery of the higher frequencies. To limit the adverse effects of potential overfitting, we propose a limited-subspace recursively weighted matrix factorization approach where the size of the row and column subspaces to construct the weight matrices is constrained. We apply our method to data collected from the Gulf of Suez, and our results show that our limited-subspace weighted recovery method significantly improves the recovery quality.},
	keywords = {algorithm, data reconstruction, frequency-domain, interpolation, Processing, SEG},
	url = {https://slim.gatech.edu/content/wavefield-recovery-limited-subspace-weighted-matrix-factorizations},
	author = {Yijun Zhang, Shashin Sharan, Oscar Lopez, Felix J. Herrmann}
}
```

```bibtex2
@conference {zhang2019SEGhfw,
	title = {High-frequency wavefield recovery with weighted matrix factorizations},
	booktitle = {SEG Technical Program Expanded Abstracts},
	year = {2019},
	note = {(SEG, San Antonio)},
	month = {09},
	pages = {3959-3963},
	abstract = {Acquired seismic data is normally not the fully sampled data we would like to work with since traces are missing due to physical constraints or budget limitations. Rank minimization is an effective way to recovering the missing trace data. Unfortunately, this technique{\textquoteright}s performance may deteriorate at higher frequency because high-frequency data can not necessarily be captured accurately by low-rank matrix factorizations albeit remedies exist such as hierarchical semi-separable matrices. As a result, recovered data often suffers from low signal to noise ratio S/Rs at the higher frequencies. To deal with this situation, we propose a weighted recovery method that improves the performance at the high frequencies by recursively using information from matrix factorizations at neighboring lower frequencies. Essentially, we include prior information from previously reconstructed frequency slices during the wavefield reconstruction. We apply our method to data collected from the Gulf of Suez, which shows that our method performs well compared to the traditional method without weighting.},
	keywords = {algorithm, data reconstruction, frequency-domain, interpolation, Processing, SEG},
	doi = {10.1190/segam2019-3215103.1},
	url = {https://slim.gatech.edu/Publications/Public/Conferences/SEG/2019/zhang2019SEGhfw/zhang2019SEGhfw.html},
	presentation = {https://slim.gatech.edu/Publications/Public/Conferences/SEG/2019/zhang2019SEGhfw/zhang2019SEGhfw_pres.pdf},
	author = {Yijun Zhang and Shashin Sharan and Felix J. Herrmann}
}
```

```bibtex3
@article {kumar2014GEOPemc,
	title = {Efficient matrix completion for seismic data reconstruction},
	journal = {Geophysics},
	volume = {80},
	number = {5},
	year = {2015},
	note = {(Geophysics)},
	month = {09},
	pages = {V97-V114},
	abstract = {Despite recent developments in improved acquisition, seismic data often remain undersampled along source and receiver coordinates, resulting in incomplete data for key applications such as migration and multiple prediction. We have interpreted the missing-trace interpolation problem in the context of matrix completion (MC), and developed three practical principles for using low-rank optimization techniques to recover seismic data. Specifically, we strive for recovery scenarios wherein the original signal is low rank and the subsampling scheme increases the singular values of the matrix. We use an optimization program that restores this low-rank structure to recover the full volume. Omitting one or more of these principles can lead to poor interpolation results, as we found experimentally. In light of this theory, we compensate for the high-rank behavior of data in the source-receiver domain by using the midpoint-offset transformation for 2D data and a source-receiver permutation for 3D data to reduce the overall singular values. Simultaneously, to work with computationally feasible algorithms for large-scale data, we use a factorization-based approach to MC, which significantly speeds up the computations compared with repeated singular value decompositions without reducing the recovery quality. In the context of our theory and experiments, we also find that windowing the data too aggressively can have adverse effects on the recovery quality. To overcome this problem, we carried out our interpolations for each frequency independently while working with the entire frequency slice. The result is a computationally efficient, theoretically motivated framework for interpolating missing-trace data. Our tests on realistic 2D and 3D seismic data sets show that our method compares favorably in terms of computational speed and recovery quality with existing curvelet- and tensor-based techniques.},
	keywords = {2D, 3D, algorithm, interpolation, low-rank, signal processing},
	doi = {10.1190/geo2014-0369.1},
	url = {https://slim.gatech.edu/Publications/Public/Journals/Geophysics/2015/kumar2014GEOPemc/kumar2014GEOPemc.pdf},
	url2 = {http://library.seg.org/doi/abs/10.1190/geo2014-0369.1},
	author = {Rajiv Kumar and Curt Da Silva and Okan Akalin and Aleksandr Y. Aravkin and Hassan Mansour and Ben Recht and Felix J. Herrmann}
}

```

## Contact

For questions or issue, please contact yzhang3198@gatech.edu.


