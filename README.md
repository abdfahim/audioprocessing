# Common functions for Spherical Harmonic decomposition

These Matlab libraries were created during my PhD research at the Australian National University in Canberra, Australia. If you would like to make reference to the libraries, you can cite any of my following works:

1. **[PhD Thesis](https://openresearch-repository.anu.edu.au/handle/1885/203506):** Fahim, Abdullah, "Spatial dissection of a soundfield using spherical harmonic decomposition," PhD thesis, The Australian National University, Canberra, Australia, 2020.
2. **[Journal 1](https://ieeexplore.ieee.org/abstract/document/8357900):** Fahim, Abdullah and Samarasinghe, Prasanga N. and Abhayapala, Thushara D., "PSD estimation and source separation in a noisy reverberant environment using a spherical microphone array," IEEE/ACM Transactions on Audio, Speech, and Language Processing, vol. 26, no. 9, pp. 1594–1607, 2018.
3. **[Journal 2](https://ieeexplore.ieee.org/abstract/document/8936994):** Fahim, Abdullah and Samarasinghe, Prasanga N. and Abhayapala, Thushara D., "Multi-Source DOA Estimation through Pattern Recognition of the Modal Coherence of a Reverberant Soundfield," IEEE/ACM Transactions on Audio, Speech, and Language Processing, vol. 28, pp. 605–618, 2019.


## File Description:
* **SHTools.m:** This is the base class for all Spherical Harmonics-related functions
* **exmaples_SHTools.m:** Demonstrates some of the usages including:
  * Get real and complex SH under various circumstances
  * Get analytical alphas and betas under various circumstances
  * Simulate sound pressure  under various circumstances
  * Perform Spherical Harmonic decomposition
  * Bessel and Hankel functions
  * Other miscellaneous useful functions
  * Get Eigenmic coordinates

* **HOAGrid.m & HOAGridData.mat:** Design spherical microphone array with weights.

* **STFTClass:** Perform efficient STFT and inverse STFT for audio processing.
* **examples_STFTClass.m:** demonstrate usage of STFTClass
