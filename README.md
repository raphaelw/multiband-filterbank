# Multiband Filterbank
Audio filterbank for real-time audio signal processing implemented in C++. You can choose:
* Number of frequency bands
* Filter order
* Stopband attenuation

*The project also contains a basic Up-/Downsampler. See `Resampling.hpp`*

## Details
Implements a filterbank that is allpass complementary meaning that it has a flat frequency response with some unnoticable phase distortion.
> **Filterbank Structure:** R. Cassidy and J. O. Smith, "A tunable, nonsubsampled, non-uniform filter bank for multi-band audition and level modification of audio signals," Conference Record of the Thirty-Eighth Asilomar Conference on Signals, Systems and Computers, 2004., 2004, pp. 2228-2232 Vol.2, doi: <https://doi.org/10.1109/ACSSC.2004.1399563>

Each filterbank stage uses a crossover filter based on IIR filters of type EMQF (Elliptic with minimal Q-factor).
> **EMQF Filter Design:** Miroslav D. Lutovac, "Filter Design for Signal Processing Using MATLAB and Mathematica", 2000, Prentice Hall, ISBN-13: 978-0201361308

## Integration
* Simple to use header-only C++ library. Just include the corresponding files from the `afx` directory.
* Note that it requires BOOST library for some special mathematical functions. `bcp` tool may be used for the following files:
    * `boost/math/special_functions/jacobi_elliptic.hpp`
    * `boost/math/special_functions/ellint_1.hpp`

## TODO
* Implement better filter structure that is more suitable for time-varying operation. (Vadim Zavalishin's TPT)
* Prevent denormal numbers driving up CPU load when filter states are decayed.
