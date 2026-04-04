# MGMR
Semi analytic calularion of radio emission from extensive airshowers that includes the effects of atmospheric electric fields

# How to install:
git clone https://github.com/galaxyvita/pyMGMR3D.git

cd pyMGMR3D/Library

./fftpack5.1d_MGMR.sh (the warnings are normal)

cd ..

cd program

make -f MGMR3D_fit-makefile-v5.make

# Run an example:

cd run

./pyMGMR3D.sh SIM000001.in


# References
* A Macroscopic Description of Coherent Geo-Magnetic Radiation from Cosmic Ray Air Showers
by O. Scholten, K. Werner, F. Rusydi (https://arxiv.org/abs/0709.2872)

* Macroscopic Treatment of Radio Emission from Cosmic Ray Air Showers based on Shower Simulations
by Klaus Werner, Olaf Scholten (https://arxiv.org/abs/0712.2517)

* Macroscopic Model of Geomagnetic-Radiation from Air Showers
by Olaf Scholten, Klaus Werner (https://arxiv.org/abs/0808.1959)

* The Lateral Distribution Function of Coherent Radio Emission from Extensive Air Showers; Determining the Chemical Composition of Cosmic Rays
by Krijn D. de Vries, Ad M. van den Berg, Olaf Scholten, Klaus Werner (https://arxiv.org/abs/1008.3308)

* The convergence of EAS radio emission models and a detailed comparison of REAS3 and MGMR simulations
by Tim Huege, Marianne Ludwig, Olaf Scholten, Krijn D. de Vries (https://arxiv.org/abs/1009.0346)

* Macroscopic Model of Geomagnetic-Radiation from Air Showers, II
by Olaf Scholten, Krijn D. de Vries, Klaus Werner (https://arxiv.org/abs/1010.5268)

* Macroscopic Geo-Magnetic Radiation Model; Polarization effects and finite volume calculations
by Krijn D. de Vries, Olaf Scholten, Klaus Werner (https://arxiv.org/abs/1010.5364)

* Coherent Cherenkov Radiation from Cosmic-Ray-Induced Air Showers
by K.D. de Vries, A.M. van den Berg, O. Scholten, K. Werner (https://arxiv.org/abs/1107.0665)

* A Realistic Treatment of Geomagnetic Cherenkov Radiation from Cosmic Ray Air Showers
by Klaus Werner, Krijn D. De Vries, Olaf Scholten (https://arxiv.org/abs/1201.4471)

* A detailed comparison of REAS3 and MGMR simulations for radio emission from EAS
by Marianne Ludwig, Tim Huege, Olaf Scholten, Krijn D. de Vries (https://arxiv.org/abs/1202.3352)

* What the radio signal tells about the cosmic-ray air shower
by Olaf Scholten, Krijn D. de Vries, Klaus Werner (https://arxiv.org/abs/1207.1874)

* The air shower maximum probed by Cherenkov effects from radio emission
by Krijn D. de Vries, Olaf Scholten, Klaus Werner (https://arxiv.org/abs/1304.1321)

* The cosmic-ray air-shower signal in Askaryan radio detectors
by Krijn D. de Vries, Stijn Buitink, Nick van Eijndhoven, Thomas Meures, Aongus O'Murchadha, Olaf Scholten (https://arxiv.org/abs/1503.02808)

* Determining Atmospheric Electric Fields using MGMR3D
by T. N. G. Trinh, O. Scholten, S. Buitink, K. D. de Vries, P. Mitra, T. Phong Nguyen, D. T. Si (https://arxiv.org/abs/2203.03134)

* Reconstructing Air Shower Parameters with MGMR3D
by P. Mitra, O. Scholten, T. N. G. Trinh, S. Buitink, J. Bhavani, A. Corstanje, M. Desmet, H. Falcke, B. M. Hare, J. R. Hörandel, T. Huege, N. Karastathis, G. K. Krampah, K. Mulrey, A. Nelles, H. Pandya, S. Thoudam, K. D. de Vries, S. ter Veen (https://arxiv.org/abs/2307.04242)

* Aperture correction for Beamforming in radiometric detection of ultra-high energy cosmic rays
by O. Scholten, T. N. G. Trinh, S. Buitink, A. Corstanje, B.M. Hare, T. Huege, B.V. Jhansi, K. Mulrey, A. Nelles, H. Schoorlemmer, S. Thoudam, P. Turekova, K. de Vries (https://arxiv.org/abs/2411.12324)






