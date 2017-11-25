Digital Signal Processing
=========
A MATLAB code for my DSP project 

Project assignment
---------

The file signal XXX.wav, where XXX is the progressive number associated to the student name,
contains a discrete-time signal of duration 28 seconds and with sampling frequency Fs = 48000 Hz.

It is a modulated signal (for double-sideband reduced-carrier transmission) of the form

y(nT ) = x1(nT ) + A1*cos(2πf1 nT ) + x2 (nT ) + A2 * cos(2πf2nT ) ,

where x1 (nT ) and x2 (nT ) are two real audio information signals with frequency band [10, 4000] Hz, f1 and
f2 are the frequencies of the two sinusoidal carriers, where 4000 Hz ≤ f1 < f2 ≤ 19000 Hz frequencies
f1 and f2 are chosen such that there is no frequency overlap between the two modulated components, and
specifically such that f2 − f1 ≥ 8500 Hz, and A1 and A2 are the amplitudes of the carriers.

Write a MATLAB procedure that finds the two carrier frequencies f1 , f2 and their amplitudes A1, A2
by means of a DFT analysis, extracts the two sinusoidal carriers by means of appropriate filtering with
two band-pass filters with very narrow bandwidths, and demodulates the two audio information signals
multiplying the modulated signal by the previously extracted carriers and filtering the results with a band-
pass filter having pass-band [10, 4000] Hz.
