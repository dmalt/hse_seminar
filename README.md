Instructions
------------

### Study phase lag influence on the algorithm performance
How does phase lag affect DICS and iDICS performance?



#### Instructions
###### Load data here:
https://drive.google.com/open?id=0BycSyZwtFDq9b3lQblZaSm5aYlk

https://drive.google.com/open?id=0BycSyZwtFDq9V002Y2VCSTh4TXc

https://drive.google.com/open?id=0BycSyZwtFDq9bUlvWTQ5M0ZOZnM

https://drive.google.com/open?id=0BycSyZwtFDq9WTdIRERfRXFYU0U
###### Seminar presentation:
https://github.com/dmalt/hse_seminar/blob/master/pres.pdf

1. Fix signal to noise ratio (InducedScale)
2. Generate artificial MEG recordings for phase shifts 0, pi/4 and pi/2 (function: SimulateData)
3. Run DICS and imaginary DICS on simulated data (function: DICS; to run DICS with imaginary part modification use the flag 'is_imag')
4. Generate Precision-Recall curves for each run and plot them on the same figure.

#### Questions
1.What changes in DICS performance when phase shift increases from 0 to pi/2? What about iDICS?

2.Which of the two algorithms shows better overall performance? Explain the observed result.
