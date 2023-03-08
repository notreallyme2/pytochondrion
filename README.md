# pytochondrion

A simulator of mitochondrial energetics using `scipy`  
This model is based on the mathematical model of mitochondrial energetics described in [Edwards *et al.*]()

## Installation
```[bash]
pip install .
```

## Useful resources
[`scipy` cookbook](https://scipy-cookbook.readthedocs.io/)  
[Modelling a Zombie Apocalypse](https://scipy-cookbook.readthedocs.io/items/Zombie_Apocalypse_ODEINT.html)

Me:
"The model is 'working' now, in as much as the numerical routine is working just fine. My problem now is that the model is not steady, the NADH (and UQH) just 'bleeds out' of the system. If I had some ballpark state 4 figures I could clamp some of the parameters and see if I couldn't figure out the problem. "


BK:
"Here are the values for state 4 (kU = 0)


NADH(0)=1643.512
UQH(0)=1150.044
C2+(0)=54.913
O2(0)=240
Hi(0)=0.044625
ATPti=15410.1516
PIti=56498.3329
ATPte=2000.0397
PIte=9937.17346
ADPte=0.132158

The resulting vRESP is: 15.1202

The first test of your method would be to pass to state 3 (kU = about 700 a.u.; do not use too great values because then ATP and ADP are converted to AMP).

In my model there is a fixed contribution of delta-psi and delta-pH to delta-p. Now we are working on a more mechanistic relationship between the components of delta-p.

I have not used the model for isolated mitochondria for years - tell me if you are able to run it now."