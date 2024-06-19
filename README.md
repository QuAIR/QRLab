# QRLab

**QRLab** is a MATLAB toolbox for exploring quantum information processing and quantum resource theory.


## Features
- **Entanglement Theory**: 

    - *Static Entanglement Measure*: Tempered Logarithmic Negativity $E_{
\mathrm{\tau}}$, Rains Bound $R$, MaxRainsEntropy $R_{
\mathrm{max}}$, Logarithmic Negativity $E_\text{N}$, $E_{
\mathrm{PPT}}$, $E_
\mathrm{eta}$

    - *Dynamic Entanglement Measure*: Max Logarithmic Negativity, Max Rains information

    - *Quantum Capacity*

- **Coherence Theory**: 

    - *Static Coherence Measure*: Robustness of Coherence

    - *Channel Simulation*: Simulating non-free operations via resource states

- **Magic Theory**: 

    - *Static Magic Measure*: Robustness of Magic (qubit), Magic Mana (qudit), max/min Thuama (qudit) 

    - *Representative Magic State Generation*


- **Quasi-Theory**: 

    - Probabilistic error cancelation 

    - Observable dependent probabilistic error cancellation 

    - Observable dependent deterministic error cancelation 

    - Circuit Knitting 

    - Virtual Recovery

- **Supermap**: 
    - Quantum Switch (both kraus and choi), Apply Quantum Switch
    - Link Product

- **Seesaw Algorithms**: Algorithms for providing sub-optimal solutions for non-linear optimization problems. 

    - CHSH game 

    - LOCC protocol

- **Extra Functions**: 

    - Conditional quantum mutual information 

- **API Documents**:
    - API documents can be found in this website https://quair.github.io/QRLab/api/.

## Requirements
1. QETLAB == 0.9
2. CVX == 2.1


## Installation
1. Clone QRLab to your local machine.
2. Download QETLAB 0.9. You could download it from https://qetlab.com/.
3. Add QRLab and QETLAB to MATLAB's path​, through command
```matlab
addpath(genpath('...\QETLAB-0.9'));
addpath(genpath('...\QRlab'));
```
4. Download and install CVX 2.1. You could download it from https://cvxr.com/cvx/.
Install CVX on Windows
```matlab
cd yourpath\cvx;​
cvx_setup;
```
Install CVX on Linux or a Mac
```matlab
cd ~/MATLAB/cvx;​
cvx_setup;
```

5. To fully unlock the magic qubit related functions, it is necessary to install channel_magic 2.0 in https://github.com/jamesrseddon/channel_magic.

## Getting Started: Use Cases

To use the package, simply call the function you need with appropriate parameters. For example:

```matlab
% To calculate the logarithmic negativity for a given state
rho = MaxEntangled(2) * MaxEntangled(2)'; % Define a quantum state
LN = LogNeg(rho);
disp('Logarithmic Negativity:');
disp(LN);
```


## More functions are coming.


## Contributing

Contributions to expand and improve this package are welcome.

## Acknowledgement

We acknowledge the use of the [CVXQUAD](https://github.com/hfawzi/cvxquad) package, a tool for MATLAB-based convex optimization. This package offers several essential functions, such as calculating the von Neumann entropy and quantum relative entropy, which have been invaluable to our research.