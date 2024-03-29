# Hamiltonian Spectrum Alignment for partial shape similarity [[Paper]](https://arxiv.org/abs/1906.06226)

This repository is the MATLAB implementation of the paper:

**Correspondence-Free Region Localization for Partial Shape Similarity via Hamiltonian Spectrum Alignment**
<br> Arianna Rampini, Irene Tallini, Maks Ovsjanikov, Alex M. Bronstein, Emanuele Rodolà <br>
*Best Paper Award* at [3DV, 2019](https://www.computer.org/csdl/proceedings-article/3dv/2019/313100a037/1ezRALztN1m).

<p align="center">
  <img src="./teaser.PNG" width="500">
</p>

To reproduce the result in this figure, run the script ```run.m``` in MATLAB.

## Optimization

In this code we used the toolbox [manopt](https://www.manopt.org/) for optimization on manifolds.

## Data

The shapes used in this example (in ```data```) are from the [SHREC’16 Partiality benchmark](https://www.dais.unive.it/~shrec2016/dataset.php).
