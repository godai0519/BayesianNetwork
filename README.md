# Bayesian Network [![Build Status](https://travis-ci.org/godai0519/BayesianNetwork.svg?branch=master)](https://travis-ci.org/godai0519/BayesianNetwork)

## What's this
This software is for Bayesian Network (BN) Library.  
A implementation of Bayesian Networks Model for pure C++11, included Loopy-BP and Likelihood Weighting.  
Development is not yet finished, but you can already use.

このソフトウェアはベイジアンネットワークのためのライブラリです．  
C++14を使用して実装しており，Loopy-BPやLikelihood Weighingといったベイジアンネットワークモデルを提供します．


## Install
This library is implemented as Header Only Library.  
Just through the `/path/to/BayesianNetwork` directory. (ex. `-I/path/to/BayesianNetwork`)

header onlyライブラリとして実装してあるので，ディレクトリへのパスを通すだけで使用可能です．  
(ex. `-I/path/to/BayesianNetwork`)


## Compilers Tested
* Linux:
    + GCC, C++11: 4.8.1, 4.9.2, 5.0.0
    + Clang, C++11: 3.5, 3.7.0 (trunk)
* Windows:
    + Visual C++: 12.0 (CTP_Nov2013), 13.0


## Feature
### Bayesian Network structure learning
* Brute Force class
    + search all patterns, very very slow.
* Greedy class
    + some patterns(not become good score) ignore, fast a little.
* Simulated Annealing class
    + as if a metal, with spending time, structure is less likely to change.
* K2 Algorithm class
    + give information that you already know some relations between nodes.
* Stepwise Structure
    + first, make some cluster that include some node, one after another cluster combine two clusters.

##### Available Algorithms in Stepwise Structure
|                   | Learning in Cluster | Learning between Cluster |
|------------------:|:-------------------:|:------------------------:|
|        Brute Force|          OK         |             OK           |
|             Greedy|          OK         |             OK           |
|Simulated Annealing|          OK         |            *NG*          |
|       K2 Algorithm|          OK         |            *NG*          |


### Inferring unobserved variables
* Loopy-BP class
    + high-speed and strong in a loop graph.
* Rejection Sampling (a.k.a. Logic Sampling) class
    + very simple, but very slowly.
* Likelihood Weighting
    + higher than Rejection Sampling, and extremely accuracy.

## Author and Contact
Feel free to contact me ;)  
Bugs and issues are reportable below:
* [GitHub Issue](//github.com/godai0519/BayesianNetwork/issues)
* [Twitter](//twitter.com/godai_0519)
* [Facebook](//www.facebook.com/godai.azuma)
* [Blog](//d.hatena.ne.jp/godai_0519/)


## Licence
Code released under [the MIT license](//github.com/godai0519/BayesianNetwork/blob/master/LICENSE).
