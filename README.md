Bayesian Network [![Build Status](https://travis-ci.org/godai0519/BayesianNetwork.svg?branch=master)](https://travis-ci.org/godai0519/BayesianNetwork)
===========

What's this
-----------
This software is for Bayesian Network (BN) Library.  
A implementation of Bayesian Networks Model for pure C++14, included Loopy-BP and Likelihood Weighting.  
Development is not yet finished, but you can already use.

このソフトウェアはベイジアンネットワークのためのライブラリです．  
C++14を使用して実装しており，Loopy-BPやLikelihood Weighingといったベイジアンネットワークモデルを提供します．


Install
-------
This library is implemented as Header Only Library.  
Just through the `/path/to/BayesianNetwork` directory. (ex. `-I/path/to/BayesianNetwork`)

header onlyライブラリとして実装してあるので，ディレクトリへのパスを通すだけで使用可能です．  
(ex. `-I/path/to/BayesianNetwork`)


Compilers Tested
----------------
* Linux:
    + GCC, C++14: 4.8.1, 4.9.2
    + Clang, C++14: 3.5, 3.7.0 (trunk)
* Windows:
    + Visual C++: 12.0 (CTP_Nov2013), 13.0


State of progress (進捗どうですか？)
------------------------------------
Implementation of Belief Propagation (BP) become stable.
Maybe, return value will be changed from 2-D array to 1-D array.

ダメです．
確率伝播(BP)法の実装が安定しました．戻り値を1次元に落とすかもしれないです．
