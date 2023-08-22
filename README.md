# [ChainConsumer](https://samreay.github.io/ChainConsumer)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/eefe9aa7d4904306877be1e17b952f39)](https://www.codacy.com/app/samuelreay/ChainConsumer?utm_source=github.com&utm_medium=referral&utm_content=Samreay/ChainConsumer&utm_campaign=badger)
[![Build Status](https://img.shields.io/travis/Samreay/ChainConsumer.svg)](https://travis-ci.org/Samreay/ChainConsumer)
[![Coverage Status](https://codecov.io/gh/Samreay/ChainConsumer/branch/master/graph/badge.svg)](https://codecov.io/gh/Samreay/ChainConsumer)
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/dessn/abc/blob/master/LICENSE)


[![PyPi](https://img.shields.io/pypi/v/ChainConsumer)](https://pypi.python.org/pypi/ChainConsumer)
[![Conda](https://anaconda.org/samreay/chainconsumer/badges/version.svg)](https://anaconda.org/samreay/chainconsumer)
[![DOI](https://zenodo.org/badge/23430/Samreay/ChainConsumer.svg)](https://zenodo.org/badge/latestdoi/23430/Samreay/ChainConsumer)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.00045/status.svg?style=flat)](http://dx.doi.org/10.21105/joss.00045)

A library to consume your fitting chains! Produce likelihood surfaces,
plot your walks to check convergence, output a LaTeX table of the
marginalised parameter distributions with uncertainties and significant
figures all done for you, or throw in a bunch of chains from different models
and perform some model selection!

[Click through to the online documentation](https://samreay.github.io/ChainConsumer)

### Installation

Install via `pip`:

    pip install chainconsumer

### Python Versions

Due to dependencies dropping support for Python 2.7, from 0.31.0 and onwards, only Python 3 will be supported.

Previous versions will still be installable for Python 2.7 environments, however I'd strong recommend upgrading.

### Contributors

I would like to thank the following people for their contribution in issues, algorithms and code snippets
which have helped improve ChainConsumer:

* Chris Davis (check out https://github.com/cpadavis/preliminize)
* Joe Zuntz
* Scott Dedelson
* Elizabeth Krause
* David Parkinson
* Caitlin Adams
* Tom McClintock
* Steven Murray
* J. Michael Burgess
* Matthew Kirby
* Michael Troxel
* Eduardo Rozo
* Warren Morningstar


### Common Issues

Users on some Linux platforms have reported issues rendering plots using ChainConsumer.
The common error states that `dvipng: not found`, and as per [StackOverflow](http://stackoverflow.com/a/32915992/3339667)
post, it can be solved by explicitly install the `matplotlib` dependency `dvipng` via `sudo apt-get install dvipng`.

If you are running on HPC or clusters where you can't install things yourself, users may run into issues where LaTeX or other optional dependencies aren't installed. In this case, set `usetex=False` in `configure` to request matplotlib not try to use TeX. If this does not work, also set `serif=False`, which has helped some uses.

