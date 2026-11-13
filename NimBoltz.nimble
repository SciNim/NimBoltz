# Package

version       = "0.1.0"
author        = "Vindaar"
description   = "Electron transport in gas mixtures: Magboltz 11.21 in Nim, on CPUs, GPUs and in the browser"
license       = "MIT"

# Dependencies

requires "nim >= 2.4.0"
requires "npeg"
requires "measuremancer >= 0.2.5"
# Rational powers of units (sqrt of any unit), exact scales: not released yet
requires "https://github.com/SciNim/Unchained#rational-powers"
requires "shell"
