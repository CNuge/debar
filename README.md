# seqDenoise
## An R package for the denoising coi5p barcode data generatd through Pacific Biosciences single molecule, real-time (SMRT) sequencing.
[![Build Status](https://travis-ci.com/CNuge/seqdenoise.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/seqdenoise)
--------------------------------------------------------

`seqDenoise` is an R package designed to denoise DNA barcode data generated through [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) on [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). The package uses taxonomic specific profile hidden Markov models (PHMMs) to identify and correct insertion and deletion errors within individual circular consensus sequences (CCS). When multiple CCSs are avaliable for a sample, the corrected sequences are then evaluated to determine the most likely true sequence.




