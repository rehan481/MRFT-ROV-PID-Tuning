# MRFT-ROV-PID-Tuning

MATLAB implementation of an MRFT-based system identification and PID tuning workflow for BlueROV2 ROV depth/attitude control.

This repository accompanies the manuscript:

**Identification and Control of ROV Attitude and Heave: A Compact Approach using Modified Relay Feedback Test**

The code is provided to support reproducibility of the MRFT-based identification and controller-tuning results reported in the paper.

## Overview

The workflow performs:

1. Conversion of experimental Modified Relay Feedback Test (MRFT) oscillation data into complex frequency-response samples.
2. Complex-domain fitting of a compact low-order depth-channel model.
3. Construction of the identified transfer function and state-space realization.
4. PID tuning and comparison using:
   - GA-based ISE tuning using the identified model
   - GM-beta MRFT tuning
   - Ziegler-Nichols tuning from MRFT
   - Tyreus-Luyben tuning from MRFT

## Main script

Run the main script in MATLAB:

```matlab
MAIN
