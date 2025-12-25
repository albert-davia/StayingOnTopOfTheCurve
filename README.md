# French Sovereign Yield Curve Modeling

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-Academic-blue.svg)](#)

## Overview

This project implements a **Cascade Affine Term Structure Model** to analyze the French sovereign yield curve during the 2024 political crisis. Using state-of-the-art filtering techniques (Unscented Kalman Filter), we estimate latent factors driving interest rate dynamics and investigate how political uncertainty propagates through the term structure.

The core methodology is based on:

> Calvet LE, Fisher AJ, Wu L. _"Staying on Top of the Curve: A Cascade Model of Term Structure Dynamics."_ **Journal of Financial and Quantitative Analysis**, 2018;53(2):937-963. [DOI:10.1017/S0022109018000030](https://doi.org/10.1017/S0022109018000030)

## Key Features

- **Multi-factor cascade model**: Captures mean-reversion at multiple frequencies through a hierarchical structure
- **Unscented Kalman Filter**: Robust state estimation for nonlinear measurement equations
- **Maximum likelihood estimation**: Rigorous parameter inference with standard error computation
- **Yield curve forecasting**: Out-of-sample prediction capabilities with recalibration options

## Quick Start

### 1. Prerequisites

- MATLAB R2020a or later
- Optimization Toolbox (for `fminunc`)

### 2. Prepare Your Data

Create a `data/` folder and place your data file there:

```
your_repo/
├── data/
│   └── nusrates_dette_cleaned.mat   # Your yield data
├── estimationCANFCPv2S.m            # Main script
├── CANFCPv2.m                       # Model specification
└── ...
```

The `.mat` file should contain:

- `rates`: T × N matrix of observed yields (T observations, N maturities)
- `swapmat`: Vector of maturities in years
- `mdate`: Vector of MATLAB datenum values

### 3. Run the Estimation

```matlab
% Open MATLAB in the repository directory
>> estimationCANFCPv2S
```

### 4. Configuration

Edit the **CONFIGURATION** section at the top of `estimationCANFCPv2S.m`:

```matlab
% --- Path Configuration ---
DATA_DIR = 'data';           % Input data directory
OUTPUT_DIR = 'output';       % Results directory
FIGURES_DIR = 'figures';     % Plots directory

% --- Data Configuration ---
dataDette = 1;               % 1 = French debt, 0 = US LIBOR/swap

% --- Estimation Settings ---
estimation = 1;              % 1 = run optimization
stderror = 1;                % 1 = compute standard errors
AttemptNumber = '31';        % Run identifier
```

## Project Structure

```
├── estimationCANFCPv2S.m    # Main estimation script
├── CANFCPv2.m               # Core model: parameters & bond pricing
├── ukf_lfnlh.m              # Unscented Kalman Filter
├── ratelikefunlf.m          # Log-likelihood computation
├── liborswap.m              # Yield measurement equations
├── SigmaPoints.m            # UKF sigma points generation
├── scalecoef01.m            # Cascade coefficient computation
├── data/                    # [Create] Place your data here
├── output/                  # [Generated] Estimation results
└── figures/                 # [Generated] Plots
```

## Methodology

### The Cascade Affine Model

The model specifies latent state dynamics under the physical measure:

$$dX_t = \kappa(\theta - X_t)dt + \sigma dW_t$$

where the mean-reversion matrix $\kappa$ exhibits a cascade structure with frequencies scaling by a factor $b$:

$$\kappa_j = \kappa_1 \cdot b^{j-1}, \quad j = 1, \ldots, n$$

### Risk Pricing

The market price of risk is specified as:

$$\Lambda_t = \gamma_0 + \gamma_1 X_t$$

allowing for time-varying risk premia that depend on the state of the yield curve.

### Filtering

We employ the **Unscented Kalman Filter (UKF)** to extract latent factors from observed yields. The UKF handles the nonlinear bond pricing equations through sigma-point propagation, avoiding the linearization errors of the Extended Kalman Filter.

## Output Files

After running the estimation, you'll find:

| File                 | Description                |
| -------------------- | -------------------------- |
| `output/par_*.txt`   | Optimized model parameters |
| `output/mu_dd_*.txt` | Extracted latent factors   |
| `output/y_dd_*.txt`  | Model-fitted yields        |
| `figures/*.png`      | Visualization plots        |

## Core Components

| File                    | Description                                                                                       |
| ----------------------- | ------------------------------------------------------------------------------------------------- |
| `CANFCPv2.m`            | Defines model parameters, computes bond pricing coefficients using affine term structure formulas |
| `estimationCANFCPv2S.m` | Main script: loads data, runs optimization, computes standard errors, generates plots             |
| `ukf_lfnlh.m`           | Implements the Unscented Kalman Filter for state estimation                                       |
| `ratelikefunlf.m`       | Evaluates the log-likelihood function for parameter optimization                                  |
| `liborswap.m`           | Computes model-implied yields from latent factors                                                 |
| `SigmaPoints.m`         | Generates sigma points and weights for the UKF                                                    |

## Academic Context

This project was developed as part of a research project at **École Polytechnique** (2025), investigating the impact of political instability on French government bond markets during the 2024 political crisis.

## References

1. Calvet, L.E., Fisher, A.J., & Wu, L. (2018). Staying on Top of the Curve: A Cascade Model of Term Structure Dynamics. _Journal of Financial and Quantitative Analysis_, 53(2), 937-963.

2. Julier, S.J., & Uhlmann, J.K. (1997). A New Extension of the Kalman Filter to Nonlinear Systems. _Proceedings of AeroSense_.

