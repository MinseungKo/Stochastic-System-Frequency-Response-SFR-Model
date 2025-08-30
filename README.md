# Stochastic-System-Frequency-Response-SFR-Model
MATLAB Code for the paper "A Stochastic Framework for Voltage Dependency and Uncertainty-Integrated Frequency Stability Assessment."

## Code Overview
This code aims to enhance the system frequency response (SFR) model, which is widely used for assessing power system frequency.
Three improvements are made as follows:

### 1) Unified SFR model with IBRs
In addition to synchronous generators and fast-frequency response resources, grid-following and grid-forming inverters are added with realistic control settings.

### 2) Load bus voltage generation method
Based on the generator bus voltage, we estimate the load bus voltage to incorporate the voltage-dependent behaviors of the load and DER.
The conventional estimation method uses the fixed voltage sensitivity of the generator to the load bus, which represents the load bus voltage as a linear combination of the generator bus voltages.
Our method enhances this approach by utilizing the generator-to-generator sensitivity, which has a linear relationship with the generator-to-load sensitivity.
This enables us to obtain the time-variant sensitivity, resulting in a more accurate load bus voltage estimation.

### 3) Stochastic SFR model
We modify the conventional deterministic SFR to incorporate stochastic elements, thereby accounting for various types of power system uncertainties. Solutions can be obtained by solving stochastic differential equations.
Specifically, uncertainties implemented in this code are as follows:
  i) Steady-state uncertainties:
  Load amount, IBR output, BTM DER distribution & output
  ii) Dynamic uncertainties:
  Load model parameters, Voltage sensitivity, Voltage jump, Completely unknown uncertainties

## File Description
- Stochastic_SFR.m: main file to run stochastic SFR model
- loadbus_3_500.xlsx: PSS/E simulation output on IEEE 39-bus system. A 500 MW load increase at Bus 3 is simulated.
- volt_sensitivity.csv: Voltage sensitivity obtained by injecting reactive power perturbation signals at generator buses. Simulations are conducted by using PSS/E.
