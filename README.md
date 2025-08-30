# Stochastic-System-Frequency-Response-SFR-Model
MATLAB Code for the paper "A Stochastic Framework for Voltage Dependency and Uncertainty-Integrated Frequency Stability Assessment."

This code aims to enhance the system frequency response (SFR) model, which is widely used for assessing power system frequency.
Three improvements are made as follows:
1) Unified SFR model with IBRs: In addition to synchronous generators and fast-frequency response resources, grid-following and grid-forming inverters are added with realistic control settings.
2) Load bus voltage generation method: Based on the generator bus voltage, we estimate the load bus voltage to incorporate the voltage-dependent behaviors of the load and DER.
3) Stochastic SFR model: We modify the conventional deterministic SFR into a stochastic SFR to incorporate the various types of power system uncertainties. Solutions can be obtained by solving stochastic differential equations.
