# Actuator Fault Analysis of the Otter Ship

This repository contains a MATLAB-based analysis of actuator faults in the Otter vessel. Our main goal is to demonstrate how to detect and estimate actuator faults using an Adaptive Extended Kalman Filter (AEKF).

## Project Overview

- **Matlab Code**: The file `m4.m` implements the Otter’s dynamic model, introduces actuator faults, and applies the AEKF for fault diagnosis.
- **Report**: The PDF (`AIS4004_Individual_Portfolio.pdf`) discusses the theoretical background, methodology, and results of the simulation.

## Key Features

1. **Otter Vessel Model**: A 3-DOF representation capturing surge, sway, and yaw dynamics.
2. **Cosine-Blending Control**: Smoothly varying control inputs to reduce transient effects.
3. **Fault Injection**: Scheduled changes in actuator effectiveness to simulate faults.
4. **AEKF Implementation**: An extended Kalman filter that adaptively estimates both states and fault parameters.

## How to Run

1. **Clone this repository**:
   ```bash
   git clone https://github.com/elisha136/Actuator-fault-Analysis-of-the-otter-ship.git
   Open m4.m in MATLAB (or MATLAB Online).
   ```

Run the script to generate plots and observe fault estimation performance.

Refer to the PDF (AIS4004_Individual_Portfolio.pdf) for a comprehensive explanation of the results.

Contributing
Feel free to fork this repository and submit pull requests. Any suggestions or improvements are welcome!

License
This project is for educational purposes. Check with the repository owner for licensing details.
