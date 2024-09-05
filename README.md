# Intelligent-and-Adaptive-Control-Systems

## Direct Model Reference Adaptive Control of an Inverted Pendulum

**Author:** Nikolaos Konstas  
**Institution:** Department of Electrical & Computer Engineering, Aristotle University of Thessaloniki  
**Date:** January 2023  

## Project Overview

This project focuses on the design and simulation of Direct Model Reference Adaptive Control (DMRAC) methods to stabilize an inverted pendulum system. The study involves three different adaptive control strategies to achieve precise stabilization and tracking of the pendulum's angle using MATLAB.

### Part A: Reference Model Control
- Developed a reference model controller to stabilize the pendulum at the zero position.
- Designed control laws based on system linearization and simulated the performance for various initial conditions.

### Part B: State Feedback DMRAC
- Implemented an adaptive state feedback controller that adjusts its parameters in real-time to ensure stability.
- Simulated the controller's performance for different scenarios, showing improved tracking accuracy with varying system dynamics.

### Part C: Output Feedback DMRAC
- Developed an adaptive output feedback controller for cases where only limited measurements are available.
- Evaluated the controller's robustness and its ability to adapt to different system conditions.

### Part D: Robustness Analysis
- Performed a robustness study to evaluate the sensitivity of the three controllers to changes in free parameters.
- Analyzed how variations in controller parameters and model specifications affect system performance, identifying the optimal settings for different scenarios.

## Results

- **Reference Model Control:** Achieved stabilization of the pendulum with minimal error and ensured compliance with the desired performance specifications.
- **State Feedback DMRAC:** Demonstrated accurate tracking of reference signals with adaptive parameter adjustments, resulting in reduced output errors.
- **Output Feedback DMRAC:** Successfully controlled the system under limited feedback conditions, maintaining robust performance despite parameter variations.
- **Robustness Analysis:** Identified the ranges of parameter settings that maximize performance across different operating conditions.

## Report

For a detailed explanation of the project's objectives, methods, and results, please refer to the [full project report](Report.pdf).

## Technologies Used

- **Languages:** MATLAB  

## How to Run

1. **Download the MATLAB files** from this repository.
2. **Run the corresponding MATLAB files** for each part of the project
