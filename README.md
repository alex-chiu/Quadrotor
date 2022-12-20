# Quadrotor
Code used to model and simulate the dynamics and control of a small quadrotor drone

## Features
* High-fidelity linear and attitude dynamics models that take motor spin-up time and drag into account.
* Simple PD controller used to direct the quadrotor.
* GNSS, camera, and IMU data processing code.
* An Unscented Kalman Filter (UKF) used to estimate the quad's location, and a Wahba solver used to initialize it.
* A* path finding algorithm used to plan the quadrotor's path, and a computer vision algorithm (RANSAC) used to identify target locations.

## License
MIT

Contributors:

Alex Chiu - https://github.com/alex-chiu