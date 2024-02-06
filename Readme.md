# 2D Analytic Signal Processing for Ultrasound C-scan Images

Developed by Xiaoyu (Leo) Yang, a research fellow at Nanyang Technological University (NTU), this project is part of the Ultrasonics and NDE Research Group led by Zheng (David) Fan. It focuses on implementing 2D analytic signal processing techniques to enhance ultrasound C-scan images, leveraging FFT for efficient computation and providing tools for detailed image analysis.

## Features

- Efficient FFT-based processing for 2D images.
- Advanced convolution operations in the frequency domain.
- Analytic signal computation for in-depth image analysis.

## Dependencies

- Qt 6.4.3+ for GUI and core functionality.
- FFTW3 for Fourier transform operations.

## Setup and Build

1. Ensure Qt and FFTW3 are installed.
2. Clone the repository, create a build directory inside the project folder, and navigate into it.
3. Configure and build the project using CMake:

   ```bash
   cmake .. && cmake --build .
Usage
Run the generated executable to process sample data, showcasing the project's capabilities in ultrasound image processing.

Contributing
Contributions are welcome. Feel free to submit pull requests or open issues for discussion.

## Copyright

Copyright © 2024 Ultrasonics and NDE Research Group, Nanyang Technological University (NTU). All rights reserved.