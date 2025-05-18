=========================================
README - Programming Assignment #3: Global Placement
=========================================

Student ID: B11901188  
Name: Yi-En Wu  
Course: Physical Design for Nanometer ICs, Spring 2025  
Instructor: Yao-Wen Chang  
Date: May 18th, 2025  

-----------------------------------------
1. Compilation
-----------------------------------------
To compile the project:

1. Open a terminal and navigate to the root folder of this project.
2. Run:

    make

This will compile the source files located in `src/` and generate an executable named `place` inside the `bin/` directory.

To clean and recompile:

    make clean
    make

-----------------------------------------
2. Execution
-----------------------------------------
To run the global placer:

    ./bin/place -aux <input.aux>

Example:

    ./bin/place -aux ./input/ibm01-cu85.aux

The output will be a `.gp.pl` file in the same folder as the input, and visualizations will be generated in `plot_output/`.

-----------------------------------------
3. Description of the Implementation
-----------------------------------------

This global placement tool aims to distribute standard cells on a chip to minimize Half-Perimeter Wirelength (HPWL) and congestion. The optimization approach is based on Conjugate Gradient Descent with enhancements:

- **Objective Function:** 
  A hybrid of wirelength and density cost. The density cost is computed using a bell-shaped potential function, smoothed with a Gaussian kernel, and enhanced with level-smoothing to control bin overflow.

- **Gradient Calculation:**  
  The gradient of the density function is computed from the cost derivative w.r.t. smoothed bin densities. The bell-shaped potential function is used to model a cell’s influence on nearby bins. Gradients are accumulated accordingly.

- **Dynamic Step Size:**  
  The optimizer adjusts the step size dynamically at each iteration based on the magnitude of the gradient direction to avoid divergence or stagnation.

- **Boundary Clamping:**  
  Modules are clamped to the chip outline after each step to ensure no cell escapes the chip region.

- **Visualization:**  
  Plots include:
  - Density heatmap (Gnuplot + PNG)
  - Cell scatter distribution
  - Combined plots per iteration

  These are stored in `plot_output/`.

-----------------------------------------
4. Input Format
-----------------------------------------
- `.aux` file: Entry point file pointing to the `.nodes`, `.nets`, and `.scl` files.
- `.nodes`: Module definitions (standard cell sizes).
- `.nets`: Netlist of pins connecting cells.
- `.scl`: Chip layout and bin configuration.

-----------------------------------------
5. Output
-----------------------------------------
- `<circuit>.gp.pl`: The output global placement file.
- `plot_output/`: Contains visualization files:
    - `density_*.png`: Density maps
    - `cells_*.png`: Cell distribution per iteration
    - `combined_*.png`: Side-by-side comparison
    - `grad_vectors.txt`: Gradient vectors per cell

-----------------------------------------
6. Known Issues / Notes
-----------------------------------------
- The Gaussian kernel used for smoothing in both forward and backward passes uses a kernel size of 5 and σ = 2.0.
- Density cost normalization is applied using a target density of 1.0.
- Overflowing modules may still appear in very congested regions; additional legalizers should be used post-placement.

-----------------------------------------
7. Contact
-----------------------------------------
Yi-En Wu  
National Taiwan University  
Department of Electrical Engineering  
Email: b11901188@ntu.edu.tw
