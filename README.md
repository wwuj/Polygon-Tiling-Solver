# Polygon-Tiling-Solver
## 📌 About This Project
This project aims to explore and compute the topological and geometric properties of polygon edge-to-edge tilings. Specifically, the program enumerates and verifies all possible edge and angle configurations for the following scenarios:
1. **Single $n$-gon**: Forming an orientable surface of genus $g$, or a non-orientable surface of genus $k$ under edge-to-edge tessellation.
2. **Two congruent $n$-gons**: Forming an orientable surface of genus $g$, or a non-orientable surface of genus $k$ under edge-to-edge tessellation.

**Algorithmic Constraints:**
* The program automatically filters out angle combinations that yield no positive solutions.
* Strict geometric rule: The degree of each vertex (i.e., the number of polygons meeting at a single vertex) must be $\ge 3$.

## 💡 Developer's Note
This program was independently developed by an author with a non-Computer Science background. Due to limited programming experience, the code formatting and naming conventions might not be strictly standardized. This is a long-term project, and I will continuously update and refactor the code structure in the future. Feedback and contributions are highly welcome!

## ⚙️ Requirements & Setup
The solver is built on a **MATLAB + C++ (MEX)** hybrid architecture, utilizing C++ multithreading to significantly accelerate the parallel search process.

Please ensure the following software is installed before running:
1. **[MATLAB](https://www.mathworks.com/)**: Version R2024b or higher is recommended.
2. **[Visual Studio 2022](https://visualstudio.microsoft.com/)**: Make sure the "Desktop development with C++" workload is checked during installation to provide the C++ compiler for MATLAB.
