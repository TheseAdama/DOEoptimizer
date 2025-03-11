# 📦 DOEoptimizer

**DOEoptimizer** is an R package that provides four optimization algorithms specifically designed to optimize criteria (matrix input functions) for physical experiments.  
It offers a comprehensive set of tools and utilities for constructing optimal experimental designs using various optimization techniques, including *genetic algorithms*, *simulated annealing*, *stochastic optimization*, and *greedy algorithms*.

## 📥 Installation

You can install the latest version of the package manually or directly from GitHub.

### Option 1: Install from GitHub (Recommended)

Make sure you have the `devtools` package installed, then use:

```r
install.packages("devtools")
devtools::install_github("TheseAdama/DOEoptimizer")
```
### Option 2: Manual Installation (Download & Install ZIP)

1. **Download the ZIP or TAR.GZ file**
   Download the latest version of the package in ZIP or TAR.GZ format.

   - For **Windows**: `DOEoptimizer_x.y.z.zip`
   - For **Linux/macOS**: `DOEoptimizer_x.y.z.tar.gz`

2. **Install the package manually in R**

   Open your R session and run one of the following commands, replacing the file path with where you downloaded the archive:

   - **On Windows**:
     ```r
     install.packages("path/to/DOEoptimizer_x.y.z.zip", repos = NULL, type = "win.binary")
     ```

   - **On Linux/macOS**:
     ```r
     install.packages("path/to/DOEoptimizer_x.y.z.tar.gz", repos = NULL, type = "source")
     ```

---

After installation, load the package:

```r
library(DOEoptimizer)
```
