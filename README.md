# shinybcells R Package

This repository contains the R package `shinybcells`, which includes large data files managed with Git Large File Storage (LFS). Follow the instructions below to properly install the package and run the application.

---

## 📦 Prerequisites

- **R** (>= 4.1)
- **devtools** R package
- **Git** and **Git LFS**

---

## 🔧 Install Git and Git LFS

### macOS

Install Git and Git LFS using [Homebrew](https://brew.sh/):

```bash
brew install git git-lfs
git lfs install
```

### Linux (Debian/Ubuntu)

```bash
sudo apt update
sudo apt install git git-lfs
git lfs install
```

For other Linux distros, refer to the [Git LFS installation guide](https://git-lfs.github.com/).

---

## ⬇️ Clone the Repository (with LFS)

Use `git` to clone this repository and fetch the LFS-tracked files:

```bash
git clone https://github.com/gutierrez-arcelus-lab/shinybcells.git
cd shinybcells
git lfs pull
```

---

## 📥 Install the Package in R

Open R or RStudio and run:

```r
# Install the package from the local directory
devtools::install_local("path/to/shinybcells")
```

---

## 🚀 Run the App

After installation, launch the app with:

```r
library(shinybcells)
bcellApp()
```

---

