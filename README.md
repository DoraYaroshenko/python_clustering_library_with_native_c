# 🧠 Python Clustering Library with Native C Implementation

A high-performance clustering library that combines Python for usability with a native C implementation for efficient computation.

This project implements clustering algorithms with a focus on **performance optimization**, **modular design**, and **clean integration between Python and low-level C code**.

---

## 🚀 Overview

This library provides a Python interface for clustering algorithms whose core computations are implemented in C.

The design follows a common pattern in high-performance data processing systems:  
- Python handles usability and data processing  
- C handles computationally intensive operations  

This allows achieving significantly better performance compared to pure Python implementations.

---

## ✨ Features

- ⚡ Native C implementation for improved performance  
- 🐍 Python interface for ease of use and integration  
- 📊 Support for clustering tasks on large datasets  
- 🧩 Modular and extensible design  
- 📈 Evaluation and comparison with reference implementations  

---

## 🧠 Key Idea

Clustering algorithms often involve heavy numerical computations (distance calculations, iterations over large datasets).

To improve efficiency:
- Core algorithm logic is implemented in **C**
- Python acts as a **wrapper layer**, handling input/output and usability

This approach is widely used in production-grade libraries to balance **performance and developer experience**.

---

## 🏗️ Architecture

- Python handles:
  - Data preprocessing
  - User interface
  - Result formatting  

- C handles:
  - Core clustering computations
  - Performance-critical loops  
  - Memory-efficient operations  

---

## ⚙️ Tech Stack

- Python  
- C (native implementation)  
- NumPy  
- Pandas  
- Scikit-learn (for validation and comparison)  

---

## 📦 Installation

Clone the repository:

```bash
git clone https://github.com/DoraYaroshenko/python_clustering_library_with_native_c
cd python_clustering_library_with_native_c
```
## 📊 Evaluation

- The implementation was tested on real datasets using:
  - NumPy and Pandas for data handling
  - Scikit-learn for comparison with standard implementations

- Focus areas:
  - Performance improvements
  - Accuracy validation
  - Scalability

## 🚀 Future Improvements
- Add support for additional clustering algorithms (e.g., DBSCAN, hierarchical clustering)
- Improve memory efficiency
- Add parallelization support
- Package as a pip-installable library

## 💡 Key Takeaways
- Implemented a clustering algorithm in C for performance optimization
- Designed a Python interface for usability and flexibility
- Built a modular system combining low-level efficiency with high-level accessibility
- Gained experience in bridging Python with native code

## 👩‍💻 Author
- Dora Yaroshenko
- Geut Hadadi
