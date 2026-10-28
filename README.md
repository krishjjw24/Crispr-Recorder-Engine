# Crispr-Recorder-Engine: Hardware-Accelerated C++17 CRISPR Specificity Scorer & Cellular Lineage Recording Simulator

![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=c%2B%2B&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Build](https://img.shields.io/badge/Build-CMake%20%7C%20g%2B%2B%20--O3-success)
![License](https://img.shields.io/badge/License-MIT-blue)

A high-throughput computational bioengineering engine that combines **zero-copy C++17 genome scanning**, **nearest-neighbor RNA–DNA hybridization thermodynamics ΔG°37**, and **64-bit bitwise hardware population count (`popcount`)** off-target evaluation with a **Python stochastic CRISPR cellular lineage recording simulator**.

---

## System Architecture

```text
┌──────────────────────────┐     ┌──────────────────────────┐     ┌──────────────────────────┐     ┌──────────────────────────┐
│  1. Genomic FASTA Input  │     │  2. C++17 Bitwise Engine │     │  3. Barcode Array Sim    │     │     4. Lineage Graph     │
│ Forward & Reverse Strand │────►│ 64-Bit Integer Spacers   │────►│ Stochastic CRISPR Edits  │────►│ NetworkX Phylogenetic    │
│ Cas9 (NGG)/Cas12a (TTTV) │     │ XOR + Hardware Popcount  │     │ Accumulated Across 10    │     │ Tree Reconstruction &    │
│ 2-Bit Nucleotide Packing │     │ + Nearest-Neighbor ΔG    │     │ Cell Division Generations│     │ C++ Speed Benchmark      │
└──────────────────────────┘     └──────────────────────────┘     └──────────────────────────┘     └──────────────────────────┘
