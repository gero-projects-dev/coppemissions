# Offshore GHG Emissions Calculator 🌍⚡

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![Code style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A Python-based tool for calculating greenhouse gas (GHG) emissions from offshore oil & gas operations, with time-vectorized computations for large-scale scenarios.

---

## Table of Contents
- [Key Features](#key-features-)
- [Methodology](#methodology-)
- [Installation](#installation-)
- [Quick Start](#quick-start-)
- [Contributing](#contributing-)
- [License](#license-)

---

## Key Features ✨
- 🚀 **Time-vectorized calculations** for efficient timeseries processing
- 🔥 Multi-source emission models:
  - Flaring & venting
  - Power generation (turbines, diesel generators)
  - Combustion (turbines, diesel generators)
  - Fugitive emissions (ISO 27916 - valves, pumps, compressors)
- ✅ Customizable emission factors (based on equipment type and fuel);
- 📊 Scenario analysis (baseline vs. mitigation strategies)

---

## Methodology 📚

Implements:

    IPCC 2006 Guidelines for stationary combustion

    OGMP 2.0 reporting standards

    Monte Carlo uncertainty analysis (optional)

---

## Installation 🐍

```bash
  # Clone repo
  git clone https://github.com/yourusername/STSW-Emissions.git
  cd STSW-Emissions
  
  # Create and activate environment (Python 3.8+ required)
  python -m venv .venv
  source .venv/bin/activate  # Linux/Mac
  .\.venv\Scripts\activate  # Windows
  
  # Install with pip
  pip install -e .
  
  # Install dependencies:
  bash
  
  pip install -r requirements.txt

```

---

## License ⚖️
**Effective Date**: [01/12/2025]  
**Licensor**: [Universidade Federal do Rio de Janeiro]  
**Licensee**: [Petrobras]

**By using this Software, Licensee agrees to these terms.**  
© [2025] [Universidade Federal do Rio de Janeiro]. All rights reserved.  
