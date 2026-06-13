# Rates Pricing Engine

A C++ fixed-income pricing engine for discount curve construction, swap pricing, par rate calculation, and analytical risk sensitivities.

## Features

- Cash curve construction
- Swap curve bootstrapping
- Log-linear interpolation
- Averaged quadratic interpolation
- Interest rate swap PV calculation
- Par swap rate calculation
- Rate sensitivity calculation
- CSV input/output

## Build

```bash
g++ -std=c++17 Pricing.cpp -o pricing
```

## Run

```bash
./pricing
```

The program reads from `Input.csv` and writes results to `Output.csv`.

## Documentation

See `Documentation.pdf` for the mathematical methodology, system architecture, assumptions, and validation details.
