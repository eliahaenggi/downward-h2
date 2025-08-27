<img src="misc/images/fast-downward.svg" width="800" alt="Fast Downward">

# Master's Thesis: Efficient Implementation of h^2 in the Fast Downward Planning System

This repository contains the code developed as part of my Master's thesis titled "Efficient Implementation of h^2 in the Fast Downward Planning System". The thesis explores the computation and optimization of the h2 heuristic within the Fast Downward planner.

## Repository Overview

This project is based on a fork of [Fast Downward](https://github.com/aibasel/downward) (revision dated 24.11.2024). The following sections describe the main components and their locations in the codebase.

---

## h2 Heuristic Optimization

- Implementation: downward/src/heuristics/htwo_heuristic

This module includes the optimized implementation of the h2 heuristic, as described in Chapter 4 of the thesis. The original (unoptimized) implementation is still available in downward/src/heuristics/hm_heuristic.

An alternative version using binary operators can be found on the "binary_operators" branch within the same file. The basic heuristic framework of the binary operator approach resembles the existing h^max implementation.

**Usage**:
./fast-downward.py [DOMAIN_FILE] [PROBLEM_FILE] --search "astar(h2())"

---

## Pi^m Compilation and h^max Heuristic

- Task compilation: downward/src/tasks/pi_m_compiled_task
- Heuristic extension: downward/src/heuristics/hmax_heuristic

These modules implement the Pi^m compilation strategy, which evaluates h2 by computing hmax(Pi^2), as described in Chapter 6. The compilation transforms a SAS+ task into a Pi^2 STRIPS task. The transformation is optional and activated via a configuration parameter.

**Usage**:
./fast-downward.py [DOMAIN_FILE] [PROBLEM_FILE] --search "astar(hmax(pi_m_compilation=true))"

---

## STRIPS Duality h^2

- Task transformation: downward/src/tasks/dual_task
- Heuristic: downward/src/heuristics/dual_htwo_heuristic

These modules implement the regression simulation using STRIPS duality. The dual_htwo_heuristic class builds on the optimized h2 heuristic to compute h2(Pi^d), as discussed in Chapter 7.

Note: Heuristic values differ from those of the standard h2 heuristic due to the transformation of the task.

**Usage**:
./fast-downward.py [DOMAIN_FILE] [PROBLEM_FILE] --search "astar(dualh2())"

---

## Experimental Setup

- Experiment scripts: downward/experiments/hm/v1.py
- Dependencies: downward/experiments/hm/requirements.txt

This directory contains all experimental scripts and setup. The main script v1.py runs experiments using Downward Lab, a framework for managing and analyzing planning experiments.

