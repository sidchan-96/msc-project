# MSc Project Investigating the Power Efficiency of Numerical Software

## Description

This projects explored various methods of reducing the power consumption of a system while running Conjugate Gradient Method(CGM). The Conjugate Gradient Method is commonly used for benchmarking HPC systems, although primarily used to examine the performance of a system, in this project CGM will be used to assess the power consumption of a system and how the power consumption could be reduced. 

Some of the methods that will be attempted are trying MPI, openMP and hybrid MPI and openMP to compare their effects on performance and power consumption. Using the matrix free method to reduce the memory usage and limit the data transmission from RAM to cache and between cache levels. Finally, implementing different data storage formats and comparing their effects on the power consumption.

The experiments are primarily run on the Fulhame machine, using the tx2mon package to get power readings during the benchmark.


## High Performance Conjugate Gradient

The benchmark used for the tests is High Peformance Conjugate Gradient (HPCG), a modified version of this is in this repository with some minor changes to allow using matrices with smaller dimensions in the benchmark.

## Fulhame Code

The slurm script for running HPCG on Fulhame with tx2mon reading the power consumption. Once HPCG is setup using the HPCG/INSTALL, the script can be used to run a suite of tests using various problem sizes and power readings for each of them.

## Meeting Minutes

The minutes to project meeting with my sueprvisor, Paul Bartholomew, can be found in the project wiki.