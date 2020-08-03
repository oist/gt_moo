# Optimization and Validation of Diffusion MRI-based Fiber Tracking with Neural Tracer Data as a Reference

<p align="center">
  <img width="300" src="https://github.com/oist/gt_moo/blob/master/docs/moo_gt.png">
</p>

Here is the python implementation of the `data-driven framework` to optimize and validate parameters of dMRI-based fiber-tracking algorithms.
The examples use neural tracer data as reference, however, any reference data can be used (such as connectomes). To evaluate the goodness of the fiber estimation, comparisons against data references are implemented by multiple objective functions. Examples of objective functions are TPR, FPR, distance-weighted coverage, true/false positive ratio, projection coincidence, commissural passage, etc. We encourage the search of better objective functions within this framework. The frawework implements an evolutionary optimization processing several brains in parallel. A HPC (high-performance computing) compatible version is available, as well as a Desktop PC implementation.

This page is meant as a brief overview of the framework's functionality. For more detailed information please check our [preprint](https://arxiv.org/abs/1911.13215).

HPC version:
The method reported here was tested on a cluster HPC computer. In order to share champions settings over each iteration, we ran separated jobs for each brain, synchronized between them. There are 2 types of jobs, light (synchronized jobs) and heavy jobs. The light jobs keep running the evolutionary processes for the brains (1 job per brain), they were tested on a single core with low memory, however, they are active during the whole optimization process, which took around 4~5 weeks for global tracking. For an iteration, each light job creates and dispatch 8 heavy jobs (1 job per individual parameter setting). A heavy job uses more than 1 core and higher memory. They read source data (mask, tracer, dwi, atlas, injection region), perform global tracking (implemented in matlab), calculate objective functions, and store results (jobs details, parameters, fiber tracking, density maps, champions, connection matrices, objectives details) organized in a folder structure. Although the data source is not publicly available (eventually, it will be disclosed as part of Brain/MINDS project results), the code developed here can be adapted for different scenarios.

## HPC version
The method reported here was implemented on a cluster HPC computer for global-tracking algorithm. In order to share champion settings over each iteration, we ran separate jobs for each brain, and synchronized them. There are 2 types of jobs, light (synchronized) and heavy. Light jobs keep running the evolutionary processes for the brains (1 job per brain). 
For each iteration, a light job creates and dispatches 8 heavy jobs (1 job per individual parameter setting). Thus, optimization process parallelization is implemented at the level of individual brains and global-tracking runs.

To start the optimization process -> python run_moo.py

## Desktop PC veersion
An alternative portable implementation is made available for running on a desktop PC, targeting commonly used tractography approaches that not require important computing resources. This version uses mpi4py to run parallel evolutionary processes, one per brain, while tractography runs sequentially within each process. Champions sharing and processes synchronization are implemented as well. The iFOD2 algorithm case used this version of the code.


## Output
Results are organized in a folder structure.

## Running time 
The running time depends on the fiber-tracking algorithm, the complexity of the objective functions and the number of generated fibers.

## About Brain/MINDS data
Data sources used here are not publicly available (eventually, it will be disclosed as part of Brain/MINDS project results). For more details please visit [Brain/MINDS] (https://www.brainminds.riken.jp).

## Dependencies

deap

SimpleITK

scipy

json

pickle

global tracking ->  https://www.uniklinik-freiburg.de/mr-en/research-groups/diffperf/fibertools.html




Desktop PC verion:




