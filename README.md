# gt_moo
Multi-objective Parameter Optimization of
DWI-based Global Fiber Tracking with Neuronal
Tracer Signal as a Reference


The method reported here was implemented on a cluster HPC computer. In order to share champions settings over each iteration, we ran separated jobs for each brain, synchronized between them. The code is publicly available on github (https://github.com/oist/gt\_moo/). There are 2 types of jobs, light (synchronized jobs) and heavy jobs. The light jobs keep running the evolutionary processes for the brains (1 job per brain), they were tested on a single core with low memory, however, they are active during the whole optimization process, which took around 4$\sim $5 weeks for global tracking. For an iteration, each light job creates and dispatch 8 heavy jobs (1 job per individual parameter setting). A heavy job uses more than 1 core and higher memory. They read source data (mask, tracer, dwi, atlas, injection region), perform global tracking (implemented in matlab), calculate objective functions, and store results (jobs details, parameters, fiber tracking, density maps, champions, connection matrices, objectives details) organized in a folder structure. Although the data source is not publicly available (eventually, it will be disclosed as part of Brain/MINDS project results), the code developed here can be adapted for different scenarios.



