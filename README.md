# User’s Guide for WRF-SADLES v1.0: The Simple Actuator Disk for Large Eddy Simulation in the Weather Research and Forecast Model

Hai Bui:  hai.bui@uib.no

April 2023

## I. Introduction

WRF-SADLES is the implementation of an Actuator Disk for Large Eddy Simulation in the Weather Research and Forecast Model. The purpose of WRF-SADLES is to simulate the wind turbine wakes explicitly in Large Eddy Simulation (LES) mode. There exist a few implementations of the wake model in WRF, for example, the General Actuator Disk (GAD, Mirocha, 2014). However, we have not found it publicly published. WRF-SADLES also require much less information for easier implementation and carrying out experiments.

If you use WRF-SADLES in your paper, please cite both the SADLES accompanied paper (detail will be updated in WRF-SADLES's [GitHub's repository](https://github.com/haibuihoang/WRF-SADLES) when the paper is published):

>  Hai Bui, Mostafa Bakhoday-Paskyabi, and Mohammadreza Mohammadpour-Penchah, 2023: Implementation of a Simple Actuator Disc for Large Eddy Simulation (SADLES) in the Weather Research and Forecasting model. *Submitted to Geoscientific Model Development*.

and the software as indicated in WRF-SADLES's Zenodo record.

## II. System Requirements

**Hardware Requirements:**

To run WRF-SADLES, you will need a high-performance cluster with MPI support capable of running WRF-ARW in LES mode. This is particularly necessary for realistic downscaling experiments involving multiple wind turbines/wind farms.

**Software Requirements**:

WRF-SADLES version 1.0 is built on [WRF-ARW version 4.3.1](https://github.com/wrf-model/WRF/releases/tag/v4.3.1). Trying other versions is at your own risk.

**Data Requirements:**

SADLES requires the same information as the wind turbine parameterization scheme by Fitch (2012), which includes the locations (either in longitude, latitude, or in grid indices), turbine hub heights and radii, as well as tables of the thrust coefficients and powers at different ambient wind speeds. Some studies provide such information for realistic wind turbines, such as  [Larsén and Fischereit (2021)](https://doi.org/10.5281/zenodo.4668613).

## III. Compilation

To use WRF-SADLES, follow these steps:

1. Download [WRF-ARW version 4.3.1](https://github.com/wrf-model/WRF/releases/tag/v4.3.1) and successfully test compile it in a directory, for example, at `WRFROOT=~/WRF-4.3.1`.

2. Once you have successfully test compiled WRF, copy and replace three directories in the `$WRFROOT` folder, clean the previous compilation, and recompile WRF. This can be done for either the real case `em_real` or the idealized LES case `em_les`. For example:

```bash
cp -r Registry dyn_em phys $WRFROOT
cd $WRFROOT
./clean -a
./configure
./compile em_real
```

## IV. Running WRF-SADLES

SADLES requires additional input files for turbine locations and turbine specifications that are the same as the winturbine parameterization by Fitch (2012). Please read [README.windturbine](https://github.com/wrf-model/WRF/blob/master/doc/README.windturbine) for details on the formats of the turbine input files.

### IV.1 Idealized Simulation

To begin with, compile WRF-SADLES with `./compile em_les`, which produces two executable files: `ideal.exe` for generating the idealized initial condition and `wrf.exe` for the running the WRF-SADLES. In the `example` directory, we provide the necessary files for the idealized simulation in the accompanying papers, which include:

- `input_sounding`: initial environmental conditions required for running `ideal.exe`

- `namelist_input`: required for both `ideal.exe` and `wrf.exe`

- `myoutfields.txt`: (optional), provides additional output from WRF

- `wind-turbine-6.tbl`: information on a 5-MW wind turbine with 90-m hub-height and 116-m diameter. 

- `windturbines-ij.txt`: positions of the turbines inside the domain d02,

In the example above, note the following, notice the following (check `namelist.input`):

- Two nested domain are used (`max_dom=2`)

- Simulation runs with 512 cpu cores (`nproc_x = 32`,  `nproc_y = 16`). You can change or remove those lines for a different queue job with different cores.

- SADLES is applied in the domain d02 (`sadles_opt = 0, 1`)

- Uses grid indices (`windfarm_ij = 1`) for wind turbine location

- The coriolis and roughness length can be set using `ideal_f = 0.0001177` and `ideal_znt = 0.001, 0.001`

- Cell perturbation is applied to the inflow (western) boundary of domain d02 (cell_pert_xs = 0, 1) from the surface to level 20 (`cell_pert_k1 = 20`), where it gradually decreases up to level 40 (`cell_pert_k2 = 40`). The perturbation is applied every 24 seconds (`cell_pert_interval = 72, 24`). This interval is estimated by $8\Delta x/U$, where $\Delta x$ is the domain grid size (30 m), and $U$ is the average wind speed at hub height (10 m/s).

### IV.2 Realistic Simulation

To perform a realistic simulation using the WRF-SADLES model, it is necessary to compile `em_real` normally and prepare the inputs in the same way as for WRF. The WRF user guide should be consulted for instructions on how to set up WRF. Typically, a nested domain system is used for realistic simulations, with outer domains downscaling atmospheric conditions from global data to mesoscale, and inner domains configured in LES settings to downscale mesoscales to turbine scales. The SADLES is typically applied in the innermost domain in a similar manner to the idealized simulation, but the wind turbine locations are provided by `windturbines.txt` instead of `windturbines-ij.txt`, and the namelist item `windfarm_ij` is set to `0`, which is the default value.

In the `example` directory, we provide an example of the namelist for five domains, where three outer domains (9-km, 3-km, and 1-km) are mesoscale domains, and two inner domains (200-m and 40-m) are LES domains.

### IV.3 WRF-SADLES Output Files

In addition to the standard WRF output files, WRF-SADLES produces 4 additional text output files as follows:

- `sadles_info.dxx`: where `dxx` is the domain (e.g., d02) where SADLES is applied. This file contains the number of turbines within the domain, turbine type, and turbine location in WRF indices.

- `sadles_HubSpd.dxx`: Time series of wind speed at turbine hub for each turbine. The first column is the number of hours since the start of the simulation.

- `sadles_AmbSpd.dxx`: Time series of ambient wind speed for each turbine.

- `sadles_Power.dxx`: Time series of turbine power for each turbine.

## References

- Mirocha, J., Kosovic, B., Aitken, M., and Lundquist, J.: Implementation of a generalized actuator disk wind turbine model into the weather research and forecasting model for large-eddy simulation applications, *Journal of Renewable and Sustainable Energy*, 6, 013104, 2014

- Fitch, A. C., Olson, J. B., Lundquist, J. K., Dudhia, J., Gupta, A. K., Michalakes, J., and Barstad, I.: Local and mesoscale impacts of wind  farms as parameterized in a mesoscale NWP model, *Monthly Weather Revie*w, 140, 3017–3038, 2012

- Larsén, X. G. and Fischereit, J.: A case study of wind farm effects using two wake parameterizations in WRF (V3.7.1) in the presence of low level jets, https://doi.org/10.5281/zenodo.4668613, 2021
