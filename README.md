[![DOI](https://zenodo.org/badge/706890266.svg)](https://zenodo.org/doi/10.5281/zenodo.12798252)

# PALM_22.10_Fire
PALM 22.10 version with the temperature fix module for fire experiments.  

## Temperature overshoot/undershoot control

To control the temperature overshoot and undershoot in PALM, in PALM namelist add:
```
&temp_fix_parameters
! time start of fire is hard-coded in user_module and temp_fix module
 switch_off_module=.FALSE.,
 t_min= 291.0, ! minimum temperature allowed in the simulation
 t_max= 993.0,   ! highest prescribed fire temperature
 t_tran = 346.89, ! lowest prescribed fire temperature
/

```

## How to prescribe fire
To prescribe fire heat forcing in PALM:
### Step 1

Save an ASCII file of fire temperature with the same grid configuration as the desired simulation domain (similar to _topo input in PALM) and name it `{jobname}_fire_loc` (see example [here]([https://github.com/dongqi-DQ/PALM_22.10_Fire/blob/main/palm.iofiles](https://github.com/dongqi-DQ/PALM_22.10_Fire/blob/main/blf_flat_night_loc1_4m_fire_loc)).

### Step 2

Modify the `USER_CODE` provided in [user_module.f90](https://github.com/dongqi-DQ/PALM_22.10_Fire/blob/ddf584ac65a4892c121a9e4753fa70eccf204d54/USER_CODE/user_module.f90#L636).
```
fire_start_time = 7200.0_wp       ! time in seconds when the fire starts
fire_start_z  = 1                 ! the number of vertical level where the fire profile starts
fire_end_z    = 6                 ! the number of vertical level where the fire profile ends
fire_start_x  = 190               ! the number of grid point along x-axis where the fire profile starts
fire_end_x    = 255               ! the number of grid point along x-axis where the fire profile ends
fire_start_y  = 220               ! the number of grid point along y-xis where the fire profile starts
fire_end_y    = 255               ! the number of grid point along y-xis where the fire profile ends
fire_tt      = 993.0_wp           ! highest prescribed temperature in K
```

### Step 3
For PALM to identify the fire input file, modify `.palm.iofiles` (see example [here](https://github.com/dongqi-DQ/PALM_22.10_Fire/blob/main/palm.iofiles)):
```
FIRE_DATA                inopt:tr   d3#:d3r  $base_data/$run_identifier/INPUT          _fire_loc*
```

### Step 4

Use `compile.sh` script to compile. You may need to modify the configuration names. Then run PALM as usual. 


