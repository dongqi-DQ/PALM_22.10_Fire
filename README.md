# PALM_22.10_Fire
PALM 22.10 version with the temperature fix module for fire experiments 

This documentation is under construction...

See `USER_CODE` for more.

Also in PALM namelist add:
```
&temp_fix_parameters
! time start of fire is hard-coded in user_module and temp_fix module
 switch_off_module=.FALSE.,
 t_min= 291.0, ! minimum temperature allowed in the simulation
 t_max= 993.0,   ! highest prescribed fire temperature
 t_tran = 346.89, ! lowest prescribed fire temperature
/

```
