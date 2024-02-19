#/bin/bash

mpirun --mca orte_base_help_aggregate 0 --mca btl_base_warn_component_unused 0 -np 8 ./main_mpi.x
