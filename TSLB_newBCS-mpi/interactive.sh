#!/bin/bash

srun -N 1 --ntasks-per-node=1 --exclusive --gres=gpu:4 -p boost_usr_prod -t 00:20:00 -A IscrC_recTSLB --pty bash