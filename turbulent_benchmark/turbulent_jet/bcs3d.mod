V34 :0x24 bcs3d
26 boundary_cds_3D_module.f90 S624 0
12/01/2023  08:43:13
use vars public 0 direct
use iso_c_binding public 0 indirect
use nvf_acc_common public 0 indirect
use openacc_la public 0 direct
enduse
D 214 26 907 8 906 7
D 223 26 910 8 909 7
D 232 26 907 8 906 7
D 253 26 1003 8 1002 7
S 624 24 0 0 0 9 1 0 5013 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 bcs3d
R 906 25 7 iso_c_binding c_ptr
R 907 5 8 iso_c_binding val c_ptr
R 909 25 10 iso_c_binding c_funptr
R 910 5 11 iso_c_binding val c_funptr
R 944 6 45 iso_c_binding c_null_ptr$ac
R 946 6 47 iso_c_binding c_null_funptr$ac
R 947 26 48 iso_c_binding ==
R 949 26 50 iso_c_binding !=
R 1002 25 6 nvf_acc_common c_devptr
R 1003 5 7 nvf_acc_common cptr c_devptr
R 1009 6 13 nvf_acc_common c_null_devptr$ac
R 1047 26 51 nvf_acc_common =
S 1775 23 5 0 0 0 1776 624 14388 0 0 A 0 0 0 0 B 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_poiseuille_w_bback
S 1776 14 5 0 0 0 1 1775 14388 0 400000 A 0 0 0 0 B 0 11 0 0 0 0 0 371 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 624 0 0 0 0 bcs_poiseuille_w_bback bcs_poiseuille_w_bback 
F 1776 0
S 1777 23 5 0 0 0 1778 624 14411 0 0 A 0 0 0 0 B 0 121 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_pois_w_regularized_non_eq_extrapolation
S 1778 14 5 0 0 0 1 1777 14411 0 400000 A 0 0 0 0 B 0 121 0 0 0 0 0 372 0 0 0 0 0 0 0 0 0 0 0 0 0 121 0 624 0 0 0 0 bcs_pois_w_regularized_non_eq_extrapolation bcs_pois_w_regularized_non_eq_extrapolation 
F 1778 0
S 1779 23 5 0 0 0 1780 624 14455 0 0 A 0 0 0 0 B 0 339 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_couette_w_regularized_non_eq_extrapolation
S 1780 14 5 0 0 0 1 1779 14455 0 400000 A 0 0 0 0 B 0 339 0 0 0 0 0 373 0 0 0 0 0 0 0 0 0 0 0 0 0 339 0 624 0 0 0 0 bcs_couette_w_regularized_non_eq_extrapolation bcs_couette_w_regularized_non_eq_extrapolation 
F 1780 0
S 1781 23 5 0 0 0 1782 624 14502 0 0 A 0 0 0 0 B 0 564 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_turbulent_jet
S 1782 14 5 0 0 0 1 1781 14502 0 400000 A 0 0 0 0 B 0 564 0 0 0 0 0 374 0 0 0 0 0 0 0 0 0 0 0 0 0 564 0 624 0 0 0 0 bcs_turbulent_jet bcs_turbulent_jet 
F 1782 0
A 437 1 0 0 0 214 944 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 440 1 0 0 0 223 946 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 509 1 0 0 0 253 1009 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 133 1 1
V 437 214 7 0
S 0 214 0 0 0
A 0 6 0 0 1 2 0
J 134 1 1
V 440 223 7 0
S 0 223 0 0 0
A 0 6 0 0 1 2 0
J 36 1 1
V 509 253 7 0
S 0 253 0 0 0
A 0 232 0 0 1 437 0
Z
