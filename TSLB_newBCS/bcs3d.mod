V34 :0x24 bcs3d
26 boundary_cds_3D_module.f90 S624 0
02/14/2024  14:36:29
use vars public 0 direct
use iso_c_binding public 0 indirect
use nvf_acc_common public 0 indirect
use openacc_la public 0 direct
enduse
D 247 26 953 8 952 7
D 256 26 956 8 955 7
D 265 26 953 8 952 7
D 286 26 1049 8 1048 7
S 624 24 0 0 0 9 1 0 5013 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 bcs3d
R 952 25 7 iso_c_binding c_ptr
R 953 5 8 iso_c_binding val c_ptr
R 955 25 10 iso_c_binding c_funptr
R 956 5 11 iso_c_binding val c_funptr
R 990 6 45 iso_c_binding c_null_ptr$ac
R 992 6 47 iso_c_binding c_null_funptr$ac
R 993 26 48 iso_c_binding ==
R 995 26 50 iso_c_binding !=
R 1048 25 6 nvf_acc_common c_devptr
R 1049 5 7 nvf_acc_common cptr c_devptr
R 1055 6 13 nvf_acc_common c_null_devptr$ac
R 1093 26 51 nvf_acc_common =
S 1821 23 5 0 0 0 1822 624 14606 0 0 A 0 0 0 0 B 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_all_bback_2c
S 1822 14 5 0 0 0 1 1821 14606 0 400000 A 0 0 0 0 B 0 10 0 0 0 0 0 371 0 0 0 0 0 0 0 0 0 0 0 0 0 10 0 624 0 0 0 0 bcs_all_bback_2c bcs_all_bback_2c 
F 1822 0
S 1823 23 5 0 0 0 1824 624 14623 0 0 A 0 0 0 0 B 0 77 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_poiseuille_w_bback
S 1824 14 5 0 0 0 1 1823 14623 0 400000 A 0 0 0 0 B 0 77 0 0 0 0 0 372 0 0 0 0 0 0 0 0 0 0 0 0 0 77 0 624 0 0 0 0 bcs_poiseuille_w_bback bcs_poiseuille_w_bback 
F 1824 0
S 1825 23 5 0 0 0 1826 624 14646 0 0 A 0 0 0 0 B 0 187 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_pois_w_regularized_non_eq_extrapolation
S 1826 14 5 0 0 0 1 1825 14646 0 400000 A 0 0 0 0 B 0 187 0 0 0 0 0 373 0 0 0 0 0 0 0 0 0 0 0 0 0 187 0 624 0 0 0 0 bcs_pois_w_regularized_non_eq_extrapolation bcs_pois_w_regularized_non_eq_extrapolation 
F 1826 0
S 1827 23 5 0 0 0 1828 624 14690 0 0 A 0 0 0 0 B 0 405 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_couette_w_regularized_non_eq_extrapolation
S 1828 14 5 0 0 0 1 1827 14690 0 400000 A 0 0 0 0 B 0 405 0 0 0 0 0 374 0 0 0 0 0 0 0 0 0 0 0 0 0 405 0 624 0 0 0 0 bcs_couette_w_regularized_non_eq_extrapolation bcs_couette_w_regularized_non_eq_extrapolation 
F 1828 0
S 1829 23 5 0 0 0 1830 624 14737 0 0 A 0 0 0 0 B 0 630 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_turbulent_jet
S 1830 14 5 0 0 0 1 1829 14737 0 400000 A 0 0 0 0 B 0 630 0 0 0 0 0 375 0 0 0 0 0 0 0 0 0 0 0 0 0 630 0 624 0 0 0 0 bcs_turbulent_jet bcs_turbulent_jet 
F 1830 0
S 1831 23 5 0 0 0 1832 624 14755 0 0 A 0 0 0 0 B 0 862 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_2c_turbulent_jet
S 1832 14 5 0 0 0 1 1831 14755 0 400000 A 0 0 0 0 B 0 862 0 0 0 0 0 376 0 0 0 0 0 0 0 0 0 0 0 0 0 862 0 624 0 0 0 0 bcs_2c_turbulent_jet bcs_2c_turbulent_jet 
F 1832 0
S 1833 23 5 0 0 0 1834 624 14776 0 0 A 0 0 0 0 B 0 1121 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_tslb_turbojet
S 1834 14 5 0 0 0 1 1833 14776 0 400000 A 0 0 0 0 B 0 1121 0 0 0 0 0 377 0 0 0 0 0 0 0 0 0 0 0 0 0 1121 0 624 0 0 0 0 bcs_tslb_turbojet bcs_tslb_turbojet 
F 1834 0
S 1835 23 5 0 0 0 1836 624 14794 0 0 A 0 0 0 0 B 0 1233 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 bcs_tslb_only_z_turbojet
S 1836 14 5 0 0 0 1 1835 14794 0 400000 A 0 0 0 0 B 0 1233 0 0 0 0 0 378 0 0 0 0 0 0 0 0 0 0 0 0 0 1233 0 624 0 0 0 0 bcs_tslb_only_z_turbojet bcs_tslb_only_z_turbojet 
F 1836 0
S 1837 23 5 0 0 0 1838 624 14819 0 0 A 0 0 0 0 B 0 1282 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pbcs
S 1838 14 5 0 0 0 1 1837 14819 0 400000 A 0 0 0 0 B 0 1282 0 0 0 0 0 379 0 0 0 0 0 0 0 0 0 0 0 0 0 1282 0 624 0 0 0 0 pbcs pbcs 
F 1838 0
A 533 1 0 0 0 247 990 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 536 1 0 0 0 256 992 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 605 1 0 0 0 286 1055 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 133 1 1
V 533 247 7 0
S 0 247 0 0 0
A 0 6 0 0 1 2 0
J 134 1 1
V 536 256 7 0
S 0 256 0 0 0
A 0 6 0 0 1 2 0
J 36 1 1
V 605 286 7 0
S 0 286 0 0 0
A 0 265 0 0 1 533 0
Z
