V26 glimmer_to_dycore
21 glimmer_to_dycore.F90 S582 0
04/19/2012  14:24:12
use glimmer_log public 0 indirect
use glimmer_sparse_type public 0 indirect
use glimmer_sparse_slap public 0 indirect
use glimmer_sparse_umfpack public 0 indirect
use glimmer_global public 0 indirect
use glimmer_sparse_pardiso public 0 indirect
use glimmer_sparse public 0 indirect
use netcdf public 0 indirect
use glimmer_ncdf public 0 indirect
use isostasy_types public 0 indirect
use profile public 0 indirect
use glimmer_coordinates public 0 indirect
use xls public 0 indirect
use remap_glamutils public 0 indirect
use glide_types public 0 direct
use glimmer_map_types public 1 61 indirect
enduse
D 63 18 32
D 65 21 63 1 3 57 0 0 0 0 0
 0 57 3 3 57 57
D 68 21 63 1 2 57 0 0 0 0 0
 2 41 3 2 41 57
D 77 24 675 280 674 7
D 101 20 6
D 103 20 6
D 105 20 9
D 173 24 697 16 696 7
D 187 24 844 184 842 7
D 205 20 9
D 207 20 6
D 209 24 857 24 856 7
D 218 20 173
D 345 24 697 16 696 7
D 353 24 844 184 842 7
D 363 24 857 24 856 7
D 369 20 345
D 371 24 920 896 919 7
D 383 24 983 776 982 7
D 395 24 1048 336 1047 7
D 401 24 1053 24 1052 7
D 416 20 353
D 418 20 371
D 420 20 383
D 2066 18 2
D 2070 24 6585 864 6584 3
D 2080 24 6599 10600 6598 0
D 2090 24 6607 11520 6606 7
D 2102 20 2090
D 2104 20 2090
D 2106 24 6626 984 6625 7
D 2124 20 6
D 2126 20 2106
D 2128 20 2106
D 2144 24 6707 176 6706 7
D 2156 24 6723 544 6722 7
D 2180 20 9
D 2182 20 9
D 2184 20 9
D 2192 24 6768 5812 6767 3
D 2305 24 7074 1880 7073 7
D 2446 24 7364 40 7363 7
D 2452 24 7385 56 7365 7
D 2461 24 7392 88 7370 7
D 2470 24 7402 80 7375 7
D 2479 24 7411 88 7380 7
D 2488 20 2452
D 2490 20 2461
D 2492 20 2470
D 2494 20 2479
D 2525 24 675 280 674 7
D 2531 20 6
D 2533 20 6
D 2535 20 9
D 2567 18 32
D 2613 18 2
D 2627 24 6607 11520 6606 7
D 2637 24 6626 984 6625 7
D 2649 24 6707 176 6706 7
D 2655 24 6723 544 6722 7
D 2661 20 9
D 2663 20 9
D 2665 20 9
D 2667 24 6768 5812 6767 3
D 2679 24 7364 40 7363 7
D 2685 24 7385 56 7365 7
D 2691 24 7392 88 7370 7
D 2697 24 7402 80 7375 7
D 2703 24 7411 88 7380 7
D 2709 20 2685
D 2711 20 2691
D 2713 20 2697
D 2715 20 2703
D 2717 24 7074 1880 7073 7
D 2723 24 7476 480 7475 7
D 2753 20 8
D 2755 20 8
D 2757 20 8
D 2759 20 8
D 2761 24 7546 320 7545 3
D 2769 24 7580 1200 7577 7
D 2835 21 6 1 3 14 0 0 0 0 0
 0 14 3 3 14 14
D 2838 20 8
D 2840 20 8
D 2842 20 9
D 2844 20 9
D 2846 20 9
D 2848 20 9
D 2850 20 9
D 2852 20 6
D 2854 20 6
D 2856 20 9
D 2858 24 7659 2912 7656 7
D 3020 20 9
D 3022 20 9
D 3024 20 9
D 3026 20 9
D 3028 20 9
D 3030 20 9
D 3032 20 9
D 3034 20 9
D 3036 20 9
D 3038 20 9
D 3040 20 9
D 3042 20 9
D 3044 20 9
D 3046 20 9
D 3048 20 9
D 3050 20 9
D 3052 20 9
D 3054 20 9
D 3056 20 9
D 3058 20 9
D 3060 20 9
D 3062 20 9
D 3064 20 9
D 3066 20 9
D 3068 20 9
D 3070 20 9
D 3072 24 7843 816 7839 7
D 3114 20 9
D 3116 20 9
D 3118 20 9
D 3120 20 9
D 3122 20 9
D 3124 20 9
D 3126 24 7892 2752 7888 7
D 3270 20 9
D 3272 20 9
D 3274 20 9
D 3276 20 9
D 3278 20 9
D 3280 20 9
D 3282 20 9
D 3284 20 9
D 3286 20 9
D 3288 20 9
D 3290 20 9
D 3292 20 9
D 3294 20 9
D 3296 20 9
D 3298 20 9
D 3300 20 9
D 3302 20 9
D 3304 20 9
D 3306 20 9
D 3308 20 9
D 3310 20 6
D 3312 20 6
D 3314 20 6
D 3316 24 8059 1176 8058 7
D 3340 20 9
D 3342 20 9
D 3344 20 9
D 3346 24 8085 920 8082 7
D 3400 20 8
D 3402 20 8
D 3404 20 8
D 3406 20 8
D 3408 20 8
D 3410 20 8
D 3412 20 8
D 3414 20 16
D 3416 24 8148 1752 8144 7
D 3512 20 9
D 3514 20 9
D 3516 20 9
D 3518 20 9
D 3520 20 9
D 3522 20 9
D 3524 20 9
D 3526 20 9
D 3528 20 9
D 3530 20 9
D 3532 20 9
D 3534 20 9
D 3536 20 9
D 3538 20 9
D 3540 20 9
D 3542 24 8262 1720 8258 7
D 3614 20 9
D 3616 20 16
D 3618 20 9
D 3620 20 9
D 3622 24 8346 416 8345 7
D 3634 20 2627
D 3636 20 2637
D 3638 24 8357 392 8356 7
D 3662 20 9
D 3664 20 9
D 3666 20 9
D 3668 24 8400 336 8397 7
D 3692 20 9
D 3694 20 9
D 3696 20 9
D 3698 24 8421 888 8419 7
D 3755 21 9 1 3 14 0 0 0 0 0
 0 14 3 3 14 14
D 3758 20 9
D 3760 20 9
D 3762 20 9
D 3764 20 9
D 3766 20 9
D 3768 20 9
D 3770 20 9
D 3772 20 9
D 3774 24 8482 544 8481 7
D 3792 21 9 1 3 14 0 0 0 0 0
 0 14 3 3 14 14
D 3795 21 9 1 3 41 0 0 0 0 0
 0 41 3 3 41 41
D 3798 20 9
D 3800 20 9
D 3802 24 8501 832 8498 7
D 3856 20 9
D 3858 20 9
D 3860 20 9
D 3862 20 9
D 3864 20 9
D 3866 20 9
D 3868 20 9
D 3870 20 9
D 3872 24 8558 2296 8554 7
D 3986 21 9 1 3 39 0 0 0 0 0
 0 39 3 3 39 39
D 3989 21 9 1 3 39 0 0 0 0 0
 0 39 3 3 39 39
D 3992 21 9 1 3 16 0 0 0 0 0
 0 16 3 3 16 16
D 3998 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 4001 20 9
D 4003 20 9
D 4005 20 9
D 4007 20 9
D 4009 20 9
D 4011 20 9
D 4013 20 9
D 4015 20 9
D 4017 20 9
D 4019 20 9
D 4021 20 9
D 4023 20 9
D 4025 20 9
D 4027 20 9
D 4029 20 9
D 4031 20 9
D 4033 20 9
D 4035 20 9
D 4037 24 8698 1344 8695 7
D 4115 20 9
D 4117 20 9
D 4119 20 9
D 4121 20 9
D 4123 20 9
D 4125 20 9
D 4127 20 9
D 4129 20 9
D 4131 20 9
D 4133 20 9
D 4135 20 9
D 4137 20 9
D 4139 24 8781 112 8780 7
D 4148 21 9 1 3 39 0 0 0 0 0
 0 39 3 3 39 39
D 4151 24 8792 664 8791 7
D 4187 20 9
D 4189 20 9
D 4191 20 9
D 4193 20 9
D 4195 20 9
D 4203 24 8851 336 8848 7
D 4227 20 9
D 4229 20 9
D 4231 20 9
D 4233 24 8871 29712 8870 7
D 4239 24 8898 30824 8897 7
D 4350 21 6 1 6469 6472 1 1 0 0 1
 3 6470 3 3 6470 6471
D 4353 21 7 1 6473 6476 1 1 0 0 1
 3 6474 3 3 6474 6475
S 582 24 0 0 0 8 1 0 4658 10005 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 glimmer_to_dycore
S 585 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 586 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 592 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 610 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 611 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 612 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 613 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 614 3 0 0 0 2567 1 1 0 0 0 A 0 0 0 0 0 0 0 0 4833 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 15 2a 20 55 4e 4b 4e 4f 57 4e 20 20 20 20 20 20
S 615 3 0 0 0 2567 1 1 0 0 0 A 0 0 0 0 0 0 0 0 4849 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 15 2a 20 20 20 20 20 20 20 20 20 20 20 20 20 20
S 616 3 0 0 0 2567 1 1 0 0 0 A 0 0 0 0 0 0 0 0 4865 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 15 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20
S 617 3 0 0 0 2567 1 1 0 0 0 A 0 0 0 0 0 0 0 0 4881 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 15 2a 20 57 41 52 4e 49 4e 47 3a 20 20 20 20 20
S 618 3 0 0 0 2567 1 1 0 0 0 A 0 0 0 0 0 0 0 0 4897 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 15 2a 20 45 52 52 4f 52 3a 20 20 20 20 20 20 20
S 619 3 0 0 0 2567 1 1 0 0 0 A 0 0 0 0 0 0 0 0 4913 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 15 2a 20 46 41 54 41 4c 20 45 52 52 4f 52 20 3a
S 620 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 632 7 12 glimmer_log msg_prefix$ac
S 668 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 672 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 674 25 2 glimmer_sparse_type sparse_matrix_type
R 675 5 3 glimmer_sparse_type nonzeros sparse_matrix_type
R 676 5 4 glimmer_sparse_type order sparse_matrix_type
R 677 5 5 glimmer_sparse_type symmetric sparse_matrix_type
R 679 5 7 glimmer_sparse_type col sparse_matrix_type
R 680 5 8 glimmer_sparse_type col$sd sparse_matrix_type
R 681 5 9 glimmer_sparse_type col$p sparse_matrix_type
R 682 5 10 glimmer_sparse_type col$o sparse_matrix_type
R 685 5 13 glimmer_sparse_type row sparse_matrix_type
R 686 5 14 glimmer_sparse_type row$sd sparse_matrix_type
R 687 5 15 glimmer_sparse_type row$p sparse_matrix_type
R 688 5 16 glimmer_sparse_type row$o sparse_matrix_type
R 691 5 19 glimmer_sparse_type val sparse_matrix_type
R 692 5 20 glimmer_sparse_type val$sd sparse_matrix_type
R 693 5 21 glimmer_sparse_type val$p sparse_matrix_type
R 694 5 22 glimmer_sparse_type val$o sparse_matrix_type
R 696 25 24 glimmer_sparse_type sparse_solver_options_base
R 697 5 25 glimmer_sparse_type tolerance sparse_solver_options_base
R 698 5 26 glimmer_sparse_type maxiters sparse_solver_options_base
R 699 5 27 glimmer_sparse_type method sparse_solver_options_base
R 842 25 3 glimmer_sparse_slap slap_solver_workspace
R 844 5 5 glimmer_sparse_slap rwork slap_solver_workspace
R 845 5 6 glimmer_sparse_slap rwork$sd slap_solver_workspace
R 846 5 7 glimmer_sparse_slap rwork$p slap_solver_workspace
R 847 5 8 glimmer_sparse_slap rwork$o slap_solver_workspace
R 850 5 11 glimmer_sparse_slap iwork slap_solver_workspace
R 851 5 12 glimmer_sparse_slap iwork$sd slap_solver_workspace
R 852 5 13 glimmer_sparse_slap iwork$p slap_solver_workspace
R 853 5 14 glimmer_sparse_slap iwork$o slap_solver_workspace
R 855 5 16 glimmer_sparse_slap max_nelt slap_solver_workspace
R 856 25 17 glimmer_sparse_slap slap_solver_options
R 857 5 18 glimmer_sparse_slap itol slap_solver_options
R 858 5 19 glimmer_sparse_slap use_gmres slap_solver_options
R 859 5 20 glimmer_sparse_slap gmres_saved_vectors slap_solver_options
R 860 5 21 glimmer_sparse_slap base slap_solver_options
R 862 5 23 glimmer_sparse_slap base$p slap_solver_options
S 916 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 919 25 2 glimmer_sparse_umfpack umf_solver_workspace
R 920 5 3 glimmer_sparse_umfpack control umf_solver_workspace
R 921 5 4 glimmer_sparse_umfpack info umf_solver_workspace
R 922 5 5 glimmer_sparse_umfpack numeric umf_solver_workspace
R 923 5 6 glimmer_sparse_umfpack symbolic umf_solver_workspace
R 924 5 7 glimmer_sparse_umfpack alloc umf_solver_workspace
R 982 25 1 glimmer_sparse_pardiso pardiso_solver_workspace
R 983 5 2 glimmer_sparse_pardiso pt pardiso_solver_workspace
R 984 5 3 glimmer_sparse_pardiso error pardiso_solver_workspace
R 985 5 4 glimmer_sparse_pardiso dparm pardiso_solver_workspace
R 1047 25 3 glimmer_sparse sparse_solver_options
R 1048 5 4 glimmer_sparse base sparse_solver_options
R 1049 5 5 glimmer_sparse slap sparse_solver_options
R 1050 5 6 glimmer_sparse umf sparse_solver_options
R 1051 5 7 glimmer_sparse pardiso sparse_solver_options
R 1052 25 8 glimmer_sparse sparse_solver_workspace
R 1053 5 9 glimmer_sparse slap sparse_solver_workspace
R 1055 5 11 glimmer_sparse slap$p sparse_solver_workspace
R 1057 5 13 glimmer_sparse umf sparse_solver_workspace
R 1059 5 15 glimmer_sparse umf$p sparse_solver_workspace
R 1061 5 17 glimmer_sparse pardiso sparse_solver_workspace
R 1063 5 19 glimmer_sparse pardiso$p sparse_solver_workspace
S 1172 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 6571 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 6573 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1343554297 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 6574 3 0 0 0 16 1 1 0 0 0 A 0 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 6575 3 0 0 0 16 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
S 6576 3 0 0 0 20 1 1 0 0 0 A 0 0 0 0 0 0 0 0 23916 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 1 20
S 6578 3 0 0 0 2613 1 1 0 0 0 A 0 0 0 0 0 0 0 0 23918 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 0
R 6584 25 5 glimmer_ncdf glimmer_nc_stat
R 6585 5 6 glimmer_ncdf define_mode glimmer_nc_stat
R 6586 5 7 glimmer_ncdf just_processed glimmer_nc_stat
R 6587 5 8 glimmer_ncdf processsed_time glimmer_nc_stat
R 6588 5 9 glimmer_ncdf filename glimmer_nc_stat
R 6589 5 10 glimmer_ncdf id glimmer_nc_stat
R 6590 5 11 glimmer_ncdf nlevel glimmer_nc_stat
R 6591 5 12 glimmer_ncdf nstaglevel glimmer_nc_stat
R 6592 5 13 glimmer_ncdf nstagwbndlevel glimmer_nc_stat
R 6593 5 14 glimmer_ncdf timedim glimmer_nc_stat
R 6594 5 15 glimmer_ncdf timevar glimmer_nc_stat
R 6595 5 16 glimmer_ncdf vars glimmer_nc_stat
R 6596 5 17 glimmer_ncdf hotstart glimmer_nc_stat
R 6597 5 18 glimmer_ncdf vars_copy glimmer_nc_stat
R 6598 25 19 glimmer_ncdf glimmer_nc_meta
R 6599 5 20 glimmer_ncdf title glimmer_nc_meta
R 6600 5 21 glimmer_ncdf institution glimmer_nc_meta
R 6601 5 22 glimmer_ncdf references glimmer_nc_meta
R 6602 5 23 glimmer_ncdf source glimmer_nc_meta
R 6603 5 24 glimmer_ncdf history glimmer_nc_meta
R 6604 5 25 glimmer_ncdf comment glimmer_nc_meta
R 6605 5 26 glimmer_ncdf config glimmer_nc_meta
R 6606 25 27 glimmer_ncdf glimmer_nc_output
R 6607 5 28 glimmer_ncdf nc glimmer_nc_output
R 6608 5 29 glimmer_ncdf freq glimmer_nc_output
R 6609 5 30 glimmer_ncdf next_write glimmer_nc_output
R 6610 5 31 glimmer_ncdf end_write glimmer_nc_output
R 6611 5 32 glimmer_ncdf timecounter glimmer_nc_output
R 6612 5 33 glimmer_ncdf total_time glimmer_nc_output
R 6613 5 34 glimmer_ncdf default_xtype glimmer_nc_output
R 6614 5 35 glimmer_ncdf do_averages glimmer_nc_output
R 6615 5 36 glimmer_ncdf metadata glimmer_nc_output
R 6616 5 37 glimmer_ncdf next glimmer_nc_output
R 6618 5 39 glimmer_ncdf next$p glimmer_nc_output
R 6620 5 41 glimmer_ncdf previous glimmer_nc_output
R 6622 5 43 glimmer_ncdf previous$p glimmer_nc_output
R 6624 5 45 glimmer_ncdf append glimmer_nc_output
R 6625 25 46 glimmer_ncdf glimmer_nc_input
R 6626 5 47 glimmer_ncdf nc glimmer_nc_input
R 6628 5 49 glimmer_ncdf times glimmer_nc_input
R 6629 5 50 glimmer_ncdf times$sd glimmer_nc_input
R 6630 5 51 glimmer_ncdf times$p glimmer_nc_input
R 6631 5 52 glimmer_ncdf times$o glimmer_nc_input
R 6633 5 54 glimmer_ncdf nt glimmer_nc_input
R 6634 5 55 glimmer_ncdf current_time glimmer_nc_input
R 6635 5 56 glimmer_ncdf get_time_slice glimmer_nc_input
R 6636 5 57 glimmer_ncdf next glimmer_nc_input
R 6638 5 59 glimmer_ncdf next$p glimmer_nc_input
R 6640 5 61 glimmer_ncdf previous glimmer_nc_input
R 6642 5 63 glimmer_ncdf previous$p glimmer_nc_input
S 6699 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1744706593 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 6700 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 6703 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1165623296 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 6704 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1140457472 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
R 6706 25 2 isostasy_types isostasy_elastic
R 6707 5 3 isostasy_types d isostasy_elastic
R 6708 5 4 isostasy_types lr isostasy_elastic
R 6709 5 5 isostasy_types a isostasy_elastic
R 6710 5 6 isostasy_types c1 isostasy_elastic
R 6711 5 7 isostasy_types c2 isostasy_elastic
R 6712 5 8 isostasy_types cd3 isostasy_elastic
R 6713 5 9 isostasy_types cd4 isostasy_elastic
R 6716 5 12 isostasy_types w isostasy_elastic
R 6717 5 13 isostasy_types w$sd isostasy_elastic
R 6718 5 14 isostasy_types w$p isostasy_elastic
R 6719 5 15 isostasy_types w$o isostasy_elastic
R 6721 5 17 isostasy_types wsize isostasy_elastic
R 6722 25 18 isostasy_types isos_type
R 6723 5 19 isostasy_types do_isos isos_type
R 6724 5 20 isostasy_types lithosphere isos_type
R 6725 5 21 isostasy_types asthenosphere isos_type
R 6726 5 22 isostasy_types relaxed_tau isos_type
R 6727 5 23 isostasy_types period isos_type
R 6728 5 24 isostasy_types next_calc isos_type
R 6729 5 25 isostasy_types new_load isos_type
R 6730 5 26 isostasy_types rbel isos_type
R 6733 5 29 isostasy_types relx isos_type
R 6734 5 30 isostasy_types relx$sd isos_type
R 6735 5 31 isostasy_types relx$p isos_type
R 6736 5 32 isostasy_types relx$o isos_type
R 6740 5 36 isostasy_types load isos_type
R 6741 5 37 isostasy_types load$sd isos_type
R 6742 5 38 isostasy_types load$p isos_type
R 6743 5 39 isostasy_types load$o isos_type
R 6747 5 43 isostasy_types load_factors isos_type
R 6748 5 44 isostasy_types load_factors$sd isos_type
R 6749 5 45 isostasy_types load_factors$p isos_type
R 6750 5 46 isostasy_types load_factors$o isos_type
R 6767 25 3 profile profile_type
R 6768 5 4 profile profile_unit profile_type
R 6769 5 5 profile start_time profile_type
R 6770 5 6 profile nump profile_type
R 6771 5 7 profile pstart profile_type
R 6772 5 8 profile ptotal profile_type
R 6773 5 9 profile pname profile_type
S 6813 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 6994 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1084868608 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 6996 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1083129856 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
R 7073 25 9 remap_glamutils remap_glamutils_workspace
R 7074 5 10 remap_glamutils ewn_ir remap_glamutils_workspace
R 7075 5 11 remap_glamutils nsn_ir remap_glamutils_workspace
R 7076 5 12 remap_glamutils upn_ir remap_glamutils_workspace
R 7080 5 16 remap_glamutils thck_ir remap_glamutils_workspace
R 7081 5 17 remap_glamutils thck_ir$sd remap_glamutils_workspace
R 7082 5 18 remap_glamutils thck_ir$p remap_glamutils_workspace
R 7083 5 19 remap_glamutils thck_ir$o remap_glamutils_workspace
R 7085 5 21 remap_glamutils dew_ir remap_glamutils_workspace
R 7089 5 25 remap_glamutils dew_ir$sd remap_glamutils_workspace
R 7090 5 26 remap_glamutils dew_ir$p remap_glamutils_workspace
R 7091 5 27 remap_glamutils dew_ir$o remap_glamutils_workspace
R 7093 5 29 remap_glamutils dns_ir remap_glamutils_workspace
R 7097 5 33 remap_glamutils dns_ir$sd remap_glamutils_workspace
R 7098 5 34 remap_glamutils dns_ir$p remap_glamutils_workspace
R 7099 5 35 remap_glamutils dns_ir$o remap_glamutils_workspace
R 7101 5 37 remap_glamutils dewt_ir remap_glamutils_workspace
R 7105 5 41 remap_glamutils dewt_ir$sd remap_glamutils_workspace
R 7106 5 42 remap_glamutils dewt_ir$p remap_glamutils_workspace
R 7107 5 43 remap_glamutils dewt_ir$o remap_glamutils_workspace
R 7109 5 45 remap_glamutils dnst_ir remap_glamutils_workspace
R 7113 5 49 remap_glamutils dnst_ir$sd remap_glamutils_workspace
R 7114 5 50 remap_glamutils dnst_ir$p remap_glamutils_workspace
R 7115 5 51 remap_glamutils dnst_ir$o remap_glamutils_workspace
R 7117 5 53 remap_glamutils dewu_ir remap_glamutils_workspace
R 7121 5 57 remap_glamutils dewu_ir$sd remap_glamutils_workspace
R 7122 5 58 remap_glamutils dewu_ir$p remap_glamutils_workspace
R 7123 5 59 remap_glamutils dewu_ir$o remap_glamutils_workspace
R 7125 5 61 remap_glamutils dnsu_ir remap_glamutils_workspace
R 7129 5 65 remap_glamutils dnsu_ir$sd remap_glamutils_workspace
R 7130 5 66 remap_glamutils dnsu_ir$p remap_glamutils_workspace
R 7131 5 67 remap_glamutils dnsu_ir$o remap_glamutils_workspace
R 7133 5 69 remap_glamutils hm_ir remap_glamutils_workspace
R 7137 5 73 remap_glamutils hm_ir$sd remap_glamutils_workspace
R 7138 5 74 remap_glamutils hm_ir$p remap_glamutils_workspace
R 7139 5 75 remap_glamutils hm_ir$o remap_glamutils_workspace
R 7141 5 77 remap_glamutils tarear_ir remap_glamutils_workspace
R 7145 5 81 remap_glamutils tarear_ir$sd remap_glamutils_workspace
R 7146 5 82 remap_glamutils tarear_ir$p remap_glamutils_workspace
R 7147 5 83 remap_glamutils tarear_ir$o remap_glamutils_workspace
R 7149 5 85 remap_glamutils uvel_ir remap_glamutils_workspace
R 7153 5 89 remap_glamutils uvel_ir$sd remap_glamutils_workspace
R 7154 5 90 remap_glamutils uvel_ir$p remap_glamutils_workspace
R 7155 5 91 remap_glamutils uvel_ir$o remap_glamutils_workspace
R 7157 5 93 remap_glamutils vvel_ir remap_glamutils_workspace
R 7161 5 97 remap_glamutils vvel_ir$sd remap_glamutils_workspace
R 7162 5 98 remap_glamutils vvel_ir$p remap_glamutils_workspace
R 7163 5 99 remap_glamutils vvel_ir$o remap_glamutils_workspace
R 7169 5 105 remap_glamutils trace_ir remap_glamutils_workspace
R 7170 5 106 remap_glamutils trace_ir$sd remap_glamutils_workspace
R 7171 5 107 remap_glamutils trace_ir$p remap_glamutils_workspace
R 7172 5 108 remap_glamutils trace_ir$o remap_glamutils_workspace
R 7175 5 111 remap_glamutils dsigma_ir remap_glamutils_workspace
R 7176 5 112 remap_glamutils dsigma_ir$sd remap_glamutils_workspace
R 7177 5 113 remap_glamutils dsigma_ir$p remap_glamutils_workspace
R 7178 5 114 remap_glamutils dsigma_ir$o remap_glamutils_workspace
R 7180 5 116 remap_glamutils dt_ir remap_glamutils_workspace
R 7183 5 119 remap_glamutils mask_ir remap_glamutils_workspace
R 7184 5 120 remap_glamutils mask_ir$sd remap_glamutils_workspace
R 7185 5 121 remap_glamutils mask_ir$p remap_glamutils_workspace
R 7186 5 122 remap_glamutils mask_ir$o remap_glamutils_workspace
R 7363 25 2 glimmer_map_types glimmap_proj
R 7364 5 3 glimmer_map_types found glimmap_proj
R 7365 25 4 glimmer_map_types proj_laea
R 7366 5 5 glimmer_map_types laea glimmap_proj
R 7368 5 7 glimmer_map_types laea$p glimmap_proj
R 7370 25 9 glimmer_map_types proj_aea
R 7371 5 10 glimmer_map_types aea glimmap_proj
R 7373 5 12 glimmer_map_types aea$p glimmap_proj
R 7375 25 14 glimmer_map_types proj_lcc
R 7376 5 15 glimmer_map_types lcc glimmap_proj
R 7378 5 17 glimmer_map_types lcc$p glimmap_proj
R 7380 25 19 glimmer_map_types proj_stere
R 7381 5 20 glimmer_map_types stere glimmap_proj
R 7383 5 22 glimmer_map_types stere$p glimmap_proj
R 7385 5 24 glimmer_map_types longitude_of_central_meridian proj_laea
R 7386 5 25 glimmer_map_types latitude_of_projection_origin proj_laea
R 7387 5 26 glimmer_map_types false_easting proj_laea
R 7388 5 27 glimmer_map_types false_northing proj_laea
R 7389 5 28 glimmer_map_types sinp proj_laea
R 7390 5 29 glimmer_map_types cosp proj_laea
R 7391 5 30 glimmer_map_types pole proj_laea
R 7392 5 31 glimmer_map_types standard_parallel proj_aea
R 7393 5 32 glimmer_map_types longitude_of_central_meridian proj_aea
R 7394 5 33 glimmer_map_types latitude_of_projection_origin proj_aea
R 7395 5 34 glimmer_map_types false_easting proj_aea
R 7396 5 35 glimmer_map_types false_northing proj_aea
R 7397 5 36 glimmer_map_types rho0 proj_aea
R 7398 5 37 glimmer_map_types rho0_r proj_aea
R 7399 5 38 glimmer_map_types c proj_aea
R 7400 5 39 glimmer_map_types n proj_aea
R 7401 5 40 glimmer_map_types i_n proj_aea
R 7402 5 41 glimmer_map_types standard_parallel proj_lcc
R 7403 5 42 glimmer_map_types longitude_of_central_meridian proj_lcc
R 7404 5 43 glimmer_map_types latitude_of_projection_origin proj_lcc
R 7405 5 44 glimmer_map_types false_easting proj_lcc
R 7406 5 45 glimmer_map_types false_northing proj_lcc
R 7407 5 46 glimmer_map_types rho0 proj_lcc
R 7408 5 47 glimmer_map_types f proj_lcc
R 7409 5 48 glimmer_map_types n proj_lcc
R 7410 5 49 glimmer_map_types i_n proj_lcc
R 7411 5 50 glimmer_map_types longitude_of_central_meridian proj_stere
R 7412 5 51 glimmer_map_types latitude_of_projection_origin proj_stere
R 7413 5 52 glimmer_map_types scale_factor_at_proj_origin proj_stere
R 7414 5 53 glimmer_map_types standard_parallel proj_stere
R 7415 5 54 glimmer_map_types false_easting proj_stere
R 7416 5 55 glimmer_map_types false_northing proj_stere
R 7417 5 56 glimmer_map_types pole proj_stere
R 7418 5 57 glimmer_map_types equatorial proj_stere
R 7419 5 58 glimmer_map_types k0 proj_stere
R 7420 5 59 glimmer_map_types ik0 proj_stere
R 7421 5 60 glimmer_map_types sinp proj_stere
R 7422 5 61 glimmer_map_types cosp proj_stere
S 7450 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1064011039 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 7451 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1061997773 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 7452 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 -1 -979615744 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 7453 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1074423398 1717986918 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7454 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1184645120 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 7455 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1101004800 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 7456 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1120403456 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
S 7457 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 -1066860544 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7458 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1072273817 -1717986918 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7459 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1087604736 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7460 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 9999999 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 7461 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1074266112 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7462 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1076101120 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7463 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1070176665 -1717986918 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7464 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1065646817 1202590843 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7465 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 -1079404135 -1717986918 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7466 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1016910514 -1747416644 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7467 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1071434956 -858993459 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7468 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1072064102 1717986918 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7469 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1069463633 -343597384 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7470 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1044740494 -500134854 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7471 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1037794527 -640172613 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7472 3 0 0 0 9 1 1 0 0 0 A 0 0 0 0 0 0 0 0 1103994782 1073741824 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 7473 3 0 0 0 8 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 1101896090 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8
R 7475 25 2 glide_types glide_general
R 7476 5 3 glide_types ewn glide_general
R 7477 5 4 glide_types nsn glide_general
R 7478 5 5 glide_types upn glide_general
R 7479 5 6 glide_types ice_grid glide_general
R 7480 5 7 glide_types velo_grid glide_general
R 7482 5 9 glide_types x0 glide_general
R 7483 5 10 glide_types x0$sd glide_general
R 7484 5 11 glide_types x0$p glide_general
R 7485 5 12 glide_types x0$o glide_general
R 7488 5 15 glide_types y0 glide_general
R 7489 5 16 glide_types y0$sd glide_general
R 7490 5 17 glide_types y0$p glide_general
R 7491 5 18 glide_types y0$o glide_general
R 7494 5 21 glide_types x1 glide_general
R 7495 5 22 glide_types x1$sd glide_general
R 7496 5 23 glide_types x1$p glide_general
R 7497 5 24 glide_types x1$o glide_general
R 7500 5 27 glide_types y1 glide_general
R 7501 5 28 glide_types y1$sd glide_general
R 7502 5 29 glide_types y1$p glide_general
R 7503 5 30 glide_types y1$o glide_general
R 7545 25 72 glide_types glide_options
R 7546 5 73 glide_types whichtemp glide_options
R 7547 5 74 glide_types whichflwa glide_options
R 7548 5 75 glide_types which_bmod glide_options
R 7549 5 76 glide_types whichbwat glide_options
R 7550 5 77 glide_types whichmarn glide_options
R 7551 5 78 glide_types whichbtrc glide_options
R 7552 5 79 glide_types whichevol glide_options
R 7553 5 80 glide_types whichwvel glide_options
R 7554 5 81 glide_types whichrelaxed glide_options
R 7555 5 82 glide_types hotstart glide_options
R 7556 5 83 glide_types which_ho_diagnostic glide_options
R 7557 5 84 glide_types which_ho_prognostic glide_options
R 7558 5 85 glide_types which_ho_babc glide_options
R 7559 5 86 glide_types which_ho_efvs glide_options
R 7560 5 87 glide_types which_ho_resid glide_options
R 7561 5 88 glide_types which_disp glide_options
R 7562 5 89 glide_types which_bmelt glide_options
R 7563 5 90 glide_types which_ho_nonlinear glide_options
R 7564 5 91 glide_types which_ho_sparse glide_options
R 7565 5 92 glide_types which_ho_sparse_fallback glide_options
R 7566 5 93 glide_types which_ho_source glide_options
R 7567 5 94 glide_types ho_include_thinice glide_options
R 7568 5 95 glide_types periodic_ew glide_options
R 7569 5 96 glide_types periodic_ns glide_options
R 7570 5 97 glide_types gthf glide_options
R 7571 5 98 glide_types which_sigma glide_options
R 7572 5 99 glide_types which_sigma_builtin glide_options
R 7573 5 100 glide_types diagnostic_run glide_options
R 7574 5 101 glide_types basal_mbal glide_options
R 7575 5 102 glide_types external_dycore_type glide_options
R 7576 5 103 glide_types dycore_input_file glide_options
R 7577 25 104 glide_types glide_geometry
R 7580 5 107 glide_types temporary0 glide_geometry
R 7581 5 108 glide_types temporary0$sd glide_geometry
R 7582 5 109 glide_types temporary0$p glide_geometry
R 7583 5 110 glide_types temporary0$o glide_geometry
R 7587 5 114 glide_types temporary1 glide_geometry
R 7588 5 115 glide_types temporary1$sd glide_geometry
R 7589 5 116 glide_types temporary1$p glide_geometry
R 7590 5 117 glide_types temporary1$o glide_geometry
R 7594 5 121 glide_types thck glide_geometry
R 7595 5 122 glide_types thck$sd glide_geometry
R 7596 5 123 glide_types thck$p glide_geometry
R 7597 5 124 glide_types thck$o glide_geometry
R 7601 5 128 glide_types usrf glide_geometry
R 7602 5 129 glide_types usrf$sd glide_geometry
R 7603 5 130 glide_types usrf$p glide_geometry
R 7604 5 131 glide_types usrf$o glide_geometry
R 7608 5 135 glide_types lsrf glide_geometry
R 7609 5 136 glide_types lsrf$sd glide_geometry
R 7610 5 137 glide_types lsrf$p glide_geometry
R 7611 5 138 glide_types lsrf$o glide_geometry
R 7615 5 142 glide_types topg glide_geometry
R 7616 5 143 glide_types topg$sd glide_geometry
R 7617 5 144 glide_types topg$p glide_geometry
R 7618 5 145 glide_types topg$o glide_geometry
R 7623 5 150 glide_types age glide_geometry
R 7624 5 151 glide_types age$sd glide_geometry
R 7625 5 152 glide_types age$p glide_geometry
R 7626 5 153 glide_types age$o glide_geometry
R 7630 5 157 glide_types mask glide_geometry
R 7631 5 158 glide_types mask$sd glide_geometry
R 7632 5 159 glide_types mask$p glide_geometry
R 7633 5 160 glide_types mask$o glide_geometry
R 7637 5 164 glide_types thkmask glide_geometry
R 7638 5 165 glide_types thkmask$sd glide_geometry
R 7639 5 166 glide_types thkmask$p glide_geometry
R 7640 5 167 glide_types thkmask$o glide_geometry
R 7644 5 171 glide_types marine_bc_normal glide_geometry
R 7645 5 172 glide_types marine_bc_normal$sd glide_geometry
R 7646 5 173 glide_types marine_bc_normal$p glide_geometry
R 7647 5 174 glide_types marine_bc_normal$o glide_geometry
R 7649 5 176 glide_types totpts glide_geometry
R 7650 5 177 glide_types dom glide_geometry
R 7651 5 178 glide_types empty glide_geometry
R 7652 5 179 glide_types ivol glide_geometry
R 7653 5 180 glide_types iarea glide_geometry
R 7654 5 181 glide_types iareag glide_geometry
R 7655 5 182 glide_types iareaf glide_geometry
R 7656 25 183 glide_types glide_geomderv
R 7659 5 186 glide_types dthckdew glide_geomderv
R 7660 5 187 glide_types dthckdew$sd glide_geomderv
R 7661 5 188 glide_types dthckdew$p glide_geomderv
R 7662 5 189 glide_types dthckdew$o glide_geomderv
R 7666 5 193 glide_types dusrfdew glide_geomderv
R 7667 5 194 glide_types dusrfdew$sd glide_geomderv
R 7668 5 195 glide_types dusrfdew$p glide_geomderv
R 7669 5 196 glide_types dusrfdew$o glide_geomderv
R 7673 5 200 glide_types dthckdns glide_geomderv
R 7674 5 201 glide_types dthckdns$sd glide_geomderv
R 7675 5 202 glide_types dthckdns$p glide_geomderv
R 7676 5 203 glide_types dthckdns$o glide_geomderv
R 7680 5 207 glide_types dusrfdns glide_geomderv
R 7681 5 208 glide_types dusrfdns$sd glide_geomderv
R 7682 5 209 glide_types dusrfdns$p glide_geomderv
R 7683 5 210 glide_types dusrfdns$o glide_geomderv
R 7687 5 214 glide_types dlsrfdew glide_geomderv
R 7688 5 215 glide_types dlsrfdew$sd glide_geomderv
R 7689 5 216 glide_types dlsrfdew$p glide_geomderv
R 7690 5 217 glide_types dlsrfdew$o glide_geomderv
R 7694 5 221 glide_types dlsrfdns glide_geomderv
R 7695 5 222 glide_types dlsrfdns$sd glide_geomderv
R 7696 5 223 glide_types dlsrfdns$p glide_geomderv
R 7697 5 224 glide_types dlsrfdns$o glide_geomderv
R 7701 5 228 glide_types d2usrfdew2 glide_geomderv
R 7702 5 229 glide_types d2usrfdew2$sd glide_geomderv
R 7703 5 230 glide_types d2usrfdew2$p glide_geomderv
R 7704 5 231 glide_types d2usrfdew2$o glide_geomderv
R 7708 5 235 glide_types d2usrfdns2 glide_geomderv
R 7709 5 236 glide_types d2usrfdns2$sd glide_geomderv
R 7710 5 237 glide_types d2usrfdns2$p glide_geomderv
R 7711 5 238 glide_types d2usrfdns2$o glide_geomderv
R 7715 5 242 glide_types d2thckdew2 glide_geomderv
R 7716 5 243 glide_types d2thckdew2$sd glide_geomderv
R 7717 5 244 glide_types d2thckdew2$p glide_geomderv
R 7718 5 245 glide_types d2thckdew2$o glide_geomderv
R 7722 5 249 glide_types d2thckdns2 glide_geomderv
R 7723 5 250 glide_types d2thckdns2$sd glide_geomderv
R 7724 5 251 glide_types d2thckdns2$p glide_geomderv
R 7725 5 252 glide_types d2thckdns2$o glide_geomderv
R 7729 5 256 glide_types dthckdew_unstag glide_geomderv
R 7730 5 257 glide_types dthckdew_unstag$sd glide_geomderv
R 7731 5 258 glide_types dthckdew_unstag$p glide_geomderv
R 7732 5 259 glide_types dthckdew_unstag$o glide_geomderv
R 7736 5 263 glide_types dusrfdew_unstag glide_geomderv
R 7737 5 264 glide_types dusrfdew_unstag$sd glide_geomderv
R 7738 5 265 glide_types dusrfdew_unstag$p glide_geomderv
R 7739 5 266 glide_types dusrfdew_unstag$o glide_geomderv
R 7743 5 270 glide_types dthckdns_unstag glide_geomderv
R 7744 5 271 glide_types dthckdns_unstag$sd glide_geomderv
R 7745 5 272 glide_types dthckdns_unstag$p glide_geomderv
R 7746 5 273 glide_types dthckdns_unstag$o glide_geomderv
R 7750 5 277 glide_types dusrfdns_unstag glide_geomderv
R 7751 5 278 glide_types dusrfdns_unstag$sd glide_geomderv
R 7752 5 279 glide_types dusrfdns_unstag$p glide_geomderv
R 7753 5 280 glide_types dusrfdns_unstag$o glide_geomderv
R 7757 5 284 glide_types dlsrfdew_unstag glide_geomderv
R 7758 5 285 glide_types dlsrfdew_unstag$sd glide_geomderv
R 7759 5 286 glide_types dlsrfdew_unstag$p glide_geomderv
R 7760 5 287 glide_types dlsrfdew_unstag$o glide_geomderv
R 7764 5 291 glide_types dlsrfdns_unstag glide_geomderv
R 7765 5 292 glide_types dlsrfdns_unstag$sd glide_geomderv
R 7766 5 293 glide_types dlsrfdns_unstag$p glide_geomderv
R 7767 5 294 glide_types dlsrfdns_unstag$o glide_geomderv
R 7771 5 298 glide_types d2usrfdew2_unstag glide_geomderv
R 7772 5 299 glide_types d2usrfdew2_unstag$sd glide_geomderv
R 7773 5 300 glide_types d2usrfdew2_unstag$p glide_geomderv
R 7774 5 301 glide_types d2usrfdew2_unstag$o glide_geomderv
R 7778 5 305 glide_types d2usrfdns2_unstag glide_geomderv
R 7779 5 306 glide_types d2usrfdns2_unstag$sd glide_geomderv
R 7780 5 307 glide_types d2usrfdns2_unstag$p glide_geomderv
R 7781 5 308 glide_types d2usrfdns2_unstag$o glide_geomderv
R 7785 5 312 glide_types d2thckdew2_unstag glide_geomderv
R 7786 5 313 glide_types d2thckdew2_unstag$sd glide_geomderv
R 7787 5 314 glide_types d2thckdew2_unstag$p glide_geomderv
R 7788 5 315 glide_types d2thckdew2_unstag$o glide_geomderv
R 7792 5 319 glide_types d2thckdns2_unstag glide_geomderv
R 7793 5 320 glide_types d2thckdns2_unstag$sd glide_geomderv
R 7794 5 321 glide_types d2thckdns2_unstag$p glide_geomderv
R 7795 5 322 glide_types d2thckdns2_unstag$o glide_geomderv
R 7799 5 326 glide_types dthckdtm glide_geomderv
R 7800 5 327 glide_types dthckdtm$sd glide_geomderv
R 7801 5 328 glide_types dthckdtm$p glide_geomderv
R 7802 5 329 glide_types dthckdtm$o glide_geomderv
R 7806 5 333 glide_types dusrfdtm glide_geomderv
R 7807 5 334 glide_types dusrfdtm$sd glide_geomderv
R 7808 5 335 glide_types dusrfdtm$p glide_geomderv
R 7809 5 336 glide_types dusrfdtm$o glide_geomderv
R 7813 5 340 glide_types stagthck glide_geomderv
R 7814 5 341 glide_types stagthck$sd glide_geomderv
R 7815 5 342 glide_types stagthck$p glide_geomderv
R 7816 5 343 glide_types stagthck$o glide_geomderv
R 7820 5 347 glide_types stagusrf glide_geomderv
R 7821 5 348 glide_types stagusrf$sd glide_geomderv
R 7822 5 349 glide_types stagusrf$p glide_geomderv
R 7823 5 350 glide_types stagusrf$o glide_geomderv
R 7827 5 354 glide_types staglsrf glide_geomderv
R 7828 5 355 glide_types staglsrf$sd glide_geomderv
R 7829 5 356 glide_types staglsrf$p glide_geomderv
R 7830 5 357 glide_types staglsrf$o glide_geomderv
R 7834 5 361 glide_types stagtopg glide_geomderv
R 7835 5 362 glide_types stagtopg$sd glide_geomderv
R 7836 5 363 glide_types stagtopg$p glide_geomderv
R 7837 5 364 glide_types stagtopg$o glide_geomderv
R 7839 25 366 glide_types glide_tensor
R 7843 5 370 glide_types scalar glide_tensor
R 7844 5 371 glide_types scalar$sd glide_tensor
R 7845 5 372 glide_types scalar$p glide_tensor
R 7846 5 373 glide_types scalar$o glide_tensor
R 7851 5 378 glide_types xz glide_tensor
R 7852 5 379 glide_types xz$sd glide_tensor
R 7853 5 380 glide_types xz$p glide_tensor
R 7854 5 381 glide_types xz$o glide_tensor
R 7859 5 386 glide_types yz glide_tensor
R 7860 5 387 glide_types yz$sd glide_tensor
R 7861 5 388 glide_types yz$p glide_tensor
R 7862 5 389 glide_types yz$o glide_tensor
R 7867 5 394 glide_types xx glide_tensor
R 7868 5 395 glide_types xx$sd glide_tensor
R 7869 5 396 glide_types xx$p glide_tensor
R 7870 5 397 glide_types xx$o glide_tensor
R 7875 5 402 glide_types yy glide_tensor
R 7876 5 403 glide_types yy$sd glide_tensor
R 7877 5 404 glide_types yy$p glide_tensor
R 7878 5 405 glide_types yy$o glide_tensor
R 7883 5 410 glide_types xy glide_tensor
R 7884 5 411 glide_types xy$sd glide_tensor
R 7885 5 412 glide_types xy$p glide_tensor
R 7886 5 413 glide_types xy$o glide_tensor
R 7888 25 415 glide_types glide_velocity
R 7892 5 419 glide_types uvel glide_velocity
R 7893 5 420 glide_types uvel$sd glide_velocity
R 7894 5 421 glide_types uvel$p glide_velocity
R 7895 5 422 glide_types uvel$o glide_velocity
R 7900 5 427 glide_types vvel glide_velocity
R 7901 5 428 glide_types vvel$sd glide_velocity
R 7902 5 429 glide_types vvel$p glide_velocity
R 7903 5 430 glide_types vvel$o glide_velocity
R 7908 5 435 glide_types velnorm glide_velocity
R 7909 5 436 glide_types velnorm$sd glide_velocity
R 7910 5 437 glide_types velnorm$p glide_velocity
R 7911 5 438 glide_types velnorm$o glide_velocity
R 7916 5 443 glide_types wvel glide_velocity
R 7917 5 444 glide_types wvel$sd glide_velocity
R 7918 5 445 glide_types wvel$p glide_velocity
R 7919 5 446 glide_types wvel$o glide_velocity
R 7924 5 451 glide_types wgrd glide_velocity
R 7925 5 452 glide_types wgrd$sd glide_velocity
R 7926 5 453 glide_types wgrd$p glide_velocity
R 7927 5 454 glide_types wgrd$o glide_velocity
R 7932 5 459 glide_types surfvel glide_velocity
R 7933 5 460 glide_types surfvel$sd glide_velocity
R 7934 5 461 glide_types surfvel$p glide_velocity
R 7935 5 462 glide_types surfvel$o glide_velocity
R 7939 5 466 glide_types uflx glide_velocity
R 7940 5 467 glide_types uflx$sd glide_velocity
R 7941 5 468 glide_types uflx$p glide_velocity
R 7942 5 469 glide_types uflx$o glide_velocity
R 7946 5 473 glide_types vflx glide_velocity
R 7947 5 474 glide_types vflx$sd glide_velocity
R 7948 5 475 glide_types vflx$p glide_velocity
R 7949 5 476 glide_types vflx$o glide_velocity
R 7953 5 480 glide_types diffu glide_velocity
R 7954 5 481 glide_types diffu$sd glide_velocity
R 7955 5 482 glide_types diffu$p glide_velocity
R 7956 5 483 glide_types diffu$o glide_velocity
R 7960 5 487 glide_types diffu_x glide_velocity
R 7961 5 488 glide_types diffu_x$sd glide_velocity
R 7962 5 489 glide_types diffu_x$p glide_velocity
R 7963 5 490 glide_types diffu_x$o glide_velocity
R 7967 5 494 glide_types diffu_y glide_velocity
R 7968 5 495 glide_types diffu_y$sd glide_velocity
R 7969 5 496 glide_types diffu_y$p glide_velocity
R 7970 5 497 glide_types diffu_y$o glide_velocity
R 7974 5 501 glide_types total_diffu glide_velocity
R 7975 5 502 glide_types total_diffu$sd glide_velocity
R 7976 5 503 glide_types total_diffu$p glide_velocity
R 7977 5 504 glide_types total_diffu$o glide_velocity
R 7981 5 508 glide_types ubas glide_velocity
R 7982 5 509 glide_types ubas$sd glide_velocity
R 7983 5 510 glide_types ubas$p glide_velocity
R 7984 5 511 glide_types ubas$o glide_velocity
R 7988 5 515 glide_types ubas_tavg glide_velocity
R 7989 5 516 glide_types ubas_tavg$sd glide_velocity
R 7990 5 517 glide_types ubas_tavg$p glide_velocity
R 7991 5 518 glide_types ubas_tavg$o glide_velocity
R 7995 5 522 glide_types vbas glide_velocity
R 7996 5 523 glide_types vbas$sd glide_velocity
R 7997 5 524 glide_types vbas$p glide_velocity
R 7998 5 525 glide_types vbas$o glide_velocity
R 8002 5 529 glide_types vbas_tavg glide_velocity
R 8003 5 530 glide_types vbas_tavg$sd glide_velocity
R 8004 5 531 glide_types vbas_tavg$p glide_velocity
R 8005 5 532 glide_types vbas_tavg$o glide_velocity
R 8007 5 534 glide_types is_velocity_valid glide_velocity
R 8010 5 537 glide_types bed_softness glide_velocity
R 8011 5 538 glide_types bed_softness$sd glide_velocity
R 8012 5 539 glide_types bed_softness$p glide_velocity
R 8013 5 540 glide_types bed_softness$o glide_velocity
R 8017 5 544 glide_types btrc glide_velocity
R 8018 5 545 glide_types btrc$sd glide_velocity
R 8019 5 546 glide_types btrc$p glide_velocity
R 8020 5 547 glide_types btrc$o glide_velocity
R 8025 5 552 glide_types btraction glide_velocity
R 8026 5 553 glide_types btraction$sd glide_velocity
R 8027 5 554 glide_types btraction$p glide_velocity
R 8028 5 555 glide_types btraction$o glide_velocity
R 8032 5 559 glide_types beta glide_velocity
R 8033 5 560 glide_types beta$sd glide_velocity
R 8034 5 561 glide_types beta$p glide_velocity
R 8035 5 562 glide_types beta$o glide_velocity
R 8039 5 566 glide_types velmask glide_velocity
R 8040 5 567 glide_types velmask$sd glide_velocity
R 8041 5 568 glide_types velmask$p glide_velocity
R 8042 5 569 glide_types velmask$o glide_velocity
R 8046 5 573 glide_types kinbcmask glide_velocity
R 8047 5 574 glide_types kinbcmask$sd glide_velocity
R 8048 5 575 glide_types kinbcmask$p glide_velocity
R 8049 5 576 glide_types kinbcmask$o glide_velocity
R 8053 5 580 glide_types dynbcmask glide_velocity
R 8054 5 581 glide_types dynbcmask$sd glide_velocity
R 8055 5 582 glide_types dynbcmask$p glide_velocity
R 8056 5 583 glide_types dynbcmask$o glide_velocity
R 8058 25 585 glide_types glide_stress_t
R 8059 5 586 glide_types tau glide_stress_t
R 8063 5 590 glide_types efvs glide_stress_t
R 8064 5 591 glide_types efvs$sd glide_stress_t
R 8065 5 592 glide_types efvs$p glide_stress_t
R 8066 5 593 glide_types efvs$o glide_stress_t
R 8070 5 597 glide_types tau_x glide_stress_t
R 8071 5 598 glide_types tau_x$sd glide_stress_t
R 8072 5 599 glide_types tau_x$p glide_stress_t
R 8073 5 600 glide_types tau_x$o glide_stress_t
R 8077 5 604 glide_types tau_y glide_stress_t
R 8078 5 605 glide_types tau_y$sd glide_stress_t
R 8079 5 606 glide_types tau_y$p glide_stress_t
R 8080 5 607 glide_types tau_y$o glide_stress_t
R 8082 25 609 glide_types glide_climate
R 8085 5 612 glide_types acab glide_climate
R 8086 5 613 glide_types acab$sd glide_climate
R 8087 5 614 glide_types acab$p glide_climate
R 8088 5 615 glide_types acab$o glide_climate
R 8092 5 619 glide_types acab_tavg glide_climate
R 8093 5 620 glide_types acab_tavg$sd glide_climate
R 8094 5 621 glide_types acab_tavg$p glide_climate
R 8095 5 622 glide_types acab_tavg$o glide_climate
R 8099 5 626 glide_types artm glide_climate
R 8100 5 627 glide_types artm$sd glide_climate
R 8101 5 628 glide_types artm$p glide_climate
R 8102 5 629 glide_types artm$o glide_climate
R 8106 5 633 glide_types lati glide_climate
R 8107 5 634 glide_types lati$sd glide_climate
R 8108 5 635 glide_types lati$p glide_climate
R 8109 5 636 glide_types lati$o glide_climate
R 8113 5 640 glide_types loni glide_climate
R 8114 5 641 glide_types loni$sd glide_climate
R 8115 5 642 glide_types loni$p glide_climate
R 8116 5 643 glide_types loni$o glide_climate
R 8120 5 647 glide_types calving glide_climate
R 8121 5 648 glide_types calving$sd glide_climate
R 8122 5 649 glide_types calving$p glide_climate
R 8123 5 650 glide_types calving$o glide_climate
R 8125 5 652 glide_types eus glide_climate
R 8128 5 655 glide_types backstress glide_climate
R 8129 5 656 glide_types backstress$sd glide_climate
R 8130 5 657 glide_types backstress$p glide_climate
R 8131 5 658 glide_types backstress$o glide_climate
R 8135 5 662 glide_types backstressmap glide_climate
R 8136 5 663 glide_types backstressmap$sd glide_climate
R 8137 5 664 glide_types backstressmap$p glide_climate
R 8138 5 665 glide_types backstressmap$o glide_climate
R 8140 5 667 glide_types stressin glide_climate
R 8141 5 668 glide_types stressout glide_climate
R 8142 5 669 glide_types slidconst glide_climate
R 8143 5 670 glide_types tempanmly glide_climate
R 8144 25 671 glide_types glide_temper
R 8148 5 675 glide_types temp glide_temper
R 8149 5 676 glide_types temp$sd glide_temper
R 8150 5 677 glide_types temp$p glide_temper
R 8151 5 678 glide_types temp$o glide_temper
R 8155 5 682 glide_types bheatflx glide_temper
R 8156 5 683 glide_types bheatflx$sd glide_temper
R 8157 5 684 glide_types bheatflx$p glide_temper
R 8158 5 685 glide_types bheatflx$o glide_temper
R 8163 5 690 glide_types flwa glide_temper
R 8164 5 691 glide_types flwa$sd glide_temper
R 8165 5 692 glide_types flwa$p glide_temper
R 8166 5 693 glide_types flwa$o glide_temper
R 8170 5 697 glide_types bwat glide_temper
R 8171 5 698 glide_types bwat$sd glide_temper
R 8172 5 699 glide_types bwat$p glide_temper
R 8173 5 700 glide_types bwat$o glide_temper
R 8177 5 704 glide_types bwatflx glide_temper
R 8178 5 705 glide_types bwatflx$sd glide_temper
R 8179 5 706 glide_types bwatflx$p glide_temper
R 8180 5 707 glide_types bwatflx$o glide_temper
R 8184 5 711 glide_types stagbwat glide_temper
R 8185 5 712 glide_types stagbwat$sd glide_temper
R 8186 5 713 glide_types stagbwat$p glide_temper
R 8187 5 714 glide_types stagbwat$o glide_temper
R 8191 5 718 glide_types bmlt glide_temper
R 8192 5 719 glide_types bmlt$sd glide_temper
R 8193 5 720 glide_types bmlt$p glide_temper
R 8194 5 721 glide_types bmlt$o glide_temper
R 8198 5 725 glide_types bmlt_tavg glide_temper
R 8199 5 726 glide_types bmlt_tavg$sd glide_temper
R 8200 5 727 glide_types bmlt_tavg$p glide_temper
R 8201 5 728 glide_types bmlt_tavg$o glide_temper
R 8205 5 732 glide_types stagbtemp glide_temper
R 8206 5 733 glide_types stagbtemp$sd glide_temper
R 8207 5 734 glide_types stagbtemp$p glide_temper
R 8208 5 735 glide_types stagbtemp$o glide_temper
R 8212 5 739 glide_types bpmp glide_temper
R 8213 5 740 glide_types bpmp$sd glide_temper
R 8214 5 741 glide_types bpmp$p glide_temper
R 8215 5 742 glide_types bpmp$o glide_temper
R 8219 5 746 glide_types stagbpmp glide_temper
R 8220 5 747 glide_types stagbpmp$sd glide_temper
R 8221 5 748 glide_types stagbpmp$p glide_temper
R 8222 5 749 glide_types stagbpmp$o glide_temper
R 8226 5 753 glide_types bfricflx glide_temper
R 8227 5 754 glide_types bfricflx$sd glide_temper
R 8228 5 755 glide_types bfricflx$p glide_temper
R 8229 5 756 glide_types bfricflx$o glide_temper
R 8233 5 760 glide_types ucondflx glide_temper
R 8234 5 761 glide_types ucondflx$sd glide_temper
R 8235 5 762 glide_types ucondflx$p glide_temper
R 8236 5 763 glide_types ucondflx$o glide_temper
R 8240 5 767 glide_types lcondflx glide_temper
R 8241 5 768 glide_types lcondflx$sd glide_temper
R 8242 5 769 glide_types lcondflx$p glide_temper
R 8243 5 770 glide_types lcondflx$o glide_temper
R 8247 5 774 glide_types dissipcol glide_temper
R 8248 5 775 glide_types dissipcol$sd glide_temper
R 8249 5 776 glide_types dissipcol$p glide_temper
R 8250 5 777 glide_types dissipcol$o glide_temper
R 8252 5 779 glide_types niter glide_temper
R 8253 5 780 glide_types perturb glide_temper
R 8254 5 781 glide_types grid glide_temper
R 8255 5 782 glide_types tpt glide_temper
R 8256 5 783 glide_types first1 glide_temper
R 8257 5 784 glide_types newtemps glide_temper
R 8258 25 785 glide_types glide_lithot_type
R 8262 5 789 glide_types temp glide_lithot_type
R 8263 5 790 glide_types temp$sd glide_lithot_type
R 8264 5 791 glide_types temp$p glide_lithot_type
R 8265 5 792 glide_types temp$o glide_lithot_type
R 8269 5 796 glide_types mask glide_lithot_type
R 8270 5 797 glide_types mask$sd glide_lithot_type
R 8271 5 798 glide_types mask$p glide_lithot_type
R 8272 5 799 glide_types mask$o glide_lithot_type
R 8274 5 801 glide_types num_dim glide_lithot_type
R 8275 5 802 glide_types fd_coeff glide_lithot_type
R 8276 5 803 glide_types fd_coeff_slap glide_lithot_type
R 8277 5 804 glide_types all_bar_top glide_lithot_type
R 8279 5 806 glide_types rhs glide_lithot_type
R 8280 5 807 glide_types rhs$sd glide_lithot_type
R 8281 5 808 glide_types rhs$p glide_lithot_type
R 8282 5 809 glide_types rhs$o glide_lithot_type
R 8285 5 812 glide_types answer glide_lithot_type
R 8286 5 813 glide_types answer$sd glide_lithot_type
R 8287 5 814 glide_types answer$p glide_lithot_type
R 8288 5 815 glide_types answer$o glide_lithot_type
R 8291 5 818 glide_types supd glide_lithot_type
R 8292 5 819 glide_types supd$sd glide_lithot_type
R 8293 5 820 glide_types supd$p glide_lithot_type
R 8294 5 821 glide_types supd$o glide_lithot_type
R 8296 5 823 glide_types diag glide_lithot_type
R 8298 5 825 glide_types diag$sd glide_lithot_type
R 8299 5 826 glide_types diag$p glide_lithot_type
R 8300 5 827 glide_types diag$o glide_lithot_type
R 8302 5 829 glide_types subd glide_lithot_type
R 8304 5 831 glide_types subd$sd glide_lithot_type
R 8305 5 832 glide_types subd$p glide_lithot_type
R 8306 5 833 glide_types subd$o glide_lithot_type
R 8309 5 836 glide_types rwork glide_lithot_type
R 8310 5 837 glide_types rwork$sd glide_lithot_type
R 8311 5 838 glide_types rwork$p glide_lithot_type
R 8312 5 839 glide_types rwork$o glide_lithot_type
R 8315 5 842 glide_types iwork glide_lithot_type
R 8316 5 843 glide_types iwork$sd glide_lithot_type
R 8317 5 844 glide_types iwork$p glide_lithot_type
R 8318 5 845 glide_types iwork$o glide_lithot_type
R 8320 5 847 glide_types mxnelt glide_lithot_type
R 8322 5 849 glide_types deltaz glide_lithot_type
R 8323 5 850 glide_types deltaz$sd glide_lithot_type
R 8324 5 851 glide_types deltaz$p glide_lithot_type
R 8325 5 852 glide_types deltaz$o glide_lithot_type
R 8329 5 856 glide_types zfactors glide_lithot_type
R 8330 5 857 glide_types zfactors$sd glide_lithot_type
R 8331 5 858 glide_types zfactors$p glide_lithot_type
R 8332 5 859 glide_types zfactors$o glide_lithot_type
R 8334 5 861 glide_types xfactor glide_lithot_type
R 8335 5 862 glide_types yfactor glide_lithot_type
R 8336 5 863 glide_types surft glide_lithot_type
R 8337 5 864 glide_types mart glide_lithot_type
R 8338 5 865 glide_types nlayer glide_lithot_type
R 8339 5 866 glide_types rock_base glide_lithot_type
R 8340 5 867 glide_types numt glide_lithot_type
R 8341 5 868 glide_types rho_r glide_lithot_type
R 8342 5 869 glide_types shc_r glide_lithot_type
R 8343 5 870 glide_types con_r glide_lithot_type
R 8344 5 871 glide_types diffu glide_lithot_type
R 8345 25 872 glide_types glide_funits
R 8346 5 873 glide_types sigfile glide_funits
R 8347 5 874 glide_types ncfile glide_funits
R 8348 5 875 glide_types out_first glide_funits
R 8350 5 877 glide_types out_first$p glide_funits
R 8352 5 879 glide_types in_first glide_funits
R 8354 5 881 glide_types in_first$p glide_funits
R 8356 25 883 glide_types glide_numerics
R 8357 5 884 glide_types tstart glide_numerics
R 8358 5 885 glide_types tend glide_numerics
R 8359 5 886 glide_types time glide_numerics
R 8360 5 887 glide_types tinc glide_numerics
R 8361 5 888 glide_types ntem glide_numerics
R 8362 5 889 glide_types nvel glide_numerics
R 8363 5 890 glide_types alpha glide_numerics
R 8364 5 891 glide_types alphas glide_numerics
R 8365 5 892 glide_types thklim glide_numerics
R 8366 5 893 glide_types mlimit glide_numerics
R 8367 5 894 glide_types calving_fraction glide_numerics
R 8368 5 895 glide_types dew glide_numerics
R 8369 5 896 glide_types dns glide_numerics
R 8370 5 897 glide_types dt glide_numerics
R 8371 5 898 glide_types dttem glide_numerics
R 8372 5 899 glide_types nshlf glide_numerics
R 8373 5 900 glide_types subcyc glide_numerics
R 8374 5 901 glide_types timecounter glide_numerics
R 8376 5 903 glide_types sigma glide_numerics
R 8377 5 904 glide_types sigma$sd glide_numerics
R 8378 5 905 glide_types sigma$p glide_numerics
R 8379 5 906 glide_types sigma$o glide_numerics
R 8382 5 909 glide_types stagsigma glide_numerics
R 8383 5 910 glide_types stagsigma$sd glide_numerics
R 8384 5 911 glide_types stagsigma$p glide_numerics
R 8385 5 912 glide_types stagsigma$o glide_numerics
R 8388 5 915 glide_types stagwbndsigma glide_numerics
R 8389 5 916 glide_types stagwbndsigma$sd glide_numerics
R 8390 5 917 glide_types stagwbndsigma$p glide_numerics
R 8391 5 918 glide_types stagwbndsigma$o glide_numerics
R 8393 5 920 glide_types profile_period glide_numerics
R 8394 5 921 glide_types ndiag glide_numerics
R 8395 5 922 glide_types idiag glide_numerics
R 8396 5 923 glide_types jdiag glide_numerics
R 8397 25 924 glide_types glide_grnd
R 8400 5 927 glide_types gl_ew glide_grnd
R 8401 5 928 glide_types gl_ew$sd glide_grnd
R 8402 5 929 glide_types gl_ew$p glide_grnd
R 8403 5 930 glide_types gl_ew$o glide_grnd
R 8407 5 934 glide_types gl_ns glide_grnd
R 8408 5 935 glide_types gl_ns$sd glide_grnd
R 8409 5 936 glide_types gl_ns$p glide_grnd
R 8410 5 937 glide_types gl_ns$o glide_grnd
R 8414 5 941 glide_types gline_flux glide_grnd
R 8415 5 942 glide_types gline_flux$sd glide_grnd
R 8416 5 943 glide_types gline_flux$p glide_grnd
R 8417 5 944 glide_types gline_flux$o glide_grnd
R 8419 25 946 glide_types glide_velowk
R 8421 5 948 glide_types depth glide_velowk
R 8422 5 949 glide_types depth$sd glide_velowk
R 8423 5 950 glide_types depth$p glide_velowk
R 8424 5 951 glide_types depth$o glide_velowk
R 8427 5 954 glide_types dupsw glide_velowk
R 8428 5 955 glide_types dupsw$sd glide_velowk
R 8429 5 956 glide_types dupsw$p glide_velowk
R 8430 5 957 glide_types dupsw$o glide_velowk
R 8433 5 960 glide_types depthw glide_velowk
R 8434 5 961 glide_types depthw$sd glide_velowk
R 8435 5 962 glide_types depthw$p glide_velowk
R 8436 5 963 glide_types depthw$o glide_velowk
R 8439 5 966 glide_types suvel glide_velowk
R 8440 5 967 glide_types suvel$sd glide_velowk
R 8441 5 968 glide_types suvel$p glide_velowk
R 8442 5 969 glide_types suvel$o glide_velowk
R 8445 5 972 glide_types svvel glide_velowk
R 8446 5 973 glide_types svvel$sd glide_velowk
R 8447 5 974 glide_types svvel$p glide_velowk
R 8448 5 975 glide_types svvel$o glide_velowk
R 8452 5 979 glide_types fslip glide_velowk
R 8453 5 980 glide_types fslip$sd glide_velowk
R 8454 5 981 glide_types fslip$p glide_velowk
R 8455 5 982 glide_types fslip$o glide_velowk
R 8459 5 986 glide_types dintflwa glide_velowk
R 8460 5 987 glide_types dintflwa$sd glide_velowk
R 8461 5 988 glide_types dintflwa$p glide_velowk
R 8462 5 989 glide_types dintflwa$o glide_velowk
R 8465 5 992 glide_types dups glide_velowk
R 8466 5 993 glide_types dups$sd glide_velowk
R 8467 5 994 glide_types dups$p glide_velowk
R 8468 5 995 glide_types dups$o glide_velowk
R 8470 5 997 glide_types fact glide_velowk
R 8471 5 998 glide_types c glide_velowk
R 8472 5 999 glide_types watwd glide_velowk
R 8473 5 1000 glide_types watct glide_velowk
R 8474 5 1001 glide_types trc0 glide_velowk
R 8475 5 1002 glide_types trcmin glide_velowk
R 8476 5 1003 glide_types marine glide_velowk
R 8477 5 1004 glide_types trcmax glide_velowk
R 8478 5 1005 glide_types btrac_const glide_velowk
R 8479 5 1006 glide_types btrac_slope glide_velowk
R 8480 5 1007 glide_types btrac_max glide_velowk
R 8481 25 1008 glide_types glide_pcgdwk
R 8482 5 1009 glide_types matrix glide_pcgdwk
R 8484 5 1011 glide_types rhsd glide_pcgdwk
R 8485 5 1012 glide_types rhsd$sd glide_pcgdwk
R 8486 5 1013 glide_types rhsd$p glide_pcgdwk
R 8487 5 1014 glide_types rhsd$o glide_pcgdwk
R 8490 5 1017 glide_types answ glide_pcgdwk
R 8491 5 1018 glide_types answ$sd glide_pcgdwk
R 8492 5 1019 glide_types answ$p glide_pcgdwk
R 8493 5 1020 glide_types answ$o glide_pcgdwk
R 8495 5 1022 glide_types fc glide_pcgdwk
R 8496 5 1023 glide_types fc2 glide_pcgdwk
R 8497 5 1024 glide_types ct glide_pcgdwk
R 8498 25 1025 glide_types glide_thckwk
R 8501 5 1028 glide_types oldthck glide_thckwk
R 8502 5 1029 glide_types oldthck$sd glide_thckwk
R 8503 5 1030 glide_types oldthck$p glide_thckwk
R 8504 5 1031 glide_types oldthck$o glide_thckwk
R 8508 5 1035 glide_types oldthck2 glide_thckwk
R 8509 5 1036 glide_types oldthck2$sd glide_thckwk
R 8510 5 1037 glide_types oldthck2$p glide_thckwk
R 8511 5 1038 glide_types oldthck2$o glide_thckwk
R 8515 5 1042 glide_types float glide_thckwk
R 8516 5 1043 glide_types float$sd glide_thckwk
R 8517 5 1044 glide_types float$p glide_thckwk
R 8518 5 1045 glide_types float$o glide_thckwk
R 8523 5 1050 glide_types olds glide_thckwk
R 8524 5 1051 glide_types olds$sd glide_thckwk
R 8525 5 1052 glide_types olds$p glide_thckwk
R 8526 5 1053 glide_types olds$o glide_thckwk
R 8528 5 1055 glide_types nwhich glide_thckwk
R 8529 5 1056 glide_types oldtime glide_thckwk
R 8531 5 1058 glide_types alpha glide_thckwk
R 8532 5 1059 glide_types alpha$sd glide_thckwk
R 8533 5 1060 glide_types alpha$p glide_thckwk
R 8534 5 1061 glide_types alpha$o glide_thckwk
R 8537 5 1064 glide_types beta glide_thckwk
R 8538 5 1065 glide_types beta$sd glide_thckwk
R 8539 5 1066 glide_types beta$p glide_thckwk
R 8540 5 1067 glide_types beta$o glide_thckwk
R 8543 5 1070 glide_types gamma glide_thckwk
R 8544 5 1071 glide_types gamma$sd glide_thckwk
R 8545 5 1072 glide_types gamma$p glide_thckwk
R 8546 5 1073 glide_types gamma$o glide_thckwk
R 8549 5 1076 glide_types delta glide_thckwk
R 8550 5 1077 glide_types delta$sd glide_thckwk
R 8551 5 1078 glide_types delta$p glide_thckwk
R 8552 5 1079 glide_types delta$o glide_thckwk
R 8554 25 1081 glide_types glide_tempwk
R 8558 5 1085 glide_types inittemp glide_tempwk
R 8559 5 1086 glide_types inittemp$sd glide_tempwk
R 8560 5 1087 glide_types inittemp$p glide_tempwk
R 8561 5 1088 glide_types inittemp$o glide_tempwk
R 8566 5 1093 glide_types dissip glide_tempwk
R 8567 5 1094 glide_types dissip$sd glide_tempwk
R 8568 5 1095 glide_types dissip$p glide_tempwk
R 8569 5 1096 glide_types dissip$o glide_tempwk
R 8574 5 1101 glide_types compheat glide_tempwk
R 8575 5 1102 glide_types compheat$sd glide_tempwk
R 8576 5 1103 glide_types compheat$p glide_tempwk
R 8577 5 1104 glide_types compheat$o glide_tempwk
R 8582 5 1109 glide_types initadvt glide_tempwk
R 8583 5 1110 glide_types initadvt$sd glide_tempwk
R 8584 5 1111 glide_types initadvt$p glide_tempwk
R 8585 5 1112 glide_types initadvt$o glide_tempwk
R 8588 5 1115 glide_types dupa glide_tempwk
R 8589 5 1116 glide_types dupa$sd glide_tempwk
R 8590 5 1117 glide_types dupa$p glide_tempwk
R 8591 5 1118 glide_types dupa$o glide_tempwk
R 8594 5 1121 glide_types dupb glide_tempwk
R 8595 5 1122 glide_types dupb$sd glide_tempwk
R 8596 5 1123 glide_types dupb$p glide_tempwk
R 8597 5 1124 glide_types dupb$o glide_tempwk
R 8600 5 1127 glide_types dupc glide_tempwk
R 8601 5 1128 glide_types dupc$sd glide_tempwk
R 8602 5 1129 glide_types dupc$p glide_tempwk
R 8603 5 1130 glide_types dupc$o glide_tempwk
R 8606 5 1133 glide_types c1 glide_tempwk
R 8607 5 1134 glide_types c1$sd glide_tempwk
R 8608 5 1135 glide_types c1$p glide_tempwk
R 8609 5 1136 glide_types c1$o glide_tempwk
R 8613 5 1140 glide_types dups glide_tempwk
R 8614 5 1141 glide_types dups$sd glide_tempwk
R 8615 5 1142 glide_types dups$p glide_tempwk
R 8616 5 1143 glide_types dups$o glide_tempwk
R 8620 5 1147 glide_types wphi glide_tempwk
R 8621 5 1148 glide_types wphi$sd glide_tempwk
R 8622 5 1149 glide_types wphi$p glide_tempwk
R 8623 5 1150 glide_types wphi$o glide_tempwk
R 8627 5 1154 glide_types bwatu glide_tempwk
R 8628 5 1155 glide_types bwatu$sd glide_tempwk
R 8629 5 1156 glide_types bwatu$p glide_tempwk
R 8630 5 1157 glide_types bwatu$o glide_tempwk
R 8634 5 1161 glide_types bwatv glide_tempwk
R 8635 5 1162 glide_types bwatv$sd glide_tempwk
R 8636 5 1163 glide_types bwatv$p glide_tempwk
R 8637 5 1164 glide_types bwatv$o glide_tempwk
R 8641 5 1168 glide_types fluxew glide_tempwk
R 8642 5 1169 glide_types fluxew$sd glide_tempwk
R 8643 5 1170 glide_types fluxew$p glide_tempwk
R 8644 5 1171 glide_types fluxew$o glide_tempwk
R 8648 5 1175 glide_types fluxns glide_tempwk
R 8649 5 1176 glide_types fluxns$sd glide_tempwk
R 8650 5 1177 glide_types fluxns$p glide_tempwk
R 8651 5 1178 glide_types fluxns$o glide_tempwk
R 8655 5 1182 glide_types bint glide_tempwk
R 8656 5 1183 glide_types bint$sd glide_tempwk
R 8657 5 1184 glide_types bint$p glide_tempwk
R 8658 5 1185 glide_types bint$o glide_tempwk
R 8662 5 1189 glide_types smth glide_tempwk
R 8663 5 1190 glide_types smth$sd glide_tempwk
R 8664 5 1191 glide_types smth$p glide_tempwk
R 8665 5 1192 glide_types smth$o glide_tempwk
R 8670 5 1197 glide_types hadv_u glide_tempwk
R 8671 5 1198 glide_types hadv_u$sd glide_tempwk
R 8672 5 1199 glide_types hadv_u$p glide_tempwk
R 8673 5 1200 glide_types hadv_u$o glide_tempwk
R 8678 5 1205 glide_types hadv_v glide_tempwk
R 8679 5 1206 glide_types hadv_v$sd glide_tempwk
R 8680 5 1207 glide_types hadv_v$p glide_tempwk
R 8681 5 1208 glide_types hadv_v$o glide_tempwk
R 8683 5 1210 glide_types cons glide_tempwk
R 8684 5 1211 glide_types f glide_tempwk
R 8685 5 1212 glide_types c glide_tempwk
R 8686 5 1213 glide_types slide_f glide_tempwk
R 8687 5 1214 glide_types noflow glide_tempwk
R 8688 5 1215 glide_types advconst glide_tempwk
R 8689 5 1216 glide_types zbed glide_tempwk
R 8690 5 1217 glide_types dupn glide_tempwk
R 8691 5 1218 glide_types wmax glide_tempwk
R 8692 5 1219 glide_types dt_wat glide_tempwk
R 8693 5 1220 glide_types watvel glide_tempwk
R 8694 5 1221 glide_types nwat glide_tempwk
R 8695 25 1222 glide_types glide_gridwk
R 8698 5 1225 glide_types hte glide_gridwk
R 8699 5 1226 glide_types hte$sd glide_gridwk
R 8700 5 1227 glide_types hte$p glide_gridwk
R 8701 5 1228 glide_types hte$o glide_gridwk
R 8705 5 1232 glide_types htn glide_gridwk
R 8706 5 1233 glide_types htn$sd glide_gridwk
R 8707 5 1234 glide_types htn$p glide_gridwk
R 8708 5 1235 glide_types htn$o glide_gridwk
R 8712 5 1239 glide_types dxt glide_gridwk
R 8713 5 1240 glide_types dxt$sd glide_gridwk
R 8714 5 1241 glide_types dxt$p glide_gridwk
R 8715 5 1242 glide_types dxt$o glide_gridwk
R 8719 5 1246 glide_types dyt glide_gridwk
R 8720 5 1247 glide_types dyt$sd glide_gridwk
R 8721 5 1248 glide_types dyt$p glide_gridwk
R 8722 5 1249 glide_types dyt$o glide_gridwk
R 8726 5 1253 glide_types tarea glide_gridwk
R 8727 5 1254 glide_types tarea$sd glide_gridwk
R 8728 5 1255 glide_types tarea$p glide_gridwk
R 8729 5 1256 glide_types tarea$o glide_gridwk
R 8733 5 1260 glide_types tarear glide_gridwk
R 8734 5 1261 glide_types tarear$sd glide_gridwk
R 8735 5 1262 glide_types tarear$p glide_gridwk
R 8736 5 1263 glide_types tarear$o glide_gridwk
R 8740 5 1267 glide_types mask glide_gridwk
R 8741 5 1268 glide_types mask$sd glide_gridwk
R 8742 5 1269 glide_types mask$p glide_gridwk
R 8743 5 1270 glide_types mask$o glide_gridwk
R 8747 5 1274 glide_types xav glide_gridwk
R 8748 5 1275 glide_types xav$sd glide_gridwk
R 8749 5 1276 glide_types xav$p glide_gridwk
R 8750 5 1277 glide_types xav$o glide_gridwk
R 8754 5 1281 glide_types yav glide_gridwk
R 8755 5 1282 glide_types yav$sd glide_gridwk
R 8756 5 1283 glide_types yav$p glide_gridwk
R 8757 5 1284 glide_types yav$o glide_gridwk
R 8761 5 1288 glide_types xxav glide_gridwk
R 8762 5 1289 glide_types xxav$sd glide_gridwk
R 8763 5 1290 glide_types xxav$p glide_gridwk
R 8764 5 1291 glide_types xxav$o glide_gridwk
R 8768 5 1295 glide_types xyav glide_gridwk
R 8769 5 1296 glide_types xyav$sd glide_gridwk
R 8770 5 1297 glide_types xyav$p glide_gridwk
R 8771 5 1298 glide_types xyav$o glide_gridwk
R 8775 5 1302 glide_types yyav glide_gridwk
R 8776 5 1303 glide_types yyav$sd glide_gridwk
R 8777 5 1304 glide_types yyav$p glide_gridwk
R 8778 5 1305 glide_types yyav$o glide_gridwk
R 8780 25 1307 glide_types glide_paramets
R 8781 5 1308 glide_types bpar glide_paramets
R 8782 5 1309 glide_types btrac_const glide_paramets
R 8783 5 1310 glide_types btrac_slope glide_paramets
R 8784 5 1311 glide_types btrac_max glide_paramets
R 8785 5 1312 glide_types geot glide_paramets
R 8786 5 1313 glide_types flow_factor glide_paramets
R 8787 5 1314 glide_types slip_ratio glide_paramets
R 8788 5 1315 glide_types hydtim glide_paramets
R 8789 5 1316 glide_types bwat_smooth glide_paramets
R 8790 5 1317 glide_types default_flwa glide_paramets
R 8791 25 1318 glide_types glide_basalproc
R 8792 5 1319 glide_types fric glide_basalproc
R 8793 5 1320 glide_types etillo glide_basalproc
R 8794 5 1321 glide_types no glide_basalproc
R 8795 5 1322 glide_types comp glide_basalproc
R 8796 5 1323 glide_types cv glide_basalproc
R 8797 5 1324 glide_types kh glide_basalproc
R 8798 5 1325 glide_types zs glide_basalproc
R 8799 5 1326 glide_types aconst glide_basalproc
R 8800 5 1327 glide_types bconst glide_basalproc
R 8801 5 1328 glide_types till_hot glide_basalproc
R 8802 5 1329 glide_types tnodes glide_basalproc
R 8804 5 1331 glide_types till_dz glide_basalproc
R 8805 5 1332 glide_types till_dz$sd glide_basalproc
R 8806 5 1333 glide_types till_dz$p glide_basalproc
R 8807 5 1334 glide_types till_dz$o glide_basalproc
R 8811 5 1338 glide_types mintauf glide_basalproc
R 8812 5 1339 glide_types mintauf$sd glide_basalproc
R 8813 5 1340 glide_types mintauf$p glide_basalproc
R 8814 5 1341 glide_types mintauf$o glide_basalproc
R 8818 5 1345 glide_types hwater glide_basalproc
R 8819 5 1346 glide_types hwater$sd glide_basalproc
R 8820 5 1347 glide_types hwater$p glide_basalproc
R 8821 5 1348 glide_types hwater$o glide_basalproc
R 8826 5 1353 glide_types u glide_basalproc
R 8827 5 1354 glide_types u$sd glide_basalproc
R 8828 5 1355 glide_types u$p glide_basalproc
R 8829 5 1356 glide_types u$o glide_basalproc
R 8834 5 1361 glide_types etill glide_basalproc
R 8835 5 1362 glide_types etill$sd glide_basalproc
R 8836 5 1363 glide_types etill$p glide_basalproc
R 8837 5 1364 glide_types etill$o glide_basalproc
R 8848 25 1375 glide_types glide_phaml
R 8851 5 1378 glide_types uphaml glide_phaml
R 8852 5 1379 glide_types uphaml$sd glide_phaml
R 8853 5 1380 glide_types uphaml$p glide_phaml
R 8854 5 1381 glide_types uphaml$o glide_phaml
R 8858 5 1385 glide_types init_phaml glide_phaml
R 8859 5 1386 glide_types init_phaml$sd glide_phaml
R 8860 5 1387 glide_types init_phaml$p glide_phaml
R 8861 5 1388 glide_types init_phaml$o glide_phaml
R 8865 5 1392 glide_types rs_phaml glide_phaml
R 8866 5 1393 glide_types rs_phaml$sd glide_phaml
R 8867 5 1394 glide_types rs_phaml$p glide_phaml
R 8868 5 1395 glide_types rs_phaml$o glide_phaml
R 8870 25 1397 glide_types glide_global_type
R 8871 5 1398 glide_types model_id glide_global_type
R 8872 5 1399 glide_types general glide_global_type
R 8873 5 1400 glide_types options glide_global_type
R 8874 5 1401 glide_types geometry glide_global_type
R 8875 5 1402 glide_types geomderv glide_global_type
R 8876 5 1403 glide_types velocity glide_global_type
R 8877 5 1404 glide_types stress glide_global_type
R 8878 5 1405 glide_types climate glide_global_type
R 8879 5 1406 glide_types temper glide_global_type
R 8880 5 1407 glide_types lithot glide_global_type
R 8881 5 1408 glide_types funits glide_global_type
R 8882 5 1409 glide_types numerics glide_global_type
R 8883 5 1410 glide_types velowk glide_global_type
R 8884 5 1411 glide_types pcgdwk glide_global_type
R 8885 5 1412 glide_types thckwk glide_global_type
R 8886 5 1413 glide_types tempwk glide_global_type
R 8887 5 1414 glide_types gridwk glide_global_type
R 8888 5 1415 glide_types paramets glide_global_type
R 8889 5 1416 glide_types projection glide_global_type
R 8890 5 1417 glide_types basalproc glide_global_type
R 8891 5 1418 glide_types profile glide_global_type
R 8892 5 1419 glide_types glide_prof glide_global_type
R 8893 5 1420 glide_types isos glide_global_type
R 8894 5 1421 glide_types phaml glide_global_type
R 8895 5 1422 glide_types ground glide_global_type
R 8896 5 1423 glide_types remap_wk glide_global_type
R 8897 25 1424 glide_types pass_through
R 8898 5 1425 glide_types model pass_through
R 8901 5 1428 glide_types ui pass_through
R 8902 5 1429 glide_types ui$sd pass_through
R 8903 5 1430 glide_types ui$p pass_through
R 8904 5 1431 glide_types ui$o pass_through
R 8908 5 1435 glide_types um pass_through
R 8909 5 1436 glide_types um$sd pass_through
R 8910 5 1437 glide_types um$p pass_through
R 8911 5 1438 glide_types um$o pass_through
R 8915 5 1442 glide_types d2thckcross pass_through
R 8916 5 1443 glide_types d2thckcross$sd pass_through
R 8917 5 1444 glide_types d2thckcross$p pass_through
R 8918 5 1445 glide_types d2thckcross$o pass_through
R 8922 5 1449 glide_types d2usrfcross pass_through
R 8923 5 1450 glide_types d2usrfcross$sd pass_through
R 8924 5 1451 glide_types d2usrfcross$p pass_through
R 8925 5 1452 glide_types d2usrfcross$o pass_through
R 8927 5 1454 glide_types pcgsize pass_through
R 8929 5 1456 glide_types gxf pass_through
R 8930 5 1457 glide_types gxf$sd pass_through
R 8931 5 1458 glide_types gxf$p pass_through
R 8932 5 1459 glide_types gxf$o pass_through
R 8934 5 1461 glide_types l2norm pass_through
R 8935 5 1462 glide_types matrixa pass_through
R 8936 5 1463 glide_types matrixc pass_through
S 9000 23 5 0 0 0 9001 582 43241 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_init_dycore_interface
S 9001 14 5 0 0 0 1 9000 43241 0 400000 A 0 0 0 0 0 0 0 1891 0 0 0 0 0 0 0 0 0 0 0 0 0 15 0 582 0 0 0 0 gtd_init_dycore_interface
F 9001 0
S 9002 23 5 0 0 0 9003 582 43267 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_delete_dycore_interface
S 9003 14 5 0 0 0 1 9002 43267 0 400000 A 0 0 0 0 0 0 0 1892 0 0 0 0 0 0 0 0 0 0 0 0 0 20 0 582 0 0 0 0 gtd_delete_dycore_interface
F 9003 0
S 9004 23 5 0 0 0 9007 582 43295 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_init_dycore
S 9005 1 3 0 0 4233 1 9004 42075 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 model
S 9006 1 3 0 0 6 1 9004 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9007 14 5 0 0 0 1 9004 43295 0 400000 A 0 0 0 0 0 0 0 1893 2 0 0 0 0 0 0 0 0 0 0 0 0 24 0 582 0 0 0 0 gtd_init_dycore
F 9007 2 9005 9006
S 9008 23 5 0 0 0 9012 582 43330 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_run_dycore
S 9009 1 3 0 0 6 1 9008 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9010 1 3 0 0 8 1 9008 43345 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cur_time
S 9011 1 3 0 0 8 1 9008 43354 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 time_inc
S 9012 14 5 0 0 0 1 9008 43330 0 400000 A 0 0 0 0 0 0 0 1896 3 0 0 0 0 0 0 0 0 0 0 0 0 50 0 582 0 0 0 0 gtd_run_dycore
F 9012 3 9009 9010 9011
S 9013 23 5 0 0 0 9015 582 43363 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_delete_dycore
S 9014 1 3 0 0 6 1 9013 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9015 14 5 0 0 0 1 9013 43363 0 400000 A 0 0 0 0 0 0 0 1900 1 0 0 0 0 0 0 0 0 0 0 0 0 57 0 582 0 0 0 0 gtd_delete_dycore
F 9015 1 9014
S 9016 23 5 0 0 0 9019 582 43381 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_set_dim_info
S 9017 7 3 1 0 4350 1 9016 2869 20000004 10003000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shape
S 9018 7 3 3 0 4353 1 9016 43398 20000004 10003000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dim_info
S 9019 14 5 0 0 0 1 9016 43381 20000000 400000 A 0 0 0 0 0 0 0 1902 2 0 0 0 0 0 0 0 0 0 0 0 0 63 0 582 0 0 0 0 gtd_set_dim_info
F 9019 2 9017 9018
S 9020 6 1 0 0 6 1 9016 43407 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_3
S 9021 6 1 0 0 6 1 9016 43415 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_2
S 9022 6 1 0 0 6 1 9016 43423 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_2
S 9023 6 1 0 0 6 1 9016 43431 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_6474
S 9024 6 1 0 0 6 1 9016 43440 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_2
S 9025 6 1 0 0 6 1 9016 43448 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_2
S 9026 6 1 0 0 6 1 9016 43456 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7_2
S 9027 6 1 0 0 6 1 9016 43464 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_6481
S 9028 23 5 0 0 0 9031 582 43473 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_set_geometry_vars
S 9029 1 3 0 0 4233 1 9028 42075 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 model
S 9030 1 3 0 0 6 1 9028 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9031 14 5 0 0 0 1 9028 43473 0 400000 A 0 0 0 0 0 0 0 1905 2 0 0 0 0 0 0 0 0 0 0 0 0 72 0 582 0 0 0 0 gtd_set_geometry_vars
F 9031 2 9029 9030
S 9032 23 5 0 0 0 9035 582 43495 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_set_velocity_vars
S 9033 1 3 0 0 4233 1 9032 42075 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 model
S 9034 1 3 0 0 6 1 9032 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9035 14 5 0 0 0 1 9032 43495 0 400000 A 0 0 0 0 0 0 0 1908 2 0 0 0 0 0 0 0 0 0 0 0 0 121 0 582 0 0 0 0 gtd_set_velocity_vars
F 9035 2 9033 9034
S 9036 23 5 0 0 0 9039 582 43517 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_set_numerics_vars
S 9037 1 3 0 0 4233 1 9036 42075 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 model
S 9038 1 3 0 0 6 1 9036 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9039 14 5 0 0 0 1 9036 43517 0 400000 A 0 0 0 0 0 0 0 1911 2 0 0 0 0 0 0 0 0 0 0 0 0 160 0 582 0 0 0 0 gtd_set_numerics_vars
F 9039 2 9037 9038
S 9040 23 5 0 0 0 9043 582 43539 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_set_temper_vars
S 9041 1 3 0 0 4233 1 9040 42075 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 model
S 9042 1 3 0 0 6 1 9040 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9043 14 5 0 0 0 1 9040 43539 0 400000 A 0 0 0 0 0 0 0 1914 2 0 0 0 0 0 0 0 0 0 0 0 0 181 0 582 0 0 0 0 gtd_set_temper_vars
F 9043 2 9041 9042
S 9044 23 5 0 0 0 9047 582 43559 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gtd_set_climate_vars
S 9045 1 3 0 0 4233 1 9044 42075 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 model
S 9046 1 3 0 0 6 1 9044 43311 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dycore_model_index
S 9047 14 5 0 0 0 1 9044 43559 0 400000 A 0 0 0 0 0 0 0 1917 2 0 0 0 0 0 0 0 0 0 0 0 0 204 0 582 0 0 0 0 gtd_set_climate_vars
F 9047 2 9045 9046
A 12 2 0 0 0 6 592 0 0 0 12 0 0 0 0 0 0 0 0 0
A 14 2 0 0 0 6 585 0 0 0 14 0 0 0 0 0 0 0 0 0
A 16 2 0 0 0 6 586 0 0 0 16 0 0 0 0 0 0 0 0 0
A 32 2 0 0 0 6 613 0 0 0 32 0 0 0 0 0 0 0 0 0
A 36 2 0 0 0 6 610 0 0 0 36 0 0 0 0 0 0 0 0 0
A 39 2 0 0 0 6 611 0 0 0 39 0 0 0 0 0 0 0 0 0
A 41 2 0 0 0 6 612 0 0 0 41 0 0 0 0 0 0 0 0 0
A 43 2 0 0 0 63 614 0 0 0 43 0 0 0 0 0 0 0 0 0
A 44 2 0 0 0 63 615 0 0 0 44 0 0 0 0 0 0 0 0 0
A 45 2 0 0 0 63 616 0 0 0 45 0 0 0 0 0 0 0 0 0
A 46 2 0 0 0 63 617 0 0 0 46 0 0 0 0 0 0 0 0 0
A 47 2 0 0 0 63 618 0 0 0 47 0 0 0 0 0 0 0 0 0
A 48 2 0 0 0 63 619 0 0 0 48 0 0 0 0 0 0 0 0 0
A 57 2 0 0 0 6 620 0 0 0 57 0 0 0 0 0 0 0 0 0
A 66 1 0 1 0 68 632 0 0 0 0 0 0 0 0 0 0 0 0 0
A 68 2 0 0 0 6 672 0 0 0 68 0 0 0 0 0 0 0 0 0
A 69 2 0 0 0 6 668 0 0 0 69 0 0 0 0 0 0 0 0 0
A 183 2 0 0 0 6 916 0 0 0 183 0 0 0 0 0 0 0 0 0
A 281 2 0 0 72 6 1172 0 0 0 281 0 0 0 0 0 0 0 0 0
A 3499 2 0 0 3120 6 6571 0 0 0 3499 0 0 0 0 0 0 0 0 0
A 3503 2 0 0 3125 8 6573 0 0 0 3503 0 0 0 0 0 0 0 0 0
A 3514 2 0 0 3214 16 6574 0 0 0 3514 0 0 0 0 0 0 0 0 0
A 3515 2 0 0 3193 16 6575 0 0 0 3515 0 0 0 0 0 0 0 0 0
A 3516 2 0 0 3454 8 573 0 0 0 3516 0 0 0 0 0 0 0 0 0
A 3517 2 0 0 3195 20 6576 0 0 0 3517 0 0 0 0 0 0 0 0 0
A 3518 2 0 0 2642 2066 6578 0 0 0 3518 0 0 0 0 0 0 0 0 0
A 3519 2 0 0 3222 6 6700 0 0 0 3519 0 0 0 0 0 0 0 0 0
A 3574 2 0 0 3220 8 6699 0 0 0 3574 0 0 0 0 0 0 0 0 0
A 3575 2 0 0 3224 8 6703 0 0 0 3575 0 0 0 0 0 0 0 0 0
A 3576 2 0 0 3226 8 6704 0 0 0 3576 0 0 0 0 0 0 0 0 0
A 3631 2 0 0 1424 6 6813 0 0 0 3631 0 0 0 0 0 0 0 0 0
A 3732 2 0 0 3436 9 6994 0 0 0 3732 0 0 0 0 0 0 0 0 0
A 3736 2 0 0 3064 9 6996 0 0 0 3736 0 0 0 0 0 0 0 0 0
A 6439 2 0 0 5287 2613 6578 0 0 0 6439 0 0 0 0 0 0 0 0 0
A 6440 2 0 0 3971 8 7450 0 0 0 6440 0 0 0 0 0 0 0 0 0
A 6441 2 0 0 4978 8 7451 0 0 0 6441 0 0 0 0 0 0 0 0 0
A 6442 2 0 0 4863 8 575 0 0 0 6442 0 0 0 0 0 0 0 0 0
A 6443 2 0 0 6326 8 7452 0 0 0 6443 0 0 0 0 0 0 0 0 0
A 6444 2 0 0 3489 9 7453 0 0 0 6444 0 0 0 0 0 0 0 0 0
A 6445 2 0 0 5466 8 7454 0 0 0 6445 0 0 0 0 0 0 0 0 0
A 6446 2 0 0 6222 8 7455 0 0 0 6446 0 0 0 0 0 0 0 0 0
A 6447 2 0 0 4433 8 574 0 0 0 6447 0 0 0 0 0 0 0 0 0
A 6448 2 0 0 4637 9 580 0 0 0 6448 0 0 0 0 0 0 0 0 0
A 6449 2 0 0 5275 8 7456 0 0 0 6449 0 0 0 0 0 0 0 0 0
A 6450 2 0 0 6126 9 7457 0 0 0 6450 0 0 0 0 0 0 0 0 0
A 6451 2 0 0 6071 9 7458 0 0 0 6451 0 0 0 0 0 0 0 0 0
A 6452 2 0 0 5471 9 7459 0 0 0 6452 0 0 0 0 0 0 0 0 0
A 6453 2 0 0 5991 6 7460 0 0 0 6453 0 0 0 0 0 0 0 0 0
A 6454 2 0 0 3497 9 7461 0 0 0 6454 0 0 0 0 0 0 0 0 0
A 6455 2 0 0 5459 9 7462 0 0 0 6455 0 0 0 0 0 0 0 0 0
A 6456 2 0 0 6131 9 577 0 0 0 6456 0 0 0 0 0 0 0 0 0
A 6457 2 0 0 4547 9 578 0 0 0 6457 0 0 0 0 0 0 0 0 0
A 6458 2 0 0 5568 9 7463 0 0 0 6458 0 0 0 0 0 0 0 0 0
A 6459 2 0 0 4991 9 7464 0 0 0 6459 0 0 0 0 0 0 0 0 0
A 6460 2 0 0 5998 9 7465 0 0 0 6460 0 0 0 0 0 0 0 0 0
A 6461 2 0 0 3404 9 7466 0 0 0 6461 0 0 0 0 0 0 0 0 0
A 6462 2 0 0 5479 9 7467 0 0 0 6462 0 0 0 0 0 0 0 0 0
A 6463 2 0 0 6433 9 7468 0 0 0 6463 0 0 0 0 0 0 0 0 0
A 6464 2 0 0 6434 9 7469 0 0 0 6464 0 0 0 0 0 0 0 0 0
A 6465 2 0 0 6435 9 7470 0 0 0 6465 0 0 0 0 0 0 0 0 0
A 6466 2 0 0 6436 9 7471 0 0 0 6466 0 0 0 0 0 0 0 0 0
A 6467 2 0 0 6437 9 7472 0 0 0 6467 0 0 0 0 0 0 0 0 0
A 6468 2 0 0 6438 8 7473 0 0 0 6468 0 0 0 0 0 0 0 0 0
A 6469 1 0 0 4873 6 9022 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6470 1 0 0 5932 6 9020 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6471 1 0 0 5414 6 9023 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6472 1 0 0 5905 6 9021 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6473 1 0 0 5012 6 9026 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6474 1 0 0 6421 6 9024 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6475 1 0 0 6019 6 9027 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6476 1 0 0 6007 6 9025 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 73 1 1
V 66 68 7 0
R 0 65 0 0
A 0 63 0 0 1 43 1
A 0 63 0 0 1 44 1
A 0 63 0 0 1 44 1
A 0 63 0 0 1 45 1
A 0 63 0 0 1 46 1
A 0 63 0 0 1 47 1
A 0 63 0 0 1 48 0
T 674 77 0 3 0 0
A 681 7 101 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 103 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 105 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
T 842 187 0 3 0 0
A 846 7 205 0 1 2 1
A 847 7 0 0 1 2 1
A 845 6 0 69 1 2 1
A 852 7 207 0 1 2 1
A 853 7 0 0 1 2 1
A 851 6 0 69 1 2 0
T 856 209 0 3 0 0
A 862 7 218 0 1 2 0
T 1047 395 0 3 0 0
T 1049 363 0 3 0 0
A 862 7 369 0 1 2 0
T 1052 401 0 3 0 0
A 1055 7 416 0 1 2 1
A 1059 7 418 0 1 2 1
A 1063 7 420 0 1 2 0
T 6584 2070 0 3 0 0
A 6585 16 0 0 1 3514 1
A 6586 16 0 0 1 3515 1
A 6587 8 0 0 1 3516 1
A 6588 20 0 0 1 3517 1
A 6590 6 0 0 1 2 1
A 6591 6 0 0 1 2 1
A 6592 6 0 0 1 2 1
A 6596 16 0 0 1 3515 0
T 6598 2080 0 3 0 0
A 6599 2066 0 0 1 3518 1
A 6600 2066 0 0 1 3518 1
A 6601 2066 0 0 1 3518 1
A 6602 2066 0 0 1 3518 1
A 6603 2066 0 0 1 3518 1
A 6604 2066 0 0 1 3518 1
A 6605 2066 0 0 1 3518 0
T 6606 2090 0 3 0 0
T 6607 2070 0 3 0 1
A 6585 16 0 0 1 3514 1
A 6586 16 0 0 1 3515 1
A 6587 8 0 0 1 3516 1
A 6588 20 0 0 1 3517 1
A 6590 6 0 0 1 2 1
A 6591 6 0 0 1 2 1
A 6592 6 0 0 1 2 1
A 6596 16 0 0 1 3515 0
A 6608 6 0 0 1 68 1
A 6609 6 0 0 1 2 1
A 6610 8 0 0 1 3503 1
A 6611 6 0 0 1 3 1
A 6612 8 0 0 1 3516 1
A 6613 6 0 0 1 39 1
A 6614 16 0 0 1 3515 1
T 6615 2080 0 3 0 1
A 6599 2066 0 0 1 3518 1
A 6600 2066 0 0 1 3518 1
A 6601 2066 0 0 1 3518 1
A 6602 2066 0 0 1 3518 1
A 6603 2066 0 0 1 3518 1
A 6604 2066 0 0 1 3518 1
A 6605 2066 0 0 1 3518 0
A 6618 7 2102 0 1 2 1
A 6622 7 2104 0 1 2 1
A 6624 16 0 0 1 3515 0
T 6625 2106 0 3 0 0
T 6626 2070 0 3 0 1
A 6585 16 0 0 1 3514 1
A 6586 16 0 0 1 3515 1
A 6587 8 0 0 1 3516 1
A 6588 20 0 0 1 3517 1
A 6590 6 0 0 1 2 1
A 6591 6 0 0 1 2 1
A 6592 6 0 0 1 2 1
A 6596 16 0 0 1 3515 0
A 6630 7 2124 0 1 2 1
A 6631 7 0 0 1 2 1
A 6629 6 0 69 1 2 1
A 6634 6 0 0 1 3 1
A 6635 6 0 0 1 3 1
A 6638 7 2126 0 1 2 1
A 6642 7 2128 0 1 2 0
T 6706 2144 0 3 0 0
A 6707 8 0 0 1 3574 0
T 6722 2156 0 3 0 0
A 6723 16 0 0 1 3515 1
A 6724 6 0 0 1 2 1
A 6725 6 0 0 1 2 1
A 6726 8 0 0 1 3575 1
A 6727 8 0 0 1 3576 1
A 6729 16 0 0 1 3515 1
T 6730 2144 0 3 0 1
A 6707 8 0 0 1 3574 0
A 6735 7 2180 0 1 2 1
A 6736 7 0 0 1 2 1
A 6734 6 0 3519 1 2 1
A 6742 7 2182 0 1 2 1
A 6743 7 0 0 1 2 1
A 6741 6 0 3519 1 2 1
A 6749 7 2184 0 1 2 1
A 6750 7 0 0 1 2 1
A 6748 6 0 3519 1 2 0
T 6767 2192 0 3 0 0
A 6768 6 0 0 1 2 1
A 6770 6 0 0 1 2 0
T 7073 2305 0 3 0 0
A 7074 6 0 0 1 2 1
A 7075 6 0 0 1 2 1
A 7076 6 0 0 1 3 0
T 7363 2446 0 3 0 0
A 7364 16 0 0 1 3515 1
A 7368 7 2488 0 1 2 1
A 7373 7 2490 0 1 2 1
A 7378 7 2492 0 1 2 1
A 7383 7 2494 0 1 2 0
T 7380 2479 0 3 0 0
A 7413 8 0 0 1 3516 1
A 7414 8 0 0 1 3516 0
T 7475 2723 0 3 0 0
A 7476 6 0 0 1 2 1
A 7477 6 0 0 1 2 1
A 7478 6 0 0 1 3 1
A 7484 7 2753 0 1 2 1
A 7485 7 0 0 1 2 1
A 7483 6 0 69 1 2 1
A 7490 7 2755 0 1 2 1
A 7491 7 0 0 1 2 1
A 7489 6 0 69 1 2 1
A 7496 7 2757 0 1 2 1
A 7497 7 0 0 1 2 1
A 7495 6 0 69 1 2 1
A 7502 7 2759 0 1 2 1
A 7503 7 0 0 1 2 1
A 7501 6 0 69 1 2 0
T 7545 2761 0 3 0 0
A 7546 6 0 0 1 3 1
A 7547 6 0 0 1 2 1
A 7548 6 0 0 1 2 1
A 7549 6 0 0 1 12 1
A 7550 6 0 0 1 3 1
A 7551 6 0 0 1 2 1
A 7552 6 0 0 1 2 1
A 7553 6 0 0 1 2 1
A 7554 6 0 0 1 2 1
A 7555 6 0 0 1 2 1
A 7556 6 0 0 1 2 1
A 7557 6 0 0 1 2 1
A 7558 6 0 0 1 14 1
A 7559 6 0 0 1 2 1
A 7560 6 0 0 1 36 1
A 7561 6 0 0 1 2 1
A 7562 6 0 0 1 2 1
A 7563 6 0 0 1 2 1
A 7564 6 0 0 1 2 1
A 7565 6 0 0 1 281 1
A 7566 6 0 0 1 2 1
A 7567 16 0 0 1 3514 1
A 7568 16 0 0 1 3515 1
A 7569 16 0 0 1 3515 1
A 7570 6 0 0 1 2 1
A 7571 6 0 0 1 2 1
A 7572 6 0 0 1 2 1
A 7573 6 0 0 1 2 1
A 7574 6 0 0 1 2 1
A 7575 6 0 0 1 2 1
A 7576 2613 0 0 1 6439 0
T 7577 2769 0 3 0 0
A 7582 7 2838 0 1 2 1
A 7583 7 0 0 1 2 1
A 7581 6 0 3519 1 2 1
A 7589 7 2840 0 1 2 1
A 7590 7 0 0 1 2 1
A 7588 6 0 3519 1 2 1
A 7596 7 2842 0 1 2 1
A 7597 7 0 0 1 2 1
A 7595 6 0 3519 1 2 1
A 7603 7 2844 0 1 2 1
A 7604 7 0 0 1 2 1
A 7602 6 0 3519 1 2 1
A 7610 7 2846 0 1 2 1
A 7611 7 0 0 1 2 1
A 7609 6 0 3519 1 2 1
A 7617 7 2848 0 1 2 1
A 7618 7 0 0 1 2 1
A 7616 6 0 3519 1 2 1
A 7625 7 2850 0 1 2 1
A 7626 7 0 0 1 2 1
A 7624 6 0 3631 1 2 1
A 7632 7 2852 0 1 2 1
A 7633 7 0 0 1 2 1
A 7631 6 0 3519 1 2 1
A 7639 7 2854 0 1 2 1
A 7640 7 0 0 1 2 1
A 7638 6 0 3519 1 2 1
A 7646 7 2856 0 1 2 1
A 7647 7 0 0 1 2 1
A 7645 6 0 3519 1 2 1
A 7649 6 0 0 1 2 1
R 7650 2835 0 1
A 0 6 0 14 1 2 0
A 7651 16 0 0 1 3514 0
T 7656 2858 0 3 0 0
A 7661 7 3020 0 1 2 1
A 7662 7 0 0 1 2 1
A 7660 6 0 3519 1 2 1
A 7668 7 3022 0 1 2 1
A 7669 7 0 0 1 2 1
A 7667 6 0 3519 1 2 1
A 7675 7 3024 0 1 2 1
A 7676 7 0 0 1 2 1
A 7674 6 0 3519 1 2 1
A 7682 7 3026 0 1 2 1
A 7683 7 0 0 1 2 1
A 7681 6 0 3519 1 2 1
A 7689 7 3028 0 1 2 1
A 7690 7 0 0 1 2 1
A 7688 6 0 3519 1 2 1
A 7696 7 3030 0 1 2 1
A 7697 7 0 0 1 2 1
A 7695 6 0 3519 1 2 1
A 7703 7 3032 0 1 2 1
A 7704 7 0 0 1 2 1
A 7702 6 0 3519 1 2 1
A 7710 7 3034 0 1 2 1
A 7711 7 0 0 1 2 1
A 7709 6 0 3519 1 2 1
A 7717 7 3036 0 1 2 1
A 7718 7 0 0 1 2 1
A 7716 6 0 3519 1 2 1
A 7724 7 3038 0 1 2 1
A 7725 7 0 0 1 2 1
A 7723 6 0 3519 1 2 1
A 7731 7 3040 0 1 2 1
A 7732 7 0 0 1 2 1
A 7730 6 0 3519 1 2 1
A 7738 7 3042 0 1 2 1
A 7739 7 0 0 1 2 1
A 7737 6 0 3519 1 2 1
A 7745 7 3044 0 1 2 1
A 7746 7 0 0 1 2 1
A 7744 6 0 3519 1 2 1
A 7752 7 3046 0 1 2 1
A 7753 7 0 0 1 2 1
A 7751 6 0 3519 1 2 1
A 7759 7 3048 0 1 2 1
A 7760 7 0 0 1 2 1
A 7758 6 0 3519 1 2 1
A 7766 7 3050 0 1 2 1
A 7767 7 0 0 1 2 1
A 7765 6 0 3519 1 2 1
A 7773 7 3052 0 1 2 1
A 7774 7 0 0 1 2 1
A 7772 6 0 3519 1 2 1
A 7780 7 3054 0 1 2 1
A 7781 7 0 0 1 2 1
A 7779 6 0 3519 1 2 1
A 7787 7 3056 0 1 2 1
A 7788 7 0 0 1 2 1
A 7786 6 0 3519 1 2 1
A 7794 7 3058 0 1 2 1
A 7795 7 0 0 1 2 1
A 7793 6 0 3519 1 2 1
A 7801 7 3060 0 1 2 1
A 7802 7 0 0 1 2 1
A 7800 6 0 3519 1 2 1
A 7808 7 3062 0 1 2 1
A 7809 7 0 0 1 2 1
A 7807 6 0 3519 1 2 1
A 7815 7 3064 0 1 2 1
A 7816 7 0 0 1 2 1
A 7814 6 0 3519 1 2 1
A 7822 7 3066 0 1 2 1
A 7823 7 0 0 1 2 1
A 7821 6 0 3519 1 2 1
A 7829 7 3068 0 1 2 1
A 7830 7 0 0 1 2 1
A 7828 6 0 3519 1 2 1
A 7836 7 3070 0 1 2 1
A 7837 7 0 0 1 2 1
A 7835 6 0 3519 1 2 0
T 7839 3072 0 3 0 0
A 7845 7 3114 0 1 2 1
A 7846 7 0 0 1 2 1
A 7844 6 0 3631 1 2 1
A 7853 7 3116 0 1 2 1
A 7854 7 0 0 1 2 1
A 7852 6 0 3631 1 2 1
A 7861 7 3118 0 1 2 1
A 7862 7 0 0 1 2 1
A 7860 6 0 3631 1 2 1
A 7869 7 3120 0 1 2 1
A 7870 7 0 0 1 2 1
A 7868 6 0 3631 1 2 1
A 7877 7 3122 0 1 2 1
A 7878 7 0 0 1 2 1
A 7876 6 0 3631 1 2 1
A 7885 7 3124 0 1 2 1
A 7886 7 0 0 1 2 1
A 7884 6 0 3631 1 2 0
T 7888 3126 0 3 0 0
A 7894 7 3270 0 1 2 1
A 7895 7 0 0 1 2 1
A 7893 6 0 3631 1 2 1
A 7902 7 3272 0 1 2 1
A 7903 7 0 0 1 2 1
A 7901 6 0 3631 1 2 1
A 7910 7 3274 0 1 2 1
A 7911 7 0 0 1 2 1
A 7909 6 0 3631 1 2 1
A 7918 7 3276 0 1 2 1
A 7919 7 0 0 1 2 1
A 7917 6 0 3631 1 2 1
A 7926 7 3278 0 1 2 1
A 7927 7 0 0 1 2 1
A 7925 6 0 3631 1 2 1
A 7934 7 3280 0 1 2 1
A 7935 7 0 0 1 2 1
A 7933 6 0 3631 1 2 1
A 7941 7 3282 0 1 2 1
A 7942 7 0 0 1 2 1
A 7940 6 0 3519 1 2 1
A 7948 7 3284 0 1 2 1
A 7949 7 0 0 1 2 1
A 7947 6 0 3519 1 2 1
A 7955 7 3286 0 1 2 1
A 7956 7 0 0 1 2 1
A 7954 6 0 3519 1 2 1
A 7962 7 3288 0 1 2 1
A 7963 7 0 0 1 2 1
A 7961 6 0 3519 1 2 1
A 7969 7 3290 0 1 2 1
A 7970 7 0 0 1 2 1
A 7968 6 0 3519 1 2 1
A 7976 7 3292 0 1 2 1
A 7977 7 0 0 1 2 1
A 7975 6 0 3519 1 2 1
A 7983 7 3294 0 1 2 1
A 7984 7 0 0 1 2 1
A 7982 6 0 3519 1 2 1
A 7990 7 3296 0 1 2 1
A 7991 7 0 0 1 2 1
A 7989 6 0 3519 1 2 1
A 7997 7 3298 0 1 2 1
A 7998 7 0 0 1 2 1
A 7996 6 0 3519 1 2 1
A 8004 7 3300 0 1 2 1
A 8005 7 0 0 1 2 1
A 8003 6 0 3519 1 2 1
A 8007 16 0 0 1 3515 1
A 8012 7 3302 0 1 2 1
A 8013 7 0 0 1 2 1
A 8011 6 0 3519 1 2 1
A 8019 7 3304 0 1 2 1
A 8020 7 0 0 1 2 1
A 8018 6 0 3519 1 2 1
A 8027 7 3306 0 1 2 1
A 8028 7 0 0 1 2 1
A 8026 6 0 3631 1 2 1
A 8034 7 3308 0 1 2 1
A 8035 7 0 0 1 2 1
A 8033 6 0 3519 1 2 1
A 8041 7 3310 0 1 2 1
A 8042 7 0 0 1 2 1
A 8040 6 0 3519 1 2 1
A 8048 7 3312 0 1 2 1
A 8049 7 0 0 1 2 1
A 8047 6 0 3519 1 2 1
A 8055 7 3314 0 1 2 1
A 8056 7 0 0 1 2 1
A 8054 6 0 3519 1 2 0
T 8058 3316 0 3 0 0
T 8059 3072 0 3 0 1
A 7845 7 3114 0 1 2 1
A 7846 7 0 0 1 2 1
A 7844 6 0 3631 1 2 1
A 7853 7 3116 0 1 2 1
A 7854 7 0 0 1 2 1
A 7852 6 0 3631 1 2 1
A 7861 7 3118 0 1 2 1
A 7862 7 0 0 1 2 1
A 7860 6 0 3631 1 2 1
A 7869 7 3120 0 1 2 1
A 7870 7 0 0 1 2 1
A 7868 6 0 3631 1 2 1
A 7877 7 3122 0 1 2 1
A 7878 7 0 0 1 2 1
A 7876 6 0 3631 1 2 1
A 7885 7 3124 0 1 2 1
A 7886 7 0 0 1 2 1
A 7884 6 0 3631 1 2 0
A 8065 7 3340 0 1 2 1
A 8066 7 0 0 1 2 1
A 8064 6 0 3631 1 2 1
A 8072 7 3342 0 1 2 1
A 8073 7 0 0 1 2 1
A 8071 6 0 3519 1 2 1
A 8079 7 3344 0 1 2 1
A 8080 7 0 0 1 2 1
A 8078 6 0 3519 1 2 0
T 8082 3346 0 3 0 0
A 8087 7 3400 0 1 2 1
A 8088 7 0 0 1 2 1
A 8086 6 0 3519 1 2 1
A 8094 7 3402 0 1 2 1
A 8095 7 0 0 1 2 1
A 8093 6 0 3519 1 2 1
A 8101 7 3404 0 1 2 1
A 8102 7 0 0 1 2 1
A 8100 6 0 3519 1 2 1
A 8108 7 3406 0 1 2 1
A 8109 7 0 0 1 2 1
A 8107 6 0 3519 1 2 1
A 8115 7 3408 0 1 2 1
A 8116 7 0 0 1 2 1
A 8114 6 0 3519 1 2 1
A 8122 7 3410 0 1 2 1
A 8123 7 0 0 1 2 1
A 8121 6 0 3519 1 2 1
A 8125 8 0 0 1 3516 1
A 8130 7 3412 0 1 2 1
A 8131 7 0 0 1 2 1
A 8129 6 0 3519 1 2 1
A 8137 7 3414 0 1 2 1
A 8138 7 0 0 1 2 1
A 8136 6 0 3519 1 2 1
A 8140 8 0 0 1 6440 1
A 8141 8 0 0 1 6441 1
A 8142 8 0 0 1 3516 0
T 8144 3416 0 3 0 0
A 8150 7 3512 0 1 2 1
A 8151 7 0 0 1 2 1
A 8149 6 0 3631 1 2 1
A 8157 7 3514 0 1 2 1
A 8158 7 0 0 1 2 1
A 8156 6 0 3519 1 2 1
A 8165 7 3516 0 1 2 1
A 8166 7 0 0 1 2 1
A 8164 6 0 3631 1 2 1
A 8172 7 3518 0 1 2 1
A 8173 7 0 0 1 2 1
A 8171 6 0 3519 1 2 1
A 8179 7 3520 0 1 2 1
A 8180 7 0 0 1 2 1
A 8178 6 0 3519 1 2 1
A 8186 7 3522 0 1 2 1
A 8187 7 0 0 1 2 1
A 8185 6 0 3519 1 2 1
A 8193 7 3524 0 1 2 1
A 8194 7 0 0 1 2 1
A 8192 6 0 3519 1 2 1
A 8200 7 3526 0 1 2 1
A 8201 7 0 0 1 2 1
A 8199 6 0 3519 1 2 1
A 8207 7 3528 0 1 2 1
A 8208 7 0 0 1 2 1
A 8206 6 0 3519 1 2 1
A 8214 7 3530 0 1 2 1
A 8215 7 0 0 1 2 1
A 8213 6 0 3519 1 2 1
A 8221 7 3532 0 1 2 1
A 8222 7 0 0 1 2 1
A 8220 6 0 3519 1 2 1
A 8228 7 3534 0 1 2 1
A 8229 7 0 0 1 2 1
A 8227 6 0 3519 1 2 1
A 8235 7 3536 0 1 2 1
A 8236 7 0 0 1 2 1
A 8234 6 0 3519 1 2 1
A 8242 7 3538 0 1 2 1
A 8243 7 0 0 1 2 1
A 8241 6 0 3519 1 2 1
A 8249 7 3540 0 1 2 1
A 8250 7 0 0 1 2 1
A 8248 6 0 3519 1 2 1
A 8252 6 0 0 1 2 1
A 8253 8 0 0 1 3516 1
A 8254 8 0 0 1 3516 1
A 8255 6 0 0 1 2 1
A 8256 16 0 0 1 3514 1
A 8257 16 0 0 1 3515 0
T 8258 3542 0 3 0 0
A 8264 7 3614 0 1 2 1
A 8265 7 0 0 1 2 1
A 8263 6 0 3631 1 2 1
A 8271 7 3616 0 1 2 1
A 8272 7 0 0 1 2 1
A 8270 6 0 3519 1 2 1
A 8274 6 0 0 1 3 1
T 8275 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
T 8276 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
A 8324 7 3618 0 1 2 1
A 8325 7 0 0 1 2 1
A 8323 6 0 69 1 2 1
A 8331 7 3620 0 1 2 1
A 8332 7 0 0 1 2 1
A 8330 6 0 3519 1 2 1
A 8336 8 0 0 1 6442 1
A 8337 8 0 0 1 6442 1
A 8338 6 0 0 1 183 1
A 8339 8 0 0 1 6443 1
A 8340 6 0 0 1 2 1
A 8341 9 0 0 1 3732 1
A 8342 9 0 0 1 3736 1
A 8343 9 0 0 1 6444 1
A 8344 8 0 0 1 3516 0
T 8345 3622 0 3 0 0
A 8346 2613 0 0 1 6439 1
A 8347 2613 0 0 1 6439 1
A 8350 7 3634 0 1 2 1
A 8354 7 3636 0 1 2 0
T 8356 3638 0 3 0 0
A 8357 8 0 0 1 3516 1
A 8358 8 0 0 1 6445 1
A 8359 8 0 0 1 3516 1
A 8360 8 0 0 1 6446 1
A 8361 8 0 0 1 6447 1
A 8362 8 0 0 1 6447 1
A 8363 9 0 0 1 6448 1
A 8364 9 0 0 1 6448 1
A 8365 8 0 0 1 6449 1
A 8366 9 0 0 1 6450 1
A 8367 9 0 0 1 6451 1
A 8368 9 0 0 1 6452 1
A 8369 9 0 0 1 6452 1
A 8370 8 0 0 1 3516 1
A 8371 8 0 0 1 3516 1
A 8372 8 0 0 1 3516 1
A 8373 6 0 0 1 3 1
A 8374 6 0 0 1 2 1
A 8378 7 3662 0 1 2 1
A 8379 7 0 0 1 2 1
A 8377 6 0 69 1 2 1
A 8384 7 3664 0 1 2 1
A 8385 7 0 0 1 2 1
A 8383 6 0 69 1 2 1
A 8390 7 3666 0 1 2 1
A 8391 7 0 0 1 2 1
A 8389 6 0 69 1 2 1
A 8393 6 0 0 1 3499 1
A 8394 6 0 0 1 6453 1
A 8395 6 0 0 1 3 1
A 8396 6 0 0 1 3 0
T 8397 3668 0 3 0 0
A 8402 7 3692 0 1 2 1
A 8403 7 0 0 1 2 1
A 8401 6 0 3519 1 2 1
A 8409 7 3694 0 1 2 1
A 8410 7 0 0 1 2 1
A 8408 6 0 3519 1 2 1
A 8416 7 3696 0 1 2 1
A 8417 7 0 0 1 2 1
A 8415 6 0 3519 1 2 0
T 8419 3698 0 3 0 0
A 8423 7 3758 0 1 2 1
A 8424 7 0 0 1 2 1
A 8422 6 0 69 1 2 1
A 8429 7 3760 0 1 2 1
A 8430 7 0 0 1 2 1
A 8428 6 0 69 1 2 1
A 8435 7 3762 0 1 2 1
A 8436 7 0 0 1 2 1
A 8434 6 0 69 1 2 1
A 8441 7 3764 0 1 2 1
A 8442 7 0 0 1 2 1
A 8440 6 0 69 1 2 1
A 8447 7 3766 0 1 2 1
A 8448 7 0 0 1 2 1
A 8446 6 0 69 1 2 1
A 8454 7 3768 0 1 2 1
A 8455 7 0 0 1 2 1
A 8453 6 0 3519 1 2 1
A 8461 7 3770 0 1 2 1
A 8462 7 0 0 1 2 1
A 8460 6 0 3519 1 2 1
A 8467 7 3772 0 1 2 1
A 8468 7 0 0 1 2 1
A 8466 6 0 69 1 2 1
R 8471 3755 0 1
A 0 8 0 14 1 3516 0
A 8472 9 0 0 1 6454 1
A 8473 9 0 0 1 6455 1
A 8474 8 0 0 1 3516 1
A 8475 9 0 0 1 6456 1
A 8476 9 0 0 1 6457 1
A 8477 9 0 0 1 6455 1
A 8478 9 0 0 1 6456 1
A 8479 9 0 0 1 6456 1
A 8480 9 0 0 1 6456 0
T 8481 3774 0 3 0 0
T 8482 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
A 8486 7 3798 0 1 2 1
A 8487 7 0 0 1 2 1
A 8485 6 0 69 1 2 1
A 8492 7 3800 0 1 2 1
A 8493 7 0 0 1 2 1
A 8491 6 0 69 1 2 1
R 8495 3792 0 1
A 0 8 0 14 1 3516 0
R 8496 3795 0 1
A 0 8 0 41 1 3516 0
A 8497 6 0 0 1 2 0
T 8498 3802 0 3 0 0
A 8503 7 3856 0 1 2 1
A 8504 7 0 0 1 2 1
A 8502 6 0 3519 1 2 1
A 8510 7 3858 0 1 2 1
A 8511 7 0 0 1 2 1
A 8509 6 0 3519 1 2 1
A 8517 7 3860 0 1 2 1
A 8518 7 0 0 1 2 1
A 8516 6 0 3519 1 2 1
A 8525 7 3862 0 1 2 1
A 8526 7 0 0 1 2 1
A 8524 6 0 3631 1 2 1
A 8528 6 0 0 1 12 1
A 8529 8 0 0 1 3516 1
A 8533 7 3864 0 1 2 1
A 8534 7 0 0 1 2 1
A 8532 6 0 69 1 2 1
A 8539 7 3866 0 1 2 1
A 8540 7 0 0 1 2 1
A 8538 6 0 69 1 2 1
A 8545 7 3868 0 1 2 1
A 8546 7 0 0 1 2 1
A 8544 6 0 69 1 2 1
A 8551 7 3870 0 1 2 1
A 8552 7 0 0 1 2 1
A 8550 6 0 69 1 2 0
T 8554 3872 0 3 0 0
A 8560 7 4001 0 1 2 1
A 8561 7 0 0 1 2 1
A 8559 6 0 3631 1 2 1
A 8568 7 4003 0 1 2 1
A 8569 7 0 0 1 2 1
A 8567 6 0 3631 1 2 1
A 8576 7 4005 0 1 2 1
A 8577 7 0 0 1 2 1
A 8575 6 0 3631 1 2 1
A 8584 7 4007 0 1 2 1
A 8585 7 0 0 1 2 1
A 8583 6 0 3631 1 2 1
A 8590 7 4009 0 1 2 1
A 8591 7 0 0 1 2 1
A 8589 6 0 69 1 2 1
A 8596 7 4011 0 1 2 1
A 8597 7 0 0 1 2 1
A 8595 6 0 69 1 2 1
A 8602 7 4013 0 1 2 1
A 8603 7 0 0 1 2 1
A 8601 6 0 69 1 2 1
A 8608 7 4015 0 1 2 1
A 8609 7 0 0 1 2 1
A 8607 6 0 69 1 2 1
A 8615 7 4017 0 1 2 1
A 8616 7 0 0 1 2 1
A 8614 6 0 3519 1 2 1
A 8622 7 4019 0 1 2 1
A 8623 7 0 0 1 2 1
A 8621 6 0 3519 1 2 1
A 8629 7 4021 0 1 2 1
A 8630 7 0 0 1 2 1
A 8628 6 0 3519 1 2 1
A 8636 7 4023 0 1 2 1
A 8637 7 0 0 1 2 1
A 8635 6 0 3519 1 2 1
A 8643 7 4025 0 1 2 1
A 8644 7 0 0 1 2 1
A 8642 6 0 3519 1 2 1
A 8650 7 4027 0 1 2 1
A 8651 7 0 0 1 2 1
A 8649 6 0 3519 1 2 1
A 8657 7 4029 0 1 2 1
A 8658 7 0 0 1 2 1
A 8656 6 0 3519 1 2 1
A 8664 7 4031 0 1 2 1
A 8665 7 0 0 1 2 1
A 8663 6 0 3519 1 2 1
A 8672 7 4033 0 1 2 1
A 8673 7 0 0 1 2 1
A 8671 6 0 3631 1 2 1
A 8680 7 4035 0 1 2 1
A 8681 7 0 0 1 2 1
A 8679 6 0 3631 1 2 1
R 8683 3986 0 1
A 0 8 0 39 1 3516 0
R 8684 3989 0 1
A 0 8 0 39 1 3516 0
R 8685 3992 0 1
A 0 8 0 16 1 3516 0
A 8687 6 0 0 1 281 1
R 8688 3998 0 1
A 0 8 0 12 1 3516 0
A 8689 8 0 0 1 3516 1
A 8690 8 0 0 1 3516 1
A 8691 8 0 0 1 3516 1
A 8692 8 0 0 1 3516 1
A 8693 8 0 0 1 3516 1
A 8694 6 0 0 1 2 0
T 8695 4037 0 3 0 0
A 8700 7 4115 0 1 2 1
A 8701 7 0 0 1 2 1
A 8699 6 0 3519 1 2 1
A 8707 7 4117 0 1 2 1
A 8708 7 0 0 1 2 1
A 8706 6 0 3519 1 2 1
A 8714 7 4119 0 1 2 1
A 8715 7 0 0 1 2 1
A 8713 6 0 3519 1 2 1
A 8721 7 4121 0 1 2 1
A 8722 7 0 0 1 2 1
A 8720 6 0 3519 1 2 1
A 8728 7 4123 0 1 2 1
A 8729 7 0 0 1 2 1
A 8727 6 0 3519 1 2 1
A 8735 7 4125 0 1 2 1
A 8736 7 0 0 1 2 1
A 8734 6 0 3519 1 2 1
A 8742 7 4127 0 1 2 1
A 8743 7 0 0 1 2 1
A 8741 6 0 3519 1 2 1
A 8749 7 4129 0 1 2 1
A 8750 7 0 0 1 2 1
A 8748 6 0 3519 1 2 1
A 8756 7 4131 0 1 2 1
A 8757 7 0 0 1 2 1
A 8755 6 0 3519 1 2 1
A 8763 7 4133 0 1 2 1
A 8764 7 0 0 1 2 1
A 8762 6 0 3519 1 2 1
A 8770 7 4135 0 1 2 1
A 8771 7 0 0 1 2 1
A 8769 6 0 3519 1 2 1
A 8777 7 4137 0 1 2 1
A 8778 7 0 0 1 2 1
A 8776 6 0 3519 1 2 0
T 8780 4139 0 3 0 0
R 8781 4148 0 1
A 0 9 0 0 1 6458 1
A 0 9 0 0 1 6448 1
A 0 9 0 0 1 6456 1
A 0 9 0 0 1 6459 1
A 0 9 0 0 1 6457 0
A 8782 9 0 0 1 6456 1
A 8783 9 0 0 1 6456 1
A 8784 9 0 0 1 6456 1
A 8785 9 0 0 1 6460 1
A 8786 9 0 0 1 6454 1
A 8787 9 0 0 1 6457 1
A 8788 9 0 0 1 3736 1
A 8789 9 0 0 1 6459 1
A 8790 9 0 0 1 6461 0
T 8791 4151 0 3 0 0
A 8792 9 0 0 1 6462 1
A 8793 9 0 0 1 6463 1
A 8794 9 0 0 1 3736 1
A 8795 9 0 0 1 6464 1
A 8796 9 0 0 1 6465 1
A 8797 9 0 0 1 6466 1
A 8798 9 0 0 1 6454 1
A 8799 9 0 0 1 6467 1
A 8800 8 0 0 1 6468 1
A 8801 6 0 0 1 2 1
A 8802 6 0 0 1 39 1
A 8806 7 4187 0 1 2 1
A 8807 7 0 0 1 2 1
A 8805 6 0 69 1 2 1
A 8813 7 4189 0 1 2 1
A 8814 7 0 0 1 2 1
A 8812 6 0 3519 1 2 1
A 8820 7 4191 0 1 2 1
A 8821 7 0 0 1 2 1
A 8819 6 0 3519 1 2 1
A 8828 7 4193 0 1 2 1
A 8829 7 0 0 1 2 1
A 8827 6 0 3631 1 2 1
A 8836 7 4195 0 1 2 1
A 8837 7 0 0 1 2 1
A 8835 6 0 3631 1 2 0
T 8848 4203 0 3 0 0
A 8853 7 4227 0 1 2 1
A 8854 7 0 0 1 2 1
A 8852 6 0 3519 1 2 1
A 8860 7 4229 0 1 2 1
A 8861 7 0 0 1 2 1
A 8859 6 0 3519 1 2 1
A 8867 7 4231 0 1 2 1
A 8868 7 0 0 1 2 1
A 8866 6 0 3519 1 2 0
T 8870 4233 0 3 0 0
T 8872 2723 0 3 0 1
A 7476 6 0 0 1 2 1
A 7477 6 0 0 1 2 1
A 7478 6 0 0 1 3 1
A 7484 7 2753 0 1 2 1
A 7485 7 0 0 1 2 1
A 7483 6 0 69 1 2 1
A 7490 7 2755 0 1 2 1
A 7491 7 0 0 1 2 1
A 7489 6 0 69 1 2 1
A 7496 7 2757 0 1 2 1
A 7497 7 0 0 1 2 1
A 7495 6 0 69 1 2 1
A 7502 7 2759 0 1 2 1
A 7503 7 0 0 1 2 1
A 7501 6 0 69 1 2 0
T 8873 2761 0 3 0 1
A 7546 6 0 0 1 3 1
A 7547 6 0 0 1 2 1
A 7548 6 0 0 1 2 1
A 7549 6 0 0 1 12 1
A 7550 6 0 0 1 3 1
A 7551 6 0 0 1 2 1
A 7552 6 0 0 1 2 1
A 7553 6 0 0 1 2 1
A 7554 6 0 0 1 2 1
A 7555 6 0 0 1 2 1
A 7556 6 0 0 1 2 1
A 7557 6 0 0 1 2 1
A 7558 6 0 0 1 14 1
A 7559 6 0 0 1 2 1
A 7560 6 0 0 1 36 1
A 7561 6 0 0 1 2 1
A 7562 6 0 0 1 2 1
A 7563 6 0 0 1 2 1
A 7564 6 0 0 1 2 1
A 7565 6 0 0 1 281 1
A 7566 6 0 0 1 2 1
A 7567 16 0 0 1 3514 1
A 7568 16 0 0 1 3515 1
A 7569 16 0 0 1 3515 1
A 7570 6 0 0 1 2 1
A 7571 6 0 0 1 2 1
A 7572 6 0 0 1 2 1
A 7573 6 0 0 1 2 1
A 7574 6 0 0 1 2 1
A 7575 6 0 0 1 2 1
A 7576 2613 0 0 1 6439 0
T 8874 2769 0 3 0 1
A 7582 7 2838 0 1 2 1
A 7583 7 0 0 1 2 1
A 7581 6 0 3519 1 2 1
A 7589 7 2840 0 1 2 1
A 7590 7 0 0 1 2 1
A 7588 6 0 3519 1 2 1
A 7596 7 2842 0 1 2 1
A 7597 7 0 0 1 2 1
A 7595 6 0 3519 1 2 1
A 7603 7 2844 0 1 2 1
A 7604 7 0 0 1 2 1
A 7602 6 0 3519 1 2 1
A 7610 7 2846 0 1 2 1
A 7611 7 0 0 1 2 1
A 7609 6 0 3519 1 2 1
A 7617 7 2848 0 1 2 1
A 7618 7 0 0 1 2 1
A 7616 6 0 3519 1 2 1
A 7625 7 2850 0 1 2 1
A 7626 7 0 0 1 2 1
A 7624 6 0 3631 1 2 1
A 7632 7 2852 0 1 2 1
A 7633 7 0 0 1 2 1
A 7631 6 0 3519 1 2 1
A 7639 7 2854 0 1 2 1
A 7640 7 0 0 1 2 1
A 7638 6 0 3519 1 2 1
A 7646 7 2856 0 1 2 1
A 7647 7 0 0 1 2 1
A 7645 6 0 3519 1 2 1
A 7649 6 0 0 1 2 1
R 7650 2835 0 1
A 0 6 0 14 1 2 0
A 7651 16 0 0 1 3514 0
T 8875 2858 0 3 0 1
A 7661 7 3020 0 1 2 1
A 7662 7 0 0 1 2 1
A 7660 6 0 3519 1 2 1
A 7668 7 3022 0 1 2 1
A 7669 7 0 0 1 2 1
A 7667 6 0 3519 1 2 1
A 7675 7 3024 0 1 2 1
A 7676 7 0 0 1 2 1
A 7674 6 0 3519 1 2 1
A 7682 7 3026 0 1 2 1
A 7683 7 0 0 1 2 1
A 7681 6 0 3519 1 2 1
A 7689 7 3028 0 1 2 1
A 7690 7 0 0 1 2 1
A 7688 6 0 3519 1 2 1
A 7696 7 3030 0 1 2 1
A 7697 7 0 0 1 2 1
A 7695 6 0 3519 1 2 1
A 7703 7 3032 0 1 2 1
A 7704 7 0 0 1 2 1
A 7702 6 0 3519 1 2 1
A 7710 7 3034 0 1 2 1
A 7711 7 0 0 1 2 1
A 7709 6 0 3519 1 2 1
A 7717 7 3036 0 1 2 1
A 7718 7 0 0 1 2 1
A 7716 6 0 3519 1 2 1
A 7724 7 3038 0 1 2 1
A 7725 7 0 0 1 2 1
A 7723 6 0 3519 1 2 1
A 7731 7 3040 0 1 2 1
A 7732 7 0 0 1 2 1
A 7730 6 0 3519 1 2 1
A 7738 7 3042 0 1 2 1
A 7739 7 0 0 1 2 1
A 7737 6 0 3519 1 2 1
A 7745 7 3044 0 1 2 1
A 7746 7 0 0 1 2 1
A 7744 6 0 3519 1 2 1
A 7752 7 3046 0 1 2 1
A 7753 7 0 0 1 2 1
A 7751 6 0 3519 1 2 1
A 7759 7 3048 0 1 2 1
A 7760 7 0 0 1 2 1
A 7758 6 0 3519 1 2 1
A 7766 7 3050 0 1 2 1
A 7767 7 0 0 1 2 1
A 7765 6 0 3519 1 2 1
A 7773 7 3052 0 1 2 1
A 7774 7 0 0 1 2 1
A 7772 6 0 3519 1 2 1
A 7780 7 3054 0 1 2 1
A 7781 7 0 0 1 2 1
A 7779 6 0 3519 1 2 1
A 7787 7 3056 0 1 2 1
A 7788 7 0 0 1 2 1
A 7786 6 0 3519 1 2 1
A 7794 7 3058 0 1 2 1
A 7795 7 0 0 1 2 1
A 7793 6 0 3519 1 2 1
A 7801 7 3060 0 1 2 1
A 7802 7 0 0 1 2 1
A 7800 6 0 3519 1 2 1
A 7808 7 3062 0 1 2 1
A 7809 7 0 0 1 2 1
A 7807 6 0 3519 1 2 1
A 7815 7 3064 0 1 2 1
A 7816 7 0 0 1 2 1
A 7814 6 0 3519 1 2 1
A 7822 7 3066 0 1 2 1
A 7823 7 0 0 1 2 1
A 7821 6 0 3519 1 2 1
A 7829 7 3068 0 1 2 1
A 7830 7 0 0 1 2 1
A 7828 6 0 3519 1 2 1
A 7836 7 3070 0 1 2 1
A 7837 7 0 0 1 2 1
A 7835 6 0 3519 1 2 0
T 8876 3126 0 3 0 1
A 7894 7 3270 0 1 2 1
A 7895 7 0 0 1 2 1
A 7893 6 0 3631 1 2 1
A 7902 7 3272 0 1 2 1
A 7903 7 0 0 1 2 1
A 7901 6 0 3631 1 2 1
A 7910 7 3274 0 1 2 1
A 7911 7 0 0 1 2 1
A 7909 6 0 3631 1 2 1
A 7918 7 3276 0 1 2 1
A 7919 7 0 0 1 2 1
A 7917 6 0 3631 1 2 1
A 7926 7 3278 0 1 2 1
A 7927 7 0 0 1 2 1
A 7925 6 0 3631 1 2 1
A 7934 7 3280 0 1 2 1
A 7935 7 0 0 1 2 1
A 7933 6 0 3631 1 2 1
A 7941 7 3282 0 1 2 1
A 7942 7 0 0 1 2 1
A 7940 6 0 3519 1 2 1
A 7948 7 3284 0 1 2 1
A 7949 7 0 0 1 2 1
A 7947 6 0 3519 1 2 1
A 7955 7 3286 0 1 2 1
A 7956 7 0 0 1 2 1
A 7954 6 0 3519 1 2 1
A 7962 7 3288 0 1 2 1
A 7963 7 0 0 1 2 1
A 7961 6 0 3519 1 2 1
A 7969 7 3290 0 1 2 1
A 7970 7 0 0 1 2 1
A 7968 6 0 3519 1 2 1
A 7976 7 3292 0 1 2 1
A 7977 7 0 0 1 2 1
A 7975 6 0 3519 1 2 1
A 7983 7 3294 0 1 2 1
A 7984 7 0 0 1 2 1
A 7982 6 0 3519 1 2 1
A 7990 7 3296 0 1 2 1
A 7991 7 0 0 1 2 1
A 7989 6 0 3519 1 2 1
A 7997 7 3298 0 1 2 1
A 7998 7 0 0 1 2 1
A 7996 6 0 3519 1 2 1
A 8004 7 3300 0 1 2 1
A 8005 7 0 0 1 2 1
A 8003 6 0 3519 1 2 1
A 8007 16 0 0 1 3515 1
A 8012 7 3302 0 1 2 1
A 8013 7 0 0 1 2 1
A 8011 6 0 3519 1 2 1
A 8019 7 3304 0 1 2 1
A 8020 7 0 0 1 2 1
A 8018 6 0 3519 1 2 1
A 8027 7 3306 0 1 2 1
A 8028 7 0 0 1 2 1
A 8026 6 0 3631 1 2 1
A 8034 7 3308 0 1 2 1
A 8035 7 0 0 1 2 1
A 8033 6 0 3519 1 2 1
A 8041 7 3310 0 1 2 1
A 8042 7 0 0 1 2 1
A 8040 6 0 3519 1 2 1
A 8048 7 3312 0 1 2 1
A 8049 7 0 0 1 2 1
A 8047 6 0 3519 1 2 1
A 8055 7 3314 0 1 2 1
A 8056 7 0 0 1 2 1
A 8054 6 0 3519 1 2 0
T 8877 3316 0 3 0 1
T 8059 3072 0 3 0 1
A 7845 7 3114 0 1 2 1
A 7846 7 0 0 1 2 1
A 7844 6 0 3631 1 2 1
A 7853 7 3116 0 1 2 1
A 7854 7 0 0 1 2 1
A 7852 6 0 3631 1 2 1
A 7861 7 3118 0 1 2 1
A 7862 7 0 0 1 2 1
A 7860 6 0 3631 1 2 1
A 7869 7 3120 0 1 2 1
A 7870 7 0 0 1 2 1
A 7868 6 0 3631 1 2 1
A 7877 7 3122 0 1 2 1
A 7878 7 0 0 1 2 1
A 7876 6 0 3631 1 2 1
A 7885 7 3124 0 1 2 1
A 7886 7 0 0 1 2 1
A 7884 6 0 3631 1 2 0
A 8065 7 3340 0 1 2 1
A 8066 7 0 0 1 2 1
A 8064 6 0 3631 1 2 1
A 8072 7 3342 0 1 2 1
A 8073 7 0 0 1 2 1
A 8071 6 0 3519 1 2 1
A 8079 7 3344 0 1 2 1
A 8080 7 0 0 1 2 1
A 8078 6 0 3519 1 2 0
T 8878 3346 0 3 0 1
A 8087 7 3400 0 1 2 1
A 8088 7 0 0 1 2 1
A 8086 6 0 3519 1 2 1
A 8094 7 3402 0 1 2 1
A 8095 7 0 0 1 2 1
A 8093 6 0 3519 1 2 1
A 8101 7 3404 0 1 2 1
A 8102 7 0 0 1 2 1
A 8100 6 0 3519 1 2 1
A 8108 7 3406 0 1 2 1
A 8109 7 0 0 1 2 1
A 8107 6 0 3519 1 2 1
A 8115 7 3408 0 1 2 1
A 8116 7 0 0 1 2 1
A 8114 6 0 3519 1 2 1
A 8122 7 3410 0 1 2 1
A 8123 7 0 0 1 2 1
A 8121 6 0 3519 1 2 1
A 8125 8 0 0 1 3516 1
A 8130 7 3412 0 1 2 1
A 8131 7 0 0 1 2 1
A 8129 6 0 3519 1 2 1
A 8137 7 3414 0 1 2 1
A 8138 7 0 0 1 2 1
A 8136 6 0 3519 1 2 1
A 8140 8 0 0 1 6440 1
A 8141 8 0 0 1 6441 1
A 8142 8 0 0 1 3516 0
T 8879 3416 0 3 0 1
A 8150 7 3512 0 1 2 1
A 8151 7 0 0 1 2 1
A 8149 6 0 3631 1 2 1
A 8157 7 3514 0 1 2 1
A 8158 7 0 0 1 2 1
A 8156 6 0 3519 1 2 1
A 8165 7 3516 0 1 2 1
A 8166 7 0 0 1 2 1
A 8164 6 0 3631 1 2 1
A 8172 7 3518 0 1 2 1
A 8173 7 0 0 1 2 1
A 8171 6 0 3519 1 2 1
A 8179 7 3520 0 1 2 1
A 8180 7 0 0 1 2 1
A 8178 6 0 3519 1 2 1
A 8186 7 3522 0 1 2 1
A 8187 7 0 0 1 2 1
A 8185 6 0 3519 1 2 1
A 8193 7 3524 0 1 2 1
A 8194 7 0 0 1 2 1
A 8192 6 0 3519 1 2 1
A 8200 7 3526 0 1 2 1
A 8201 7 0 0 1 2 1
A 8199 6 0 3519 1 2 1
A 8207 7 3528 0 1 2 1
A 8208 7 0 0 1 2 1
A 8206 6 0 3519 1 2 1
A 8214 7 3530 0 1 2 1
A 8215 7 0 0 1 2 1
A 8213 6 0 3519 1 2 1
A 8221 7 3532 0 1 2 1
A 8222 7 0 0 1 2 1
A 8220 6 0 3519 1 2 1
A 8228 7 3534 0 1 2 1
A 8229 7 0 0 1 2 1
A 8227 6 0 3519 1 2 1
A 8235 7 3536 0 1 2 1
A 8236 7 0 0 1 2 1
A 8234 6 0 3519 1 2 1
A 8242 7 3538 0 1 2 1
A 8243 7 0 0 1 2 1
A 8241 6 0 3519 1 2 1
A 8249 7 3540 0 1 2 1
A 8250 7 0 0 1 2 1
A 8248 6 0 3519 1 2 1
A 8252 6 0 0 1 2 1
A 8253 8 0 0 1 3516 1
A 8254 8 0 0 1 3516 1
A 8255 6 0 0 1 2 1
A 8256 16 0 0 1 3514 1
A 8257 16 0 0 1 3515 0
T 8880 3542 0 3 0 1
A 8264 7 3614 0 1 2 1
A 8265 7 0 0 1 2 1
A 8263 6 0 3631 1 2 1
A 8271 7 3616 0 1 2 1
A 8272 7 0 0 1 2 1
A 8270 6 0 3519 1 2 1
A 8274 6 0 0 1 3 1
T 8275 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
T 8276 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
A 8324 7 3618 0 1 2 1
A 8325 7 0 0 1 2 1
A 8323 6 0 69 1 2 1
A 8331 7 3620 0 1 2 1
A 8332 7 0 0 1 2 1
A 8330 6 0 3519 1 2 1
A 8336 8 0 0 1 6442 1
A 8337 8 0 0 1 6442 1
A 8338 6 0 0 1 183 1
A 8339 8 0 0 1 6443 1
A 8340 6 0 0 1 2 1
A 8341 9 0 0 1 3732 1
A 8342 9 0 0 1 3736 1
A 8343 9 0 0 1 6444 1
A 8344 8 0 0 1 3516 0
T 8881 3622 0 3 0 1
A 8346 2613 0 0 1 6439 1
A 8347 2613 0 0 1 6439 1
A 8350 7 3634 0 1 2 1
A 8354 7 3636 0 1 2 0
T 8882 3638 0 3 0 1
A 8357 8 0 0 1 3516 1
A 8358 8 0 0 1 6445 1
A 8359 8 0 0 1 3516 1
A 8360 8 0 0 1 6446 1
A 8361 8 0 0 1 6447 1
A 8362 8 0 0 1 6447 1
A 8363 9 0 0 1 6448 1
A 8364 9 0 0 1 6448 1
A 8365 8 0 0 1 6449 1
A 8366 9 0 0 1 6450 1
A 8367 9 0 0 1 6451 1
A 8368 9 0 0 1 6452 1
A 8369 9 0 0 1 6452 1
A 8370 8 0 0 1 3516 1
A 8371 8 0 0 1 3516 1
A 8372 8 0 0 1 3516 1
A 8373 6 0 0 1 3 1
A 8374 6 0 0 1 2 1
A 8378 7 3662 0 1 2 1
A 8379 7 0 0 1 2 1
A 8377 6 0 69 1 2 1
A 8384 7 3664 0 1 2 1
A 8385 7 0 0 1 2 1
A 8383 6 0 69 1 2 1
A 8390 7 3666 0 1 2 1
A 8391 7 0 0 1 2 1
A 8389 6 0 69 1 2 1
A 8393 6 0 0 1 3499 1
A 8394 6 0 0 1 6453 1
A 8395 6 0 0 1 3 1
A 8396 6 0 0 1 3 0
T 8883 3698 0 3 0 1
A 8423 7 3758 0 1 2 1
A 8424 7 0 0 1 2 1
A 8422 6 0 69 1 2 1
A 8429 7 3760 0 1 2 1
A 8430 7 0 0 1 2 1
A 8428 6 0 69 1 2 1
A 8435 7 3762 0 1 2 1
A 8436 7 0 0 1 2 1
A 8434 6 0 69 1 2 1
A 8441 7 3764 0 1 2 1
A 8442 7 0 0 1 2 1
A 8440 6 0 69 1 2 1
A 8447 7 3766 0 1 2 1
A 8448 7 0 0 1 2 1
A 8446 6 0 69 1 2 1
A 8454 7 3768 0 1 2 1
A 8455 7 0 0 1 2 1
A 8453 6 0 3519 1 2 1
A 8461 7 3770 0 1 2 1
A 8462 7 0 0 1 2 1
A 8460 6 0 3519 1 2 1
A 8467 7 3772 0 1 2 1
A 8468 7 0 0 1 2 1
A 8466 6 0 69 1 2 1
R 8471 3755 0 1
A 0 8 0 14 1 3516 0
A 8472 9 0 0 1 6454 1
A 8473 9 0 0 1 6455 1
A 8474 8 0 0 1 3516 1
A 8475 9 0 0 1 6456 1
A 8476 9 0 0 1 6457 1
A 8477 9 0 0 1 6455 1
A 8478 9 0 0 1 6456 1
A 8479 9 0 0 1 6456 1
A 8480 9 0 0 1 6456 0
T 8884 3774 0 3 0 1
T 8482 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
A 8486 7 3798 0 1 2 1
A 8487 7 0 0 1 2 1
A 8485 6 0 69 1 2 1
A 8492 7 3800 0 1 2 1
A 8493 7 0 0 1 2 1
A 8491 6 0 69 1 2 1
R 8495 3792 0 1
A 0 8 0 14 1 3516 0
R 8496 3795 0 1
A 0 8 0 41 1 3516 0
A 8497 6 0 0 1 2 0
T 8885 3802 0 3 0 1
A 8503 7 3856 0 1 2 1
A 8504 7 0 0 1 2 1
A 8502 6 0 3519 1 2 1
A 8510 7 3858 0 1 2 1
A 8511 7 0 0 1 2 1
A 8509 6 0 3519 1 2 1
A 8517 7 3860 0 1 2 1
A 8518 7 0 0 1 2 1
A 8516 6 0 3519 1 2 1
A 8525 7 3862 0 1 2 1
A 8526 7 0 0 1 2 1
A 8524 6 0 3631 1 2 1
A 8528 6 0 0 1 12 1
A 8529 8 0 0 1 3516 1
A 8533 7 3864 0 1 2 1
A 8534 7 0 0 1 2 1
A 8532 6 0 69 1 2 1
A 8539 7 3866 0 1 2 1
A 8540 7 0 0 1 2 1
A 8538 6 0 69 1 2 1
A 8545 7 3868 0 1 2 1
A 8546 7 0 0 1 2 1
A 8544 6 0 69 1 2 1
A 8551 7 3870 0 1 2 1
A 8552 7 0 0 1 2 1
A 8550 6 0 69 1 2 0
T 8886 3872 0 3 0 1
A 8560 7 4001 0 1 2 1
A 8561 7 0 0 1 2 1
A 8559 6 0 3631 1 2 1
A 8568 7 4003 0 1 2 1
A 8569 7 0 0 1 2 1
A 8567 6 0 3631 1 2 1
A 8576 7 4005 0 1 2 1
A 8577 7 0 0 1 2 1
A 8575 6 0 3631 1 2 1
A 8584 7 4007 0 1 2 1
A 8585 7 0 0 1 2 1
A 8583 6 0 3631 1 2 1
A 8590 7 4009 0 1 2 1
A 8591 7 0 0 1 2 1
A 8589 6 0 69 1 2 1
A 8596 7 4011 0 1 2 1
A 8597 7 0 0 1 2 1
A 8595 6 0 69 1 2 1
A 8602 7 4013 0 1 2 1
A 8603 7 0 0 1 2 1
A 8601 6 0 69 1 2 1
A 8608 7 4015 0 1 2 1
A 8609 7 0 0 1 2 1
A 8607 6 0 69 1 2 1
A 8615 7 4017 0 1 2 1
A 8616 7 0 0 1 2 1
A 8614 6 0 3519 1 2 1
A 8622 7 4019 0 1 2 1
A 8623 7 0 0 1 2 1
A 8621 6 0 3519 1 2 1
A 8629 7 4021 0 1 2 1
A 8630 7 0 0 1 2 1
A 8628 6 0 3519 1 2 1
A 8636 7 4023 0 1 2 1
A 8637 7 0 0 1 2 1
A 8635 6 0 3519 1 2 1
A 8643 7 4025 0 1 2 1
A 8644 7 0 0 1 2 1
A 8642 6 0 3519 1 2 1
A 8650 7 4027 0 1 2 1
A 8651 7 0 0 1 2 1
A 8649 6 0 3519 1 2 1
A 8657 7 4029 0 1 2 1
A 8658 7 0 0 1 2 1
A 8656 6 0 3519 1 2 1
A 8664 7 4031 0 1 2 1
A 8665 7 0 0 1 2 1
A 8663 6 0 3519 1 2 1
A 8672 7 4033 0 1 2 1
A 8673 7 0 0 1 2 1
A 8671 6 0 3631 1 2 1
A 8680 7 4035 0 1 2 1
A 8681 7 0 0 1 2 1
A 8679 6 0 3631 1 2 1
R 8683 3986 0 1
A 0 8 0 39 1 3516 0
R 8684 3989 0 1
A 0 8 0 39 1 3516 0
R 8685 3992 0 1
A 0 8 0 16 1 3516 0
A 8687 6 0 0 1 281 1
R 8688 3998 0 1
A 0 8 0 12 1 3516 0
A 8689 8 0 0 1 3516 1
A 8690 8 0 0 1 3516 1
A 8691 8 0 0 1 3516 1
A 8692 8 0 0 1 3516 1
A 8693 8 0 0 1 3516 1
A 8694 6 0 0 1 2 0
T 8887 4037 0 3 0 1
A 8700 7 4115 0 1 2 1
A 8701 7 0 0 1 2 1
A 8699 6 0 3519 1 2 1
A 8707 7 4117 0 1 2 1
A 8708 7 0 0 1 2 1
A 8706 6 0 3519 1 2 1
A 8714 7 4119 0 1 2 1
A 8715 7 0 0 1 2 1
A 8713 6 0 3519 1 2 1
A 8721 7 4121 0 1 2 1
A 8722 7 0 0 1 2 1
A 8720 6 0 3519 1 2 1
A 8728 7 4123 0 1 2 1
A 8729 7 0 0 1 2 1
A 8727 6 0 3519 1 2 1
A 8735 7 4125 0 1 2 1
A 8736 7 0 0 1 2 1
A 8734 6 0 3519 1 2 1
A 8742 7 4127 0 1 2 1
A 8743 7 0 0 1 2 1
A 8741 6 0 3519 1 2 1
A 8749 7 4129 0 1 2 1
A 8750 7 0 0 1 2 1
A 8748 6 0 3519 1 2 1
A 8756 7 4131 0 1 2 1
A 8757 7 0 0 1 2 1
A 8755 6 0 3519 1 2 1
A 8763 7 4133 0 1 2 1
A 8764 7 0 0 1 2 1
A 8762 6 0 3519 1 2 1
A 8770 7 4135 0 1 2 1
A 8771 7 0 0 1 2 1
A 8769 6 0 3519 1 2 1
A 8777 7 4137 0 1 2 1
A 8778 7 0 0 1 2 1
A 8776 6 0 3519 1 2 0
T 8888 4139 0 3 0 1
R 8781 4148 0 1
A 0 9 0 0 1 6458 1
A 0 9 0 0 1 6448 1
A 0 9 0 0 1 6456 1
A 0 9 0 0 1 6459 1
A 0 9 0 0 1 6457 0
A 8782 9 0 0 1 6456 1
A 8783 9 0 0 1 6456 1
A 8784 9 0 0 1 6456 1
A 8785 9 0 0 1 6460 1
A 8786 9 0 0 1 6454 1
A 8787 9 0 0 1 6457 1
A 8788 9 0 0 1 3736 1
A 8789 9 0 0 1 6459 1
A 8790 9 0 0 1 6461 0
T 8889 2679 0 3 0 1
A 7364 16 0 0 1 3515 1
A 7368 7 2709 0 1 2 1
A 7373 7 2711 0 1 2 1
A 7378 7 2713 0 1 2 1
A 7383 7 2715 0 1 2 0
T 8890 4151 0 3 0 1
A 8792 9 0 0 1 6462 1
A 8793 9 0 0 1 6463 1
A 8794 9 0 0 1 3736 1
A 8795 9 0 0 1 6464 1
A 8796 9 0 0 1 6465 1
A 8797 9 0 0 1 6466 1
A 8798 9 0 0 1 6454 1
A 8799 9 0 0 1 6467 1
A 8800 8 0 0 1 6468 1
A 8801 6 0 0 1 2 1
A 8802 6 0 0 1 39 1
A 8806 7 4187 0 1 2 1
A 8807 7 0 0 1 2 1
A 8805 6 0 69 1 2 1
A 8813 7 4189 0 1 2 1
A 8814 7 0 0 1 2 1
A 8812 6 0 3519 1 2 1
A 8820 7 4191 0 1 2 1
A 8821 7 0 0 1 2 1
A 8819 6 0 3519 1 2 1
A 8828 7 4193 0 1 2 1
A 8829 7 0 0 1 2 1
A 8827 6 0 3631 1 2 1
A 8836 7 4195 0 1 2 1
A 8837 7 0 0 1 2 1
A 8835 6 0 3631 1 2 0
T 8891 2667 0 3 0 1
A 6768 6 0 0 1 2 1
A 6770 6 0 0 1 2 0
T 8893 2655 0 3 0 1
A 6723 16 0 0 1 3515 1
A 6724 6 0 0 1 2 1
A 6725 6 0 0 1 2 1
A 6726 8 0 0 1 3575 1
A 6727 8 0 0 1 3576 1
A 6729 16 0 0 1 3515 1
T 6730 2649 0 3 0 1
A 6707 8 0 0 1 3574 0
A 6735 7 2661 0 1 2 1
A 6736 7 0 0 1 2 1
A 6734 6 0 3519 1 2 1
A 6742 7 2663 0 1 2 1
A 6743 7 0 0 1 2 1
A 6741 6 0 3519 1 2 1
A 6749 7 2665 0 1 2 1
A 6750 7 0 0 1 2 1
A 6748 6 0 3519 1 2 0
T 8894 4203 0 3 0 1
A 8853 7 4227 0 1 2 1
A 8854 7 0 0 1 2 1
A 8852 6 0 3519 1 2 1
A 8860 7 4229 0 1 2 1
A 8861 7 0 0 1 2 1
A 8859 6 0 3519 1 2 1
A 8867 7 4231 0 1 2 1
A 8868 7 0 0 1 2 1
A 8866 6 0 3519 1 2 0
T 8895 3668 0 3 0 1
A 8402 7 3692 0 1 2 1
A 8403 7 0 0 1 2 1
A 8401 6 0 3519 1 2 1
A 8409 7 3694 0 1 2 1
A 8410 7 0 0 1 2 1
A 8408 6 0 3519 1 2 1
A 8416 7 3696 0 1 2 1
A 8417 7 0 0 1 2 1
A 8415 6 0 3519 1 2 0
T 8896 2717 0 3 0 0
A 7074 6 0 0 1 2 1
A 7075 6 0 0 1 2 1
A 7076 6 0 0 1 3 0
T 8897 4239 0 3 0 0
T 8898 4233 0 3 0 1
T 8872 2723 0 3 0 1
A 7476 6 0 0 1 2 1
A 7477 6 0 0 1 2 1
A 7478 6 0 0 1 3 1
A 7484 7 2753 0 1 2 1
A 7485 7 0 0 1 2 1
A 7483 6 0 69 1 2 1
A 7490 7 2755 0 1 2 1
A 7491 7 0 0 1 2 1
A 7489 6 0 69 1 2 1
A 7496 7 2757 0 1 2 1
A 7497 7 0 0 1 2 1
A 7495 6 0 69 1 2 1
A 7502 7 2759 0 1 2 1
A 7503 7 0 0 1 2 1
A 7501 6 0 69 1 2 0
T 8873 2761 0 3 0 1
A 7546 6 0 0 1 3 1
A 7547 6 0 0 1 2 1
A 7548 6 0 0 1 2 1
A 7549 6 0 0 1 12 1
A 7550 6 0 0 1 3 1
A 7551 6 0 0 1 2 1
A 7552 6 0 0 1 2 1
A 7553 6 0 0 1 2 1
A 7554 6 0 0 1 2 1
A 7555 6 0 0 1 2 1
A 7556 6 0 0 1 2 1
A 7557 6 0 0 1 2 1
A 7558 6 0 0 1 14 1
A 7559 6 0 0 1 2 1
A 7560 6 0 0 1 36 1
A 7561 6 0 0 1 2 1
A 7562 6 0 0 1 2 1
A 7563 6 0 0 1 2 1
A 7564 6 0 0 1 2 1
A 7565 6 0 0 1 281 1
A 7566 6 0 0 1 2 1
A 7567 16 0 0 1 3514 1
A 7568 16 0 0 1 3515 1
A 7569 16 0 0 1 3515 1
A 7570 6 0 0 1 2 1
A 7571 6 0 0 1 2 1
A 7572 6 0 0 1 2 1
A 7573 6 0 0 1 2 1
A 7574 6 0 0 1 2 1
A 7575 6 0 0 1 2 1
A 7576 2613 0 0 1 6439 0
T 8874 2769 0 3 0 1
A 7582 7 2838 0 1 2 1
A 7583 7 0 0 1 2 1
A 7581 6 0 3519 1 2 1
A 7589 7 2840 0 1 2 1
A 7590 7 0 0 1 2 1
A 7588 6 0 3519 1 2 1
A 7596 7 2842 0 1 2 1
A 7597 7 0 0 1 2 1
A 7595 6 0 3519 1 2 1
A 7603 7 2844 0 1 2 1
A 7604 7 0 0 1 2 1
A 7602 6 0 3519 1 2 1
A 7610 7 2846 0 1 2 1
A 7611 7 0 0 1 2 1
A 7609 6 0 3519 1 2 1
A 7617 7 2848 0 1 2 1
A 7618 7 0 0 1 2 1
A 7616 6 0 3519 1 2 1
A 7625 7 2850 0 1 2 1
A 7626 7 0 0 1 2 1
A 7624 6 0 3631 1 2 1
A 7632 7 2852 0 1 2 1
A 7633 7 0 0 1 2 1
A 7631 6 0 3519 1 2 1
A 7639 7 2854 0 1 2 1
A 7640 7 0 0 1 2 1
A 7638 6 0 3519 1 2 1
A 7646 7 2856 0 1 2 1
A 7647 7 0 0 1 2 1
A 7645 6 0 3519 1 2 1
A 7649 6 0 0 1 2 1
R 7650 2835 0 1
A 0 6 0 14 1 2 0
A 7651 16 0 0 1 3514 0
T 8875 2858 0 3 0 1
A 7661 7 3020 0 1 2 1
A 7662 7 0 0 1 2 1
A 7660 6 0 3519 1 2 1
A 7668 7 3022 0 1 2 1
A 7669 7 0 0 1 2 1
A 7667 6 0 3519 1 2 1
A 7675 7 3024 0 1 2 1
A 7676 7 0 0 1 2 1
A 7674 6 0 3519 1 2 1
A 7682 7 3026 0 1 2 1
A 7683 7 0 0 1 2 1
A 7681 6 0 3519 1 2 1
A 7689 7 3028 0 1 2 1
A 7690 7 0 0 1 2 1
A 7688 6 0 3519 1 2 1
A 7696 7 3030 0 1 2 1
A 7697 7 0 0 1 2 1
A 7695 6 0 3519 1 2 1
A 7703 7 3032 0 1 2 1
A 7704 7 0 0 1 2 1
A 7702 6 0 3519 1 2 1
A 7710 7 3034 0 1 2 1
A 7711 7 0 0 1 2 1
A 7709 6 0 3519 1 2 1
A 7717 7 3036 0 1 2 1
A 7718 7 0 0 1 2 1
A 7716 6 0 3519 1 2 1
A 7724 7 3038 0 1 2 1
A 7725 7 0 0 1 2 1
A 7723 6 0 3519 1 2 1
A 7731 7 3040 0 1 2 1
A 7732 7 0 0 1 2 1
A 7730 6 0 3519 1 2 1
A 7738 7 3042 0 1 2 1
A 7739 7 0 0 1 2 1
A 7737 6 0 3519 1 2 1
A 7745 7 3044 0 1 2 1
A 7746 7 0 0 1 2 1
A 7744 6 0 3519 1 2 1
A 7752 7 3046 0 1 2 1
A 7753 7 0 0 1 2 1
A 7751 6 0 3519 1 2 1
A 7759 7 3048 0 1 2 1
A 7760 7 0 0 1 2 1
A 7758 6 0 3519 1 2 1
A 7766 7 3050 0 1 2 1
A 7767 7 0 0 1 2 1
A 7765 6 0 3519 1 2 1
A 7773 7 3052 0 1 2 1
A 7774 7 0 0 1 2 1
A 7772 6 0 3519 1 2 1
A 7780 7 3054 0 1 2 1
A 7781 7 0 0 1 2 1
A 7779 6 0 3519 1 2 1
A 7787 7 3056 0 1 2 1
A 7788 7 0 0 1 2 1
A 7786 6 0 3519 1 2 1
A 7794 7 3058 0 1 2 1
A 7795 7 0 0 1 2 1
A 7793 6 0 3519 1 2 1
A 7801 7 3060 0 1 2 1
A 7802 7 0 0 1 2 1
A 7800 6 0 3519 1 2 1
A 7808 7 3062 0 1 2 1
A 7809 7 0 0 1 2 1
A 7807 6 0 3519 1 2 1
A 7815 7 3064 0 1 2 1
A 7816 7 0 0 1 2 1
A 7814 6 0 3519 1 2 1
A 7822 7 3066 0 1 2 1
A 7823 7 0 0 1 2 1
A 7821 6 0 3519 1 2 1
A 7829 7 3068 0 1 2 1
A 7830 7 0 0 1 2 1
A 7828 6 0 3519 1 2 1
A 7836 7 3070 0 1 2 1
A 7837 7 0 0 1 2 1
A 7835 6 0 3519 1 2 0
T 8876 3126 0 3 0 1
A 7894 7 3270 0 1 2 1
A 7895 7 0 0 1 2 1
A 7893 6 0 3631 1 2 1
A 7902 7 3272 0 1 2 1
A 7903 7 0 0 1 2 1
A 7901 6 0 3631 1 2 1
A 7910 7 3274 0 1 2 1
A 7911 7 0 0 1 2 1
A 7909 6 0 3631 1 2 1
A 7918 7 3276 0 1 2 1
A 7919 7 0 0 1 2 1
A 7917 6 0 3631 1 2 1
A 7926 7 3278 0 1 2 1
A 7927 7 0 0 1 2 1
A 7925 6 0 3631 1 2 1
A 7934 7 3280 0 1 2 1
A 7935 7 0 0 1 2 1
A 7933 6 0 3631 1 2 1
A 7941 7 3282 0 1 2 1
A 7942 7 0 0 1 2 1
A 7940 6 0 3519 1 2 1
A 7948 7 3284 0 1 2 1
A 7949 7 0 0 1 2 1
A 7947 6 0 3519 1 2 1
A 7955 7 3286 0 1 2 1
A 7956 7 0 0 1 2 1
A 7954 6 0 3519 1 2 1
A 7962 7 3288 0 1 2 1
A 7963 7 0 0 1 2 1
A 7961 6 0 3519 1 2 1
A 7969 7 3290 0 1 2 1
A 7970 7 0 0 1 2 1
A 7968 6 0 3519 1 2 1
A 7976 7 3292 0 1 2 1
A 7977 7 0 0 1 2 1
A 7975 6 0 3519 1 2 1
A 7983 7 3294 0 1 2 1
A 7984 7 0 0 1 2 1
A 7982 6 0 3519 1 2 1
A 7990 7 3296 0 1 2 1
A 7991 7 0 0 1 2 1
A 7989 6 0 3519 1 2 1
A 7997 7 3298 0 1 2 1
A 7998 7 0 0 1 2 1
A 7996 6 0 3519 1 2 1
A 8004 7 3300 0 1 2 1
A 8005 7 0 0 1 2 1
A 8003 6 0 3519 1 2 1
A 8007 16 0 0 1 3515 1
A 8012 7 3302 0 1 2 1
A 8013 7 0 0 1 2 1
A 8011 6 0 3519 1 2 1
A 8019 7 3304 0 1 2 1
A 8020 7 0 0 1 2 1
A 8018 6 0 3519 1 2 1
A 8027 7 3306 0 1 2 1
A 8028 7 0 0 1 2 1
A 8026 6 0 3631 1 2 1
A 8034 7 3308 0 1 2 1
A 8035 7 0 0 1 2 1
A 8033 6 0 3519 1 2 1
A 8041 7 3310 0 1 2 1
A 8042 7 0 0 1 2 1
A 8040 6 0 3519 1 2 1
A 8048 7 3312 0 1 2 1
A 8049 7 0 0 1 2 1
A 8047 6 0 3519 1 2 1
A 8055 7 3314 0 1 2 1
A 8056 7 0 0 1 2 1
A 8054 6 0 3519 1 2 0
T 8877 3316 0 3 0 1
T 8059 3072 0 3 0 1
A 7845 7 3114 0 1 2 1
A 7846 7 0 0 1 2 1
A 7844 6 0 3631 1 2 1
A 7853 7 3116 0 1 2 1
A 7854 7 0 0 1 2 1
A 7852 6 0 3631 1 2 1
A 7861 7 3118 0 1 2 1
A 7862 7 0 0 1 2 1
A 7860 6 0 3631 1 2 1
A 7869 7 3120 0 1 2 1
A 7870 7 0 0 1 2 1
A 7868 6 0 3631 1 2 1
A 7877 7 3122 0 1 2 1
A 7878 7 0 0 1 2 1
A 7876 6 0 3631 1 2 1
A 7885 7 3124 0 1 2 1
A 7886 7 0 0 1 2 1
A 7884 6 0 3631 1 2 0
A 8065 7 3340 0 1 2 1
A 8066 7 0 0 1 2 1
A 8064 6 0 3631 1 2 1
A 8072 7 3342 0 1 2 1
A 8073 7 0 0 1 2 1
A 8071 6 0 3519 1 2 1
A 8079 7 3344 0 1 2 1
A 8080 7 0 0 1 2 1
A 8078 6 0 3519 1 2 0
T 8878 3346 0 3 0 1
A 8087 7 3400 0 1 2 1
A 8088 7 0 0 1 2 1
A 8086 6 0 3519 1 2 1
A 8094 7 3402 0 1 2 1
A 8095 7 0 0 1 2 1
A 8093 6 0 3519 1 2 1
A 8101 7 3404 0 1 2 1
A 8102 7 0 0 1 2 1
A 8100 6 0 3519 1 2 1
A 8108 7 3406 0 1 2 1
A 8109 7 0 0 1 2 1
A 8107 6 0 3519 1 2 1
A 8115 7 3408 0 1 2 1
A 8116 7 0 0 1 2 1
A 8114 6 0 3519 1 2 1
A 8122 7 3410 0 1 2 1
A 8123 7 0 0 1 2 1
A 8121 6 0 3519 1 2 1
A 8125 8 0 0 1 3516 1
A 8130 7 3412 0 1 2 1
A 8131 7 0 0 1 2 1
A 8129 6 0 3519 1 2 1
A 8137 7 3414 0 1 2 1
A 8138 7 0 0 1 2 1
A 8136 6 0 3519 1 2 1
A 8140 8 0 0 1 6440 1
A 8141 8 0 0 1 6441 1
A 8142 8 0 0 1 3516 0
T 8879 3416 0 3 0 1
A 8150 7 3512 0 1 2 1
A 8151 7 0 0 1 2 1
A 8149 6 0 3631 1 2 1
A 8157 7 3514 0 1 2 1
A 8158 7 0 0 1 2 1
A 8156 6 0 3519 1 2 1
A 8165 7 3516 0 1 2 1
A 8166 7 0 0 1 2 1
A 8164 6 0 3631 1 2 1
A 8172 7 3518 0 1 2 1
A 8173 7 0 0 1 2 1
A 8171 6 0 3519 1 2 1
A 8179 7 3520 0 1 2 1
A 8180 7 0 0 1 2 1
A 8178 6 0 3519 1 2 1
A 8186 7 3522 0 1 2 1
A 8187 7 0 0 1 2 1
A 8185 6 0 3519 1 2 1
A 8193 7 3524 0 1 2 1
A 8194 7 0 0 1 2 1
A 8192 6 0 3519 1 2 1
A 8200 7 3526 0 1 2 1
A 8201 7 0 0 1 2 1
A 8199 6 0 3519 1 2 1
A 8207 7 3528 0 1 2 1
A 8208 7 0 0 1 2 1
A 8206 6 0 3519 1 2 1
A 8214 7 3530 0 1 2 1
A 8215 7 0 0 1 2 1
A 8213 6 0 3519 1 2 1
A 8221 7 3532 0 1 2 1
A 8222 7 0 0 1 2 1
A 8220 6 0 3519 1 2 1
A 8228 7 3534 0 1 2 1
A 8229 7 0 0 1 2 1
A 8227 6 0 3519 1 2 1
A 8235 7 3536 0 1 2 1
A 8236 7 0 0 1 2 1
A 8234 6 0 3519 1 2 1
A 8242 7 3538 0 1 2 1
A 8243 7 0 0 1 2 1
A 8241 6 0 3519 1 2 1
A 8249 7 3540 0 1 2 1
A 8250 7 0 0 1 2 1
A 8248 6 0 3519 1 2 1
A 8252 6 0 0 1 2 1
A 8253 8 0 0 1 3516 1
A 8254 8 0 0 1 3516 1
A 8255 6 0 0 1 2 1
A 8256 16 0 0 1 3514 1
A 8257 16 0 0 1 3515 0
T 8880 3542 0 3 0 1
A 8264 7 3614 0 1 2 1
A 8265 7 0 0 1 2 1
A 8263 6 0 3631 1 2 1
A 8271 7 3616 0 1 2 1
A 8272 7 0 0 1 2 1
A 8270 6 0 3519 1 2 1
A 8274 6 0 0 1 3 1
T 8275 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
T 8276 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
A 8324 7 3618 0 1 2 1
A 8325 7 0 0 1 2 1
A 8323 6 0 69 1 2 1
A 8331 7 3620 0 1 2 1
A 8332 7 0 0 1 2 1
A 8330 6 0 3519 1 2 1
A 8336 8 0 0 1 6442 1
A 8337 8 0 0 1 6442 1
A 8338 6 0 0 1 183 1
A 8339 8 0 0 1 6443 1
A 8340 6 0 0 1 2 1
A 8341 9 0 0 1 3732 1
A 8342 9 0 0 1 3736 1
A 8343 9 0 0 1 6444 1
A 8344 8 0 0 1 3516 0
T 8881 3622 0 3 0 1
A 8346 2613 0 0 1 6439 1
A 8347 2613 0 0 1 6439 1
A 8350 7 3634 0 1 2 1
A 8354 7 3636 0 1 2 0
T 8882 3638 0 3 0 1
A 8357 8 0 0 1 3516 1
A 8358 8 0 0 1 6445 1
A 8359 8 0 0 1 3516 1
A 8360 8 0 0 1 6446 1
A 8361 8 0 0 1 6447 1
A 8362 8 0 0 1 6447 1
A 8363 9 0 0 1 6448 1
A 8364 9 0 0 1 6448 1
A 8365 8 0 0 1 6449 1
A 8366 9 0 0 1 6450 1
A 8367 9 0 0 1 6451 1
A 8368 9 0 0 1 6452 1
A 8369 9 0 0 1 6452 1
A 8370 8 0 0 1 3516 1
A 8371 8 0 0 1 3516 1
A 8372 8 0 0 1 3516 1
A 8373 6 0 0 1 3 1
A 8374 6 0 0 1 2 1
A 8378 7 3662 0 1 2 1
A 8379 7 0 0 1 2 1
A 8377 6 0 69 1 2 1
A 8384 7 3664 0 1 2 1
A 8385 7 0 0 1 2 1
A 8383 6 0 69 1 2 1
A 8390 7 3666 0 1 2 1
A 8391 7 0 0 1 2 1
A 8389 6 0 69 1 2 1
A 8393 6 0 0 1 3499 1
A 8394 6 0 0 1 6453 1
A 8395 6 0 0 1 3 1
A 8396 6 0 0 1 3 0
T 8883 3698 0 3 0 1
A 8423 7 3758 0 1 2 1
A 8424 7 0 0 1 2 1
A 8422 6 0 69 1 2 1
A 8429 7 3760 0 1 2 1
A 8430 7 0 0 1 2 1
A 8428 6 0 69 1 2 1
A 8435 7 3762 0 1 2 1
A 8436 7 0 0 1 2 1
A 8434 6 0 69 1 2 1
A 8441 7 3764 0 1 2 1
A 8442 7 0 0 1 2 1
A 8440 6 0 69 1 2 1
A 8447 7 3766 0 1 2 1
A 8448 7 0 0 1 2 1
A 8446 6 0 69 1 2 1
A 8454 7 3768 0 1 2 1
A 8455 7 0 0 1 2 1
A 8453 6 0 3519 1 2 1
A 8461 7 3770 0 1 2 1
A 8462 7 0 0 1 2 1
A 8460 6 0 3519 1 2 1
A 8467 7 3772 0 1 2 1
A 8468 7 0 0 1 2 1
A 8466 6 0 69 1 2 1
R 8471 3755 0 1
A 0 8 0 14 1 3516 0
A 8472 9 0 0 1 6454 1
A 8473 9 0 0 1 6455 1
A 8474 8 0 0 1 3516 1
A 8475 9 0 0 1 6456 1
A 8476 9 0 0 1 6457 1
A 8477 9 0 0 1 6455 1
A 8478 9 0 0 1 6456 1
A 8479 9 0 0 1 6456 1
A 8480 9 0 0 1 6456 0
T 8884 3774 0 3 0 1
T 8482 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
A 8486 7 3798 0 1 2 1
A 8487 7 0 0 1 2 1
A 8485 6 0 69 1 2 1
A 8492 7 3800 0 1 2 1
A 8493 7 0 0 1 2 1
A 8491 6 0 69 1 2 1
R 8495 3792 0 1
A 0 8 0 14 1 3516 0
R 8496 3795 0 1
A 0 8 0 41 1 3516 0
A 8497 6 0 0 1 2 0
T 8885 3802 0 3 0 1
A 8503 7 3856 0 1 2 1
A 8504 7 0 0 1 2 1
A 8502 6 0 3519 1 2 1
A 8510 7 3858 0 1 2 1
A 8511 7 0 0 1 2 1
A 8509 6 0 3519 1 2 1
A 8517 7 3860 0 1 2 1
A 8518 7 0 0 1 2 1
A 8516 6 0 3519 1 2 1
A 8525 7 3862 0 1 2 1
A 8526 7 0 0 1 2 1
A 8524 6 0 3631 1 2 1
A 8528 6 0 0 1 12 1
A 8529 8 0 0 1 3516 1
A 8533 7 3864 0 1 2 1
A 8534 7 0 0 1 2 1
A 8532 6 0 69 1 2 1
A 8539 7 3866 0 1 2 1
A 8540 7 0 0 1 2 1
A 8538 6 0 69 1 2 1
A 8545 7 3868 0 1 2 1
A 8546 7 0 0 1 2 1
A 8544 6 0 69 1 2 1
A 8551 7 3870 0 1 2 1
A 8552 7 0 0 1 2 1
A 8550 6 0 69 1 2 0
T 8886 3872 0 3 0 1
A 8560 7 4001 0 1 2 1
A 8561 7 0 0 1 2 1
A 8559 6 0 3631 1 2 1
A 8568 7 4003 0 1 2 1
A 8569 7 0 0 1 2 1
A 8567 6 0 3631 1 2 1
A 8576 7 4005 0 1 2 1
A 8577 7 0 0 1 2 1
A 8575 6 0 3631 1 2 1
A 8584 7 4007 0 1 2 1
A 8585 7 0 0 1 2 1
A 8583 6 0 3631 1 2 1
A 8590 7 4009 0 1 2 1
A 8591 7 0 0 1 2 1
A 8589 6 0 69 1 2 1
A 8596 7 4011 0 1 2 1
A 8597 7 0 0 1 2 1
A 8595 6 0 69 1 2 1
A 8602 7 4013 0 1 2 1
A 8603 7 0 0 1 2 1
A 8601 6 0 69 1 2 1
A 8608 7 4015 0 1 2 1
A 8609 7 0 0 1 2 1
A 8607 6 0 69 1 2 1
A 8615 7 4017 0 1 2 1
A 8616 7 0 0 1 2 1
A 8614 6 0 3519 1 2 1
A 8622 7 4019 0 1 2 1
A 8623 7 0 0 1 2 1
A 8621 6 0 3519 1 2 1
A 8629 7 4021 0 1 2 1
A 8630 7 0 0 1 2 1
A 8628 6 0 3519 1 2 1
A 8636 7 4023 0 1 2 1
A 8637 7 0 0 1 2 1
A 8635 6 0 3519 1 2 1
A 8643 7 4025 0 1 2 1
A 8644 7 0 0 1 2 1
A 8642 6 0 3519 1 2 1
A 8650 7 4027 0 1 2 1
A 8651 7 0 0 1 2 1
A 8649 6 0 3519 1 2 1
A 8657 7 4029 0 1 2 1
A 8658 7 0 0 1 2 1
A 8656 6 0 3519 1 2 1
A 8664 7 4031 0 1 2 1
A 8665 7 0 0 1 2 1
A 8663 6 0 3519 1 2 1
A 8672 7 4033 0 1 2 1
A 8673 7 0 0 1 2 1
A 8671 6 0 3631 1 2 1
A 8680 7 4035 0 1 2 1
A 8681 7 0 0 1 2 1
A 8679 6 0 3631 1 2 1
R 8683 3986 0 1
A 0 8 0 39 1 3516 0
R 8684 3989 0 1
A 0 8 0 39 1 3516 0
R 8685 3992 0 1
A 0 8 0 16 1 3516 0
A 8687 6 0 0 1 281 1
R 8688 3998 0 1
A 0 8 0 12 1 3516 0
A 8689 8 0 0 1 3516 1
A 8690 8 0 0 1 3516 1
A 8691 8 0 0 1 3516 1
A 8692 8 0 0 1 3516 1
A 8693 8 0 0 1 3516 1
A 8694 6 0 0 1 2 0
T 8887 4037 0 3 0 1
A 8700 7 4115 0 1 2 1
A 8701 7 0 0 1 2 1
A 8699 6 0 3519 1 2 1
A 8707 7 4117 0 1 2 1
A 8708 7 0 0 1 2 1
A 8706 6 0 3519 1 2 1
A 8714 7 4119 0 1 2 1
A 8715 7 0 0 1 2 1
A 8713 6 0 3519 1 2 1
A 8721 7 4121 0 1 2 1
A 8722 7 0 0 1 2 1
A 8720 6 0 3519 1 2 1
A 8728 7 4123 0 1 2 1
A 8729 7 0 0 1 2 1
A 8727 6 0 3519 1 2 1
A 8735 7 4125 0 1 2 1
A 8736 7 0 0 1 2 1
A 8734 6 0 3519 1 2 1
A 8742 7 4127 0 1 2 1
A 8743 7 0 0 1 2 1
A 8741 6 0 3519 1 2 1
A 8749 7 4129 0 1 2 1
A 8750 7 0 0 1 2 1
A 8748 6 0 3519 1 2 1
A 8756 7 4131 0 1 2 1
A 8757 7 0 0 1 2 1
A 8755 6 0 3519 1 2 1
A 8763 7 4133 0 1 2 1
A 8764 7 0 0 1 2 1
A 8762 6 0 3519 1 2 1
A 8770 7 4135 0 1 2 1
A 8771 7 0 0 1 2 1
A 8769 6 0 3519 1 2 1
A 8777 7 4137 0 1 2 1
A 8778 7 0 0 1 2 1
A 8776 6 0 3519 1 2 0
T 8888 4139 0 3 0 1
R 8781 4148 0 1
A 0 9 0 0 1 6458 1
A 0 9 0 0 1 6448 1
A 0 9 0 0 1 6456 1
A 0 9 0 0 1 6459 1
A 0 9 0 0 1 6457 0
A 8782 9 0 0 1 6456 1
A 8783 9 0 0 1 6456 1
A 8784 9 0 0 1 6456 1
A 8785 9 0 0 1 6460 1
A 8786 9 0 0 1 6454 1
A 8787 9 0 0 1 6457 1
A 8788 9 0 0 1 3736 1
A 8789 9 0 0 1 6459 1
A 8790 9 0 0 1 6461 0
T 8889 2679 0 3 0 1
A 7364 16 0 0 1 3515 1
A 7368 7 2709 0 1 2 1
A 7373 7 2711 0 1 2 1
A 7378 7 2713 0 1 2 1
A 7383 7 2715 0 1 2 0
T 8890 4151 0 3 0 1
A 8792 9 0 0 1 6462 1
A 8793 9 0 0 1 6463 1
A 8794 9 0 0 1 3736 1
A 8795 9 0 0 1 6464 1
A 8796 9 0 0 1 6465 1
A 8797 9 0 0 1 6466 1
A 8798 9 0 0 1 6454 1
A 8799 9 0 0 1 6467 1
A 8800 8 0 0 1 6468 1
A 8801 6 0 0 1 2 1
A 8802 6 0 0 1 39 1
A 8806 7 4187 0 1 2 1
A 8807 7 0 0 1 2 1
A 8805 6 0 69 1 2 1
A 8813 7 4189 0 1 2 1
A 8814 7 0 0 1 2 1
A 8812 6 0 3519 1 2 1
A 8820 7 4191 0 1 2 1
A 8821 7 0 0 1 2 1
A 8819 6 0 3519 1 2 1
A 8828 7 4193 0 1 2 1
A 8829 7 0 0 1 2 1
A 8827 6 0 3631 1 2 1
A 8836 7 4195 0 1 2 1
A 8837 7 0 0 1 2 1
A 8835 6 0 3631 1 2 0
T 8891 2667 0 3 0 1
A 6768 6 0 0 1 2 1
A 6770 6 0 0 1 2 0
T 8893 2655 0 3 0 1
A 6723 16 0 0 1 3515 1
A 6724 6 0 0 1 2 1
A 6725 6 0 0 1 2 1
A 6726 8 0 0 1 3575 1
A 6727 8 0 0 1 3576 1
A 6729 16 0 0 1 3515 1
T 6730 2649 0 3 0 1
A 6707 8 0 0 1 3574 0
A 6735 7 2661 0 1 2 1
A 6736 7 0 0 1 2 1
A 6734 6 0 3519 1 2 1
A 6742 7 2663 0 1 2 1
A 6743 7 0 0 1 2 1
A 6741 6 0 3519 1 2 1
A 6749 7 2665 0 1 2 1
A 6750 7 0 0 1 2 1
A 6748 6 0 3519 1 2 0
T 8894 4203 0 3 0 1
A 8853 7 4227 0 1 2 1
A 8854 7 0 0 1 2 1
A 8852 6 0 3519 1 2 1
A 8860 7 4229 0 1 2 1
A 8861 7 0 0 1 2 1
A 8859 6 0 3519 1 2 1
A 8867 7 4231 0 1 2 1
A 8868 7 0 0 1 2 1
A 8866 6 0 3519 1 2 0
T 8895 3668 0 3 0 1
A 8402 7 3692 0 1 2 1
A 8403 7 0 0 1 2 1
A 8401 6 0 3519 1 2 1
A 8409 7 3694 0 1 2 1
A 8410 7 0 0 1 2 1
A 8408 6 0 3519 1 2 1
A 8416 7 3696 0 1 2 1
A 8417 7 0 0 1 2 1
A 8415 6 0 3519 1 2 0
T 8896 2717 0 3 0 0
A 7074 6 0 0 1 2 1
A 7075 6 0 0 1 2 1
A 7076 6 0 0 1 3 0
T 8935 2525 0 3 0 1
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
T 8936 2525 0 3 0 0
A 681 7 2531 0 1 2 1
A 682 7 0 0 1 2 1
A 680 6 0 69 1 2 1
A 687 7 2533 0 1 2 1
A 688 7 0 0 1 2 1
A 686 6 0 69 1 2 1
A 693 7 2535 0 1 2 1
A 694 7 0 0 1 2 1
A 692 6 0 69 1 2 0
Z
