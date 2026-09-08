# results

Each `.txt` is the stdout of the script named by its prefix, produced on a GitHub Actions runner (or, under
`seam_gpu/`, on a Kaggle T4). `found/` holds the adversarial searcher's initial fields (`u`, 3 x N^3, with `Z0` and `T`).

- **adversarial** (3): `adversarial_ic.txt`, `adversarial_ic_z0tg_128.txt`, `adversarial_ic_z0tg_96.txt`
- **budgets** (3): `budgets3d_128.txt`, `budgets3d_64.txt`, `budgets3d_96.txt`
- **budgets.txt** (1): `budgets.txt`
- **burgers.txt** (1): `burgers.txt`
- **ckn** (3): `ckn_found_48.txt`, `ckn_found_64.txt`, `ckn_tg_64.txt`
- **criteria** (3): `criteria3d_abc_64.txt`, `criteria3d_kp_64.txt`, `criteria3d_tg_64.txt`
- **dyadic** (3): `dyadic_inviscid.txt`, `dyadic_nu1e-4.txt`, `dyadic_nu1e-6.txt`
- **frequency** (2): `frequency_matching_1d.txt`, `frequency_matching_2d.txt`
- **friction** (1): `friction_function_48.txt`
- **gpe** (3): `gpe_seam_2d.txt`, `gpe_seam_3d_64.txt`, `gpe_seam_3d_96.txt`
- **helicity** (5): `helicity_proj_0.0.txt`, `helicity_proj_0.25.txt`, `helicity_proj_0.5.txt`, `helicity_proj_0.75.txt`, `helicity_proj_0.9.txt`
- **help** (3): `help_fades_found_128.txt`, `help_fades_found_96.txt`, `help_fades_kp_96.txt`
- **jacobi** (11): `jacobi_kp_32.txt`, `jacobi_kp_48.txt`, `jacobi_kp_64.txt`, `jacobi_reverse_found_64.txt`, `jacobi_reverse_found_96.txt`, `jacobi_reverse_kp_64.txt`, `jacobi_reverse_kp_96.txt`, `jacobi_reverse_tg_64.txt`, `jacobi_tg_32.txt`, `jacobi_tg_48.txt`, `jacobi_tg_64.txt`
- **kolmogorov** (2): `kolmogorov_v2.txt`, `kolmogorov_v3.txt`
- **kolmogorov.txt** (1): `kolmogorov.txt`
- **kp** (4): `kp_128.txt`, `kp_64.txt`, `kp_64_nu1e-3.txt`, `kp_96.txt`
- **leashed** (4): `leashed64_dmin030_128.txt`, `leashed64_dmin030_192.txt`, `leashed_dmin030_128.txt`, `leashed_dmin045_128.txt`
- **liouville** (1): `liouville_taoclass.txt`
- **liouville.txt** (1): `liouville.txt`
- **lyapunov** (24): `lyapunov2d_64.txt`, `lyapunov2d_64_euler.txt`, `lyapunov2d_64_euler_v2.txt`, `lyapunov2d_64_v2.txt`, `lyapunov3d_v2_g1.txt`, `lyapunov3d_v2_g1_attack.txt`, `lyapunov3d_v3_relu.txt`, `lyapunov3d_v3_relu_attack.txt`, `lyapunov3d_v3_soft.txt`, `lyapunov3d_v3_soft_attack.txt`, `lyapunov_24_ci.txt`, `lyapunov_24_heads6_local.txt`, `lyapunov_24_local.txt`, `lyapunov_24_phess_local.txt`, `lyapunov_32_b1_attack.txt`, `lyapunov_32_ci.txt`, `lyapunov_32_p0_attack.txt`, `lyapunov_32_p1_attack.txt`, `lyapunov_32_r20_attack.txt`, `lyapunov_32_r20_p0.txt`, `lyapunov_32_r8_b1.txt`, `lyapunov_32_r8_g1.txt`, `lyapunov_32_r8_p0.txt`, `lyapunov_32_r8_p1.txt`
- **mandelbrot** (1): `mandelbrot_fraction_48.txt`
- **mechanism** (2): `mechanism3d_abc_64.txt`, `mechanism3d_tg_64.txt`
- **memory** (11): `memory3d_32_ppc2_nu0.txt`, `memory3d_32_ppc3_nu0.txt`, `memory3d_32_ppc3_nu2e-3.txt`, `memory3d_48_ppc2_nu0.txt`, `memory3d_48_ppc2_nu2e-3.txt`, `memory_paths_128_nu0.txt`, `memory_paths_128_nu1e-2.txt`, `memory_paths_128_ppc2.txt`, `memory_paths_128_ppc4.txt`, `memory_paths_128_ppc8.txt`, `memory_paths_256_ppc4.txt`
- **midpoint.txt** (1): `midpoint.txt`
- **minimal** (2): `minimal_datum_32.txt`, `minimal_datum_48.txt`
- **monotone** (1): `monotone_vs_dominating.txt`
- **nilpotent** (2): `nilpotent_sheet_48.txt`, `nilpotent_sheet_48b.txt`
- **ns** (3): `ns_d_1d.txt`, `ns_d_2d.txt`, `ns_d_3d.txt`
- **overnight** (5): `overnight_enstrophy_96.txt`, `overnight_helicity_0.txt`, `overnight_helicity_05.txt`, `overnight_helicity_09.txt`, `overnight_jacobi.txt`
- **pressure** (2): `pressure_share_48.txt`, `pressure_share_beta_64.txt`
- **projection** (3): `projection_price_found_48.txt`, `projection_price_kp_48.txt`, `projection_price_tg_48.txt`
- **quiet** (14): `quiet32_T0.5_pressure_0.3.txt`, `quiet32_T0.5_pressure_0.4.txt`, `quiet32_T0.5_pressure_0.5.txt`, `quiet32_T1.5_pressure_0.3.txt`, `quiet32_T1.5_pressure_0.4.txt`, `quiet32_T1.5_pressure_0.5.txt`, `quiet48_pressure_0.3.txt`, `quiet48_pressure_0.35.txt`, `quiet48_pressure_0.4.txt`, `quiet48_pressure_0.45.txt`, `quiet48_pressure_0.5.txt`, `quiet_pressure_0.3.txt`, `quiet_pressure_0.4.txt`, `quiet_pressure_0.5.txt`
- **seam** (4): `seam_race_64_nu0.txt`, `seam_race_64_nu1e-3.txt`, `seam_race_64_nu2e-3.txt`, `seam_race_96_nu2e-3.txt`
- **sheet** (1): `sheet_conjecture_64.txt`
- **spectra.txt** (1): `spectra.txt`
- **strip** (8): `strip_kp_128.txt`, `strip_kp_160.txt`, `strip_kp_192.txt`, `strip_kp_64.txt`, `strip_kp_96.txt`, `strip_tg_48.txt`, `strip_tg_64.txt`, `strip_tg_96.txt`
- **taylor** (1): `taylor_green.txt`
- **thinnest** (1): `thinnest_squeeze_48.txt`
- **twist** (4): `twist32_w10.txt`, `twist32_w3.txt`, `twist32_w30.txt`, `twist64_w10.txt`
- **type** (3): `type_found_128.txt`, `type_found_96.txt`, `type_kp_96.txt`
- **wave** (2): `wave_particle_1d.txt`, `wave_particle_2d.txt`

## seam_gpu (GPU seam race, by version)

- `seam_gpu/`: `seam_ICfound_NU0.txt`, `seam_ICfound_NU1e-2.txt`, `seam_ICfound_NU2e-3.txt`, `seam_ICpair_NU2e-3.txt`
- `seam_gpu\v2/`: `seam_ICfound_NU1e-3.txt`, `seam_ICfound_NU2e-3.txt`, `seam_ICfound_NU5e-4.txt`, `seam_ICpair_NU2e-3_T6.txt`
- `seam_gpu\v3/`: `seam_ICfound_NU1e-3_T2.0.txt`, `seam_ICfound_NU2e-3_T2.4.txt`, `seam_ICfound_NU5e-4_T1.6.txt`
- `seam_gpu\v4/`: `seam_ICfound_TAGckn_NU1e-3_T2.0.txt`, `seam_ICfound_TAGckn_NU2e-3_T2.4.txt`, `seam_ICkp_NU1e-3_T2.0.txt`, `seam_ICkp_NU2e-3_T2.4.txt`, `seam_ICpair_NU1e-3_T6.txt`, `seam_ICpair_NU2e-3_T6.txt`
- `seam_gpu\v5/`: `seam_ICfound_NU0_T1.0.txt`, `seam_ICpair_NU1e-3_T8.txt`, `seam_ICpair_NU2e-3_T8.txt`
- `seam_gpu\v5b/`: `seam_ICfound_NU0_T1.0.txt`, `seam_ICfound_NU0_T1.0_N128.txt`

- `found/`: `ckn64.npz`, `leashed64_dmin030.npz`, `minimal32.npz`, `minimal48.npz`, `quiet32_T0.5_0.3.npz`, `quiet32_T0.5_0.4.npz`, `quiet32_T0.5_0.5.npz`, `quiet32_T1.5_0.3.npz`, `quiet32_T1.5_0.4.npz`, `quiet32_T1.5_0.5.npz`, `quiet48_0.3.npz`, `quiet48_0.35.npz`, `quiet48_0.4.npz`, `quiet48_0.45.npz`, `quiet48_0.5.npz`, `quiet_0.3.npz`, `quiet_0.4.npz`, `quiet_0.5.npz`, `twist32_w10.npz`, `twist32_w3.npz`, `twist32_w30.npz`, `twist64_w10.npz`
