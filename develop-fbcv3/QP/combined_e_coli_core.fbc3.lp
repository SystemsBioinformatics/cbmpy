\ENCODING=ISO-8859-1
\Problem name: combined_e_coli_core.fbc3

Minimize
 R_BIOMASS_Ecoli_core_w_GAM_objective: [ R_BIOMASS_Ecoli_core_w_GAM_e_coli_1 ^2
                                 + R_BIOMASS_Ecoli_core_w_GAM_e_coli_2 ^2 ] / 2
Subject To
 M_13dpg_c_e_coli_1:   R_PGK_e_coli_1 + R_GAPD_e_coli_1  = 0
 M_13dpg_c_e_coli_2:   R_PGK_e_coli_2 + R_GAPD_e_coli_2  = 0
 M_2pg_c_e_coli_1:     - R_PGM_e_coli_1 - R_ENO_e_coli_1  = 0
 M_2pg_c_e_coli_2:     - R_PGM_e_coli_2 - R_ENO_e_coli_2  = 0
 M_3pg_c_e_coli_1:     - R_PGK_e_coli_1 + R_PGM_e_coli_1
                       - 1.496 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1  = 0
 M_3pg_c_e_coli_2:     - R_PGK_e_coli_2 + R_PGM_e_coli_2
                       - 1.496 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2  = 0
 M_6pgc_c_e_coli_1:    R_PGL_e_coli_1 - R_GND_e_coli_1  = 0
 M_6pgc_c_e_coli_2:    R_PGL_e_coli_2 - R_GND_e_coli_2  = 0
 M_6pgl_c_e_coli_1:    - R_PGL_e_coli_1 + R_G6PDH2r_e_coli_1  = 0
 M_6pgl_c_e_coli_2:    - R_PGL_e_coli_2 + R_G6PDH2r_e_coli_2  = 0
 M_ac_c_e_coli_1:      - R_ACKr_e_coli_1 + R_ACt2r_e_coli_1  = 0
 M_ac_c_e_coli_2:      - R_ACKr_e_coli_2 + R_ACt2r_e_coli_2  = 0
 M_ac_e:               - R_ACt2r_e_coli_1 - R_EX_ac_e - R_ACt2r_e_coli_2  = 0
 M_acald_c_e_coli_1:   - R_ACALD_e_coli_1 + R_ALCD2x_e_coli_1
                       + R_ACALDt_e_coli_1  = 0
 M_acald_c_e_coli_2:   - R_ACALD_e_coli_2 + R_ALCD2x_e_coli_2
                       + R_ACALDt_e_coli_2  = 0
 M_acald_e:            - R_ACALDt_e_coli_1 - R_EX_acald_e - R_ACALDt_e_coli_2
                        = 0
 M_accoa_c_e_coli_1:   R_PFL_e_coli_1 + R_ACALD_e_coli_1 - R_PTAr_e_coli_1
                       - 3.7478 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_CS_e_coli_1 - R_MALS_e_coli_1 + R_PDH_e_coli_1  = 0
 M_accoa_c_e_coli_2:   R_PFL_e_coli_2 + R_ACALD_e_coli_2 - R_PTAr_e_coli_2
                       - 3.7478 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_CS_e_coli_2 - R_MALS_e_coli_2 + R_PDH_e_coli_2  = 0
 M_acon_C_c_e_coli_1:  R_ACONTa_e_coli_1 - R_ACONTb_e_coli_1  = 0
 M_acon_C_c_e_coli_2:  R_ACONTa_e_coli_2 - R_ACONTb_e_coli_2  = 0
 M_actp_c_e_coli_1:    R_ACKr_e_coli_1 + R_PTAr_e_coli_1  = 0
 M_actp_c_e_coli_2:    R_ACKr_e_coli_2 + R_PTAr_e_coli_2  = 0
 M_adp_c_e_coli_1:     R_PFK_e_coli_1 + R_PGK_e_coli_1 + R_ACKr_e_coli_1
                       + R_ATPM_e_coli_1 + R_PPCK_e_coli_1 + 2 R_ADK1_e_coli_1
                       - R_ATPS4r_e_coli_1 - R_PYK_e_coli_1
                       + 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_SUCOAS_e_coli_1 + R_GLNS_e_coli_1
                       + R_GLNabc_e_coli_1  = 0
 M_adp_c_e_coli_2:     R_PFK_e_coli_2 + R_PGK_e_coli_2 + R_ACKr_e_coli_2
                       + R_ATPM_e_coli_2 + R_PPCK_e_coli_2 + 2 R_ADK1_e_coli_2
                       - R_ATPS4r_e_coli_2 - R_PYK_e_coli_2
                       + 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_SUCOAS_e_coli_2 + R_GLNS_e_coli_2
                       + R_GLNabc_e_coli_2  = 0
 M_akg_c_e_coli_1:     R_AKGt2r_e_coli_1 - R_AKGDH_e_coli_1
                       + 4.1182 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_GLUDy_e_coli_1 - R_GLUSy_e_coli_1
                       + R_ICDHyr_e_coli_1  = 0
 M_akg_c_e_coli_2:     R_AKGt2r_e_coli_2 - R_AKGDH_e_coli_2
                       + 4.1182 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_GLUDy_e_coli_2 - R_GLUSy_e_coli_2
                       + R_ICDHyr_e_coli_2  = 0
 M_akg_e:              - R_AKGt2r_e_coli_1 - R_EX_akg_e - R_AKGt2r_e_coli_2
                        = 0
 M_amp_c_e_coli_1:     R_PPS_e_coli_1 - R_ADK1_e_coli_1  = 0
 M_amp_c_e_coli_2:     R_PPS_e_coli_2 - R_ADK1_e_coli_2  = 0
 M_atp_c_e_coli_1:     - R_PFK_e_coli_1 - R_PGK_e_coli_1 - R_ACKr_e_coli_1
                       - R_ATPM_e_coli_1 - R_PPCK_e_coli_1 - R_PPS_e_coli_1
                       - R_ADK1_e_coli_1 + R_ATPS4r_e_coli_1 + R_PYK_e_coli_1
                       - 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_SUCOAS_e_coli_1 - R_GLNS_e_coli_1
                       - R_GLNabc_e_coli_1  = 0
 M_atp_c_e_coli_2:     - R_PFK_e_coli_2 - R_PGK_e_coli_2 - R_ACKr_e_coli_2
                       - R_ATPM_e_coli_2 - R_PPCK_e_coli_2 - R_PPS_e_coli_2
                       - R_ADK1_e_coli_2 + R_ATPS4r_e_coli_2 + R_PYK_e_coli_2
                       - 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_SUCOAS_e_coli_2 - R_GLNS_e_coli_2
                       - R_GLNabc_e_coli_2  = 0
 M_cit_c_e_coli_1:     - R_ACONTa_e_coli_1 + R_CS_e_coli_1  = 0
 M_cit_c_e_coli_2:     - R_ACONTa_e_coli_2 + R_CS_e_coli_2  = 0
 M_co2_c_e_coli_1:     - R_PPC_e_coli_1 + R_PPCK_e_coli_1 + R_AKGDH_e_coli_1
                       + R_CO2t_e_coli_1 + R_GND_e_coli_1 + R_ICDHyr_e_coli_1
                       + R_ME1_e_coli_1 + R_ME2_e_coli_1 + R_PDH_e_coli_1  = 0
 M_co2_c_e_coli_2:     - R_PPC_e_coli_2 + R_PPCK_e_coli_2 + R_AKGDH_e_coli_2
                       + R_CO2t_e_coli_2 + R_GND_e_coli_2 + R_ICDHyr_e_coli_2
                       + R_ME1_e_coli_2 + R_ME2_e_coli_2 + R_PDH_e_coli_2  = 0
 M_co2_e:              - R_CO2t_e_coli_1 - R_EX_co2_e - R_CO2t_e_coli_2  = 0
 M_coa_c_e_coli_1:     - R_PFL_e_coli_1 - R_ACALD_e_coli_1 - R_AKGDH_e_coli_1
                       + R_PTAr_e_coli_1
                       + 3.7478 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_CS_e_coli_1 - R_SUCOAS_e_coli_1 + R_MALS_e_coli_1
                       - R_PDH_e_coli_1  = 0
 M_coa_c_e_coli_2:     - R_PFL_e_coli_2 - R_ACALD_e_coli_2 - R_AKGDH_e_coli_2
                       + R_PTAr_e_coli_2
                       + 3.7478 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_CS_e_coli_2 - R_SUCOAS_e_coli_2 + R_MALS_e_coli_2
                       - R_PDH_e_coli_2  = 0
 M_dhap_c_e_coli_1:    - R_TPI_e_coli_1 + R_FBA_e_coli_1  = 0
 M_dhap_c_e_coli_2:    - R_TPI_e_coli_2 + R_FBA_e_coli_2  = 0
 M_e4p_c_e_coli_1:     - 0.361 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_TALA_e_coli_1 - R_TKT2_e_coli_1  = 0
 M_e4p_c_e_coli_2:     - 0.361 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_TALA_e_coli_2 - R_TKT2_e_coli_2  = 0
 M_etoh_c_e_coli_1:    - R_ALCD2x_e_coli_1 + R_ETOHt2r_e_coli_1  = 0
 M_etoh_c_e_coli_2:    - R_ALCD2x_e_coli_2 + R_ETOHt2r_e_coli_2  = 0
 M_etoh_e:             - R_ETOHt2r_e_coli_1 - R_EX_etoh_e - R_ETOHt2r_e_coli_2
                        = 0
 M_f6p_c_e_coli_1:     - R_PFK_e_coli_1 + R_PGI_e_coli_1
                       - 0.0709 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_TALA_e_coli_1 + R_TKT2_e_coli_1 + R_FBP_e_coli_1
                       + R_FRUpts2_e_coli_1  = 0
 M_f6p_c_e_coli_2:     - R_PFK_e_coli_2 + R_PGI_e_coli_2
                       - 0.0709 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_TALA_e_coli_2 + R_TKT2_e_coli_2 + R_FBP_e_coli_2
                       + R_FRUpts2_e_coli_2  = 0
 M_fdp_c_e_coli_1:     R_PFK_e_coli_1 - R_FBA_e_coli_1 - R_FBP_e_coli_1  = 0
 M_fdp_c_e_coli_2:     R_PFK_e_coli_2 - R_FBA_e_coli_2 - R_FBP_e_coli_2  = 0
 M_for_c_e_coli_1:     R_PFL_e_coli_1 + R_FORt2_e_coli_1 + R_FORt_e_coli_1  = 0
 M_for_c_e_coli_2:     R_PFL_e_coli_2 + R_FORt2_e_coli_2 + R_FORt_e_coli_2  = 0
 M_for_e:              - R_EX_for_e - R_FORt2_e_coli_1 - R_FORt_e_coli_1
                       - R_FORt2_e_coli_2 - R_FORt_e_coli_2  = 0
 M_fru_e:              - R_EX_fru_e - R_FRUpts2_e_coli_1 - R_FRUpts2_e_coli_2
                        = 0
 M_fum_c_e_coli_1:     R_SUCDi_e_coli_1 - R_FRD7_e_coli_1 - R_FUM_e_coli_1
                       + R_FUMt2_2_e_coli_1  = 0
 M_fum_c_e_coli_2:     R_SUCDi_e_coli_2 - R_FRD7_e_coli_2 - R_FUM_e_coli_2
                       + R_FUMt2_2_e_coli_2  = 0
 M_fum_e:              - R_EX_fum_e - R_FUMt2_2_e_coli_1 - R_FUMt2_2_e_coli_2
                        = 0
 M_g3p_c_e_coli_1:     - 0.129 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_TALA_e_coli_1 + R_TKT1_e_coli_1 + R_TKT2_e_coli_1
                       + R_TPI_e_coli_1 + R_FBA_e_coli_1 - R_GAPD_e_coli_1  = 0
 M_g3p_c_e_coli_2:     - 0.129 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_TALA_e_coli_2 + R_TKT1_e_coli_2 + R_TKT2_e_coli_2
                       + R_TPI_e_coli_2 + R_FBA_e_coli_2 - R_GAPD_e_coli_2  = 0
 M_g6p_c_e_coli_1:     - R_PGI_e_coli_1
                       - 0.205 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_G6PDH2r_e_coli_1 + R_GLCpts_e_coli_1  = 0
 M_g6p_c_e_coli_2:     - R_PGI_e_coli_2
                       - 0.205 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_G6PDH2r_e_coli_2 + R_GLCpts_e_coli_2  = 0
 M_glc__D_e:           - R_EX_glc__D_e - R_GLCpts_e_coli_1 - R_GLCpts_e_coli_2
                        = 0
 M_gln__L_c_e_coli_1:  - 0.2557 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_GLNS_e_coli_1 + R_GLNabc_e_coli_1 - R_GLUN_e_coli_1
                       - R_GLUSy_e_coli_1  = 0
 M_gln__L_c_e_coli_2:  - 0.2557 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_GLNS_e_coli_2 + R_GLNabc_e_coli_2 - R_GLUN_e_coli_2
                       - R_GLUSy_e_coli_2  = 0
 M_gln__L_e:           - R_EX_gln__L_e - R_GLNabc_e_coli_1 - R_GLNabc_e_coli_2
                        = 0
 M_glu__L_c_e_coli_1:  - 4.9414 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_GLNS_e_coli_1 - R_GLUDy_e_coli_1 + R_GLUN_e_coli_1
                       + 2 R_GLUSy_e_coli_1 + R_GLUt2r_e_coli_1  = 0
 M_glu__L_c_e_coli_2:  - 4.9414 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_GLNS_e_coli_2 - R_GLUDy_e_coli_2 + R_GLUN_e_coli_2
                       + 2 R_GLUSy_e_coli_2 + R_GLUt2r_e_coli_2  = 0
 M_glu__L_e:           - R_EX_glu__L_e - R_GLUt2r_e_coli_1 - R_GLUt2r_e_coli_2
                        = 0
 M_glx_c_e_coli_1:     R_ICL_e_coli_1 - R_MALS_e_coli_1  = 0
 M_glx_c_e_coli_2:     R_ICL_e_coli_2 - R_MALS_e_coli_2  = 0
 M_h2o_c_e_coli_1:     - R_PGL_e_coli_1 - R_PPC_e_coli_1 + R_ACONTa_e_coli_1
                       - R_ACONTb_e_coli_1 - R_ATPM_e_coli_1 - R_PPS_e_coli_1
                       + R_ATPS4r_e_coli_1
                       - 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_CS_e_coli_1 + R_CYTBD_e_coli_1 + R_ENO_e_coli_1
                       - R_FBP_e_coli_1 - R_FUM_e_coli_1 - R_GLNabc_e_coli_1
                       - R_GLUDy_e_coli_1 - R_GLUN_e_coli_1 + R_H2Ot_e_coli_1
                       - R_MALS_e_coli_1  = 0
 M_h2o_c_e_coli_2:     - R_PGL_e_coli_2 - R_PPC_e_coli_2 + R_ACONTa_e_coli_2
                       - R_ACONTb_e_coli_2 - R_ATPM_e_coli_2 - R_PPS_e_coli_2
                       + R_ATPS4r_e_coli_2
                       - 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_CS_e_coli_2 + R_CYTBD_e_coli_2 + R_ENO_e_coli_2
                       - R_FBP_e_coli_2 - R_FUM_e_coli_2 - R_GLNabc_e_coli_2
                       - R_GLUDy_e_coli_2 - R_GLUN_e_coli_2 + R_H2Ot_e_coli_2
                       - R_MALS_e_coli_2  = 0
 M_h2o_e:              - R_EX_h2o_e - R_H2Ot_e_coli_1 - R_H2Ot_e_coli_2  = 0
 M_h_c_e_coli_1:       R_PFK_e_coli_1 + R_PGL_e_coli_1 + R_ACALD_e_coli_1
                       + R_AKGt2r_e_coli_1 + R_PIt2r_e_coli_1
                       + R_ALCD2x_e_coli_1 + R_PPC_e_coli_1 + R_ATPM_e_coli_1
                       + R_ACt2r_e_coli_1 + 2 R_PPS_e_coli_1
                       + 3 R_ATPS4r_e_coli_1 - R_PYK_e_coli_1
                       + 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_PYRt2_e_coli_1 + R_CS_e_coli_1
                       + 2 R_SUCCt2_2_e_coli_1 - 2 R_CYTBD_e_coli_1
                       + R_D_LACt2_e_coli_1 + R_SUCCt3_e_coli_1
                       + R_ETOHt2r_e_coli_1 + 2 R_THD2_e_coli_1
                       + R_FORt2_e_coli_1 + 2 R_FUMt2_2_e_coli_1
                       + R_G6PDH2r_e_coli_1 + R_GAPD_e_coli_1 + R_GLNS_e_coli_1
                       + R_GLNabc_e_coli_1 + R_GLUDy_e_coli_1
                       - R_GLUSy_e_coli_1 + R_GLUt2r_e_coli_1
                       + R_LDH_D_e_coli_1 + R_MALS_e_coli_1
                       + 2 R_MALt2_2_e_coli_1 + R_MDH_e_coli_1
                       - 4 R_NADH16_e_coli_1  = 0
 M_h_c_e_coli_2:       R_PFK_e_coli_2 + R_PGL_e_coli_2 + R_ACALD_e_coli_2
                       + R_AKGt2r_e_coli_2 + R_PIt2r_e_coli_2
                       + R_ALCD2x_e_coli_2 + R_PPC_e_coli_2 + R_ATPM_e_coli_2
                       + R_ACt2r_e_coli_2 + 2 R_PPS_e_coli_2
                       + 3 R_ATPS4r_e_coli_2 - R_PYK_e_coli_2
                       + 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_PYRt2_e_coli_2 + R_CS_e_coli_2
                       + 2 R_SUCCt2_2_e_coli_2 - 2 R_CYTBD_e_coli_2
                       + R_D_LACt2_e_coli_2 + R_SUCCt3_e_coli_2
                       + R_ETOHt2r_e_coli_2 + 2 R_THD2_e_coli_2
                       + R_FORt2_e_coli_2 + 2 R_FUMt2_2_e_coli_2
                       + R_G6PDH2r_e_coli_2 + R_GAPD_e_coli_2 + R_GLNS_e_coli_2
                       + R_GLNabc_e_coli_2 + R_GLUDy_e_coli_2
                       - R_GLUSy_e_coli_2 + R_GLUt2r_e_coli_2
                       + R_LDH_D_e_coli_2 + R_MALS_e_coli_2
                       + 2 R_MALt2_2_e_coli_2 + R_MDH_e_coli_2
                       - 4 R_NADH16_e_coli_2  = 0
 M_h_e:                - R_AKGt2r_e_coli_1 - R_PIt2r_e_coli_1
                       - R_ACt2r_e_coli_1 - 4 R_ATPS4r_e_coli_1
                       - R_PYRt2_e_coli_1 - 2 R_SUCCt2_2_e_coli_1
                       + 2 R_CYTBD_e_coli_1 - R_D_LACt2_e_coli_1
                       - R_SUCCt3_e_coli_1 - R_ETOHt2r_e_coli_1
                       - 2 R_THD2_e_coli_1 - R_EX_h_e - R_FORt2_e_coli_1
                       - 2 R_FUMt2_2_e_coli_1 - R_GLUt2r_e_coli_1
                       - 2 R_MALt2_2_e_coli_1 + 3 R_NADH16_e_coli_1
                       - R_AKGt2r_e_coli_2 - R_PIt2r_e_coli_2
                       - R_ACt2r_e_coli_2 - 4 R_ATPS4r_e_coli_2
                       - R_PYRt2_e_coli_2 - 2 R_SUCCt2_2_e_coli_2
                       + 2 R_CYTBD_e_coli_2 - R_D_LACt2_e_coli_2
                       - R_SUCCt3_e_coli_2 - R_ETOHt2r_e_coli_2
                       - 2 R_THD2_e_coli_2 - R_FORt2_e_coli_2
                       - 2 R_FUMt2_2_e_coli_2 - R_GLUt2r_e_coli_2
                       - 2 R_MALt2_2_e_coli_2 + 3 R_NADH16_e_coli_2  = 0
 M_icit_c_e_coli_1:    R_ACONTb_e_coli_1 - R_ICDHyr_e_coli_1 - R_ICL_e_coli_1
                        = 0
 M_icit_c_e_coli_2:    R_ACONTb_e_coli_2 - R_ICDHyr_e_coli_2 - R_ICL_e_coli_2
                        = 0
 M_lac__D_c_e_coli_1:  R_D_LACt2_e_coli_1 - R_LDH_D_e_coli_1  = 0
 M_lac__D_c_e_coli_2:  R_D_LACt2_e_coli_2 - R_LDH_D_e_coli_2  = 0
 M_lac__D_e:           - R_D_LACt2_e_coli_1 - R_EX_lac__D_e
                       - R_D_LACt2_e_coli_2  = 0
 M_mal__L_c_e_coli_1:  R_FUM_e_coli_1 + R_MALS_e_coli_1 + R_MALt2_2_e_coli_1
                       - R_MDH_e_coli_1 - R_ME1_e_coli_1 - R_ME2_e_coli_1  = 0
 M_mal__L_c_e_coli_2:  R_FUM_e_coli_2 + R_MALS_e_coli_2 + R_MALt2_2_e_coli_2
                       - R_MDH_e_coli_2 - R_ME1_e_coli_2 - R_ME2_e_coli_2  = 0
 M_mal__L_e:           - R_EX_mal__L_e - R_MALt2_2_e_coli_1
                       - R_MALt2_2_e_coli_2  = 0
 M_nad_c_e_coli_1:     - R_ACALD_e_coli_1 - R_ALCD2x_e_coli_1
                       - R_AKGDH_e_coli_1
                       - 3.547 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_THD2_e_coli_1 - R_GAPD_e_coli_1 - R_LDH_D_e_coli_1
                       - R_MDH_e_coli_1 - R_ME1_e_coli_1 + R_NADH16_e_coli_1
                       - R_NADTRHD_e_coli_1 - R_PDH_e_coli_1  = 0
 M_nad_c_e_coli_2:     - R_ACALD_e_coli_2 - R_ALCD2x_e_coli_2
                       - R_AKGDH_e_coli_2
                       - 3.547 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_THD2_e_coli_2 - R_GAPD_e_coli_2 - R_LDH_D_e_coli_2
                       - R_MDH_e_coli_2 - R_ME1_e_coli_2 + R_NADH16_e_coli_2
                       - R_NADTRHD_e_coli_2 - R_PDH_e_coli_2  = 0
 M_nadh_c_e_coli_1:    R_ACALD_e_coli_1 + R_ALCD2x_e_coli_1 + R_AKGDH_e_coli_1
                       + 3.547 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_THD2_e_coli_1 + R_GAPD_e_coli_1 + R_LDH_D_e_coli_1
                       + R_MDH_e_coli_1 + R_ME1_e_coli_1 - R_NADH16_e_coli_1
                       + R_NADTRHD_e_coli_1 + R_PDH_e_coli_1  = 0
 M_nadh_c_e_coli_2:    R_ACALD_e_coli_2 + R_ALCD2x_e_coli_2 + R_AKGDH_e_coli_2
                       + 3.547 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_THD2_e_coli_2 + R_GAPD_e_coli_2 + R_LDH_D_e_coli_2
                       + R_MDH_e_coli_2 + R_ME1_e_coli_2 - R_NADH16_e_coli_2
                       + R_NADTRHD_e_coli_2 + R_PDH_e_coli_2  = 0
 M_nadp_c_e_coli_1:    13.0279 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_THD2_e_coli_1 - R_G6PDH2r_e_coli_1
                       - R_GLUDy_e_coli_1 + R_GLUSy_e_coli_1 - R_GND_e_coli_1
                       - R_ICDHyr_e_coli_1 - R_ME2_e_coli_1
                       + R_NADTRHD_e_coli_1  = 0
 M_nadp_c_e_coli_2:    13.0279 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_THD2_e_coli_2 - R_G6PDH2r_e_coli_2
                       - R_GLUDy_e_coli_2 + R_GLUSy_e_coli_2 - R_GND_e_coli_2
                       - R_ICDHyr_e_coli_2 - R_ME2_e_coli_2
                       + R_NADTRHD_e_coli_2  = 0
 M_nadph_c_e_coli_1:   - 13.0279 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_THD2_e_coli_1 + R_G6PDH2r_e_coli_1
                       + R_GLUDy_e_coli_1 - R_GLUSy_e_coli_1 + R_GND_e_coli_1
                       + R_ICDHyr_e_coli_1 + R_ME2_e_coli_1
                       - R_NADTRHD_e_coli_1  = 0
 M_nadph_c_e_coli_2:   - 13.0279 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_THD2_e_coli_2 + R_G6PDH2r_e_coli_2
                       + R_GLUDy_e_coli_2 - R_GLUSy_e_coli_2 + R_GND_e_coli_2
                       + R_ICDHyr_e_coli_2 + R_ME2_e_coli_2
                       - R_NADTRHD_e_coli_2  = 0
 M_nh4_c_e_coli_1:     - R_GLNS_e_coli_1 + R_GLUDy_e_coli_1 + R_GLUN_e_coli_1
                       + R_NH4t_e_coli_1  = 0
 M_nh4_c_e_coli_2:     - R_GLNS_e_coli_2 + R_GLUDy_e_coli_2 + R_GLUN_e_coli_2
                       + R_NH4t_e_coli_2  = 0
 M_nh4_e:              - R_EX_nh4_e - R_NH4t_e_coli_1 - R_NH4t_e_coli_2  = 0
 M_o2_c_e_coli_1:      - 0.5 R_CYTBD_e_coli_1 + R_O2t_e_coli_1  = 0
 M_o2_c_e_coli_2:      - 0.5 R_CYTBD_e_coli_2 + R_O2t_e_coli_2  = 0
 M_o2_e:               - R_EX_o2_e - R_O2t_e_coli_1 - R_O2t_e_coli_2  = 0
 M_oaa_c_e_coli_1:     R_PPC_e_coli_1 - R_PPCK_e_coli_1
                       - 1.7867 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_CS_e_coli_1 + R_MDH_e_coli_1  = 0
 M_oaa_c_e_coli_2:     R_PPC_e_coli_2 - R_PPCK_e_coli_2
                       - 1.7867 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_CS_e_coli_2 + R_MDH_e_coli_2  = 0
 M_pep_c_e_coli_1:     - R_PPC_e_coli_1 + R_PPCK_e_coli_1 + R_PPS_e_coli_1
                       - R_PYK_e_coli_1
                       - 0.5191 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_ENO_e_coli_1 - R_FRUpts2_e_coli_1
                       - R_GLCpts_e_coli_1  = 0
 M_pep_c_e_coli_2:     - R_PPC_e_coli_2 + R_PPCK_e_coli_2 + R_PPS_e_coli_2
                       - R_PYK_e_coli_2
                       - 0.5191 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_ENO_e_coli_2 - R_FRUpts2_e_coli_2
                       - R_GLCpts_e_coli_2  = 0
 M_pi_c_e_coli_1:      R_PIt2r_e_coli_1 + R_PPC_e_coli_1 + R_ATPM_e_coli_1
                       + R_PPS_e_coli_1 - R_ATPS4r_e_coli_1 - R_PTAr_e_coli_1
                       + 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_SUCOAS_e_coli_1 + R_FBP_e_coli_1 - R_GAPD_e_coli_1
                       + R_GLNS_e_coli_1 + R_GLNabc_e_coli_1  = 0
 M_pi_c_e_coli_2:      R_PIt2r_e_coli_2 + R_PPC_e_coli_2 + R_ATPM_e_coli_2
                       + R_PPS_e_coli_2 - R_ATPS4r_e_coli_2 - R_PTAr_e_coli_2
                       + 59.81 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_SUCOAS_e_coli_2 + R_FBP_e_coli_2 - R_GAPD_e_coli_2
                       + R_GLNS_e_coli_2 + R_GLNabc_e_coli_2  = 0
 M_pi_e:               - R_PIt2r_e_coli_1 - R_EX_pi_e - R_PIt2r_e_coli_2  = 0
 M_pyr_c_e_coli_1:     - R_PFL_e_coli_1 - R_PPS_e_coli_1 + R_PYK_e_coli_1
                       - 2.8328 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       + R_PYRt2_e_coli_1 + R_FRUpts2_e_coli_1
                       + R_GLCpts_e_coli_1 + R_LDH_D_e_coli_1 + R_ME1_e_coli_1
                       + R_ME2_e_coli_1 - R_PDH_e_coli_1  = 0
 M_pyr_c_e_coli_2:     - R_PFL_e_coli_2 - R_PPS_e_coli_2 + R_PYK_e_coli_2
                       - 2.8328 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       + R_PYRt2_e_coli_2 + R_FRUpts2_e_coli_2
                       + R_GLCpts_e_coli_2 + R_LDH_D_e_coli_2 + R_ME1_e_coli_2
                       + R_ME2_e_coli_2 - R_PDH_e_coli_2  = 0
 M_pyr_e:              - R_PYRt2_e_coli_1 - R_EX_pyr_e - R_PYRt2_e_coli_2  = 0
 M_q8_c_e_coli_1:      R_CYTBD_e_coli_1 - R_SUCDi_e_coli_1 + R_FRD7_e_coli_1
                       - R_NADH16_e_coli_1  = 0
 M_q8_c_e_coli_2:      R_CYTBD_e_coli_2 - R_SUCDi_e_coli_2 + R_FRD7_e_coli_2
                       - R_NADH16_e_coli_2  = 0
 M_q8h2_c_e_coli_1:    - R_CYTBD_e_coli_1 + R_SUCDi_e_coli_1 - R_FRD7_e_coli_1
                       + R_NADH16_e_coli_1  = 0
 M_q8h2_c_e_coli_2:    - R_CYTBD_e_coli_2 + R_SUCDi_e_coli_2 - R_FRD7_e_coli_2
                       + R_NADH16_e_coli_2  = 0
 M_r5p_c_e_coli_1:     - 0.8977 R_BIOMASS_Ecoli_core_w_GAM_e_coli_1
                       - R_RPI_e_coli_1 - R_TKT1_e_coli_1  = 0
 M_r5p_c_e_coli_2:     - 0.8977 R_BIOMASS_Ecoli_core_w_GAM_e_coli_2
                       - R_RPI_e_coli_2 - R_TKT1_e_coli_2  = 0
 M_ru5p__D_c_e_coli_1: - R_RPE_e_coli_1 + R_RPI_e_coli_1 + R_GND_e_coli_1  = 0
 M_ru5p__D_c_e_coli_2: - R_RPE_e_coli_2 + R_RPI_e_coli_2 + R_GND_e_coli_2  = 0
 M_s7p_c_e_coli_1:     - R_TALA_e_coli_1 + R_TKT1_e_coli_1  = 0
 M_s7p_c_e_coli_2:     - R_TALA_e_coli_2 + R_TKT1_e_coli_2  = 0
 M_succ_c_e_coli_1:    R_SUCCt2_2_e_coli_1 - R_SUCCt3_e_coli_1
                       - R_SUCDi_e_coli_1 - R_SUCOAS_e_coli_1 + R_FRD7_e_coli_1
                       + R_ICL_e_coli_1  = 0
 M_succ_c_e_coli_2:    R_SUCCt2_2_e_coli_2 - R_SUCCt3_e_coli_2
                       - R_SUCDi_e_coli_2 - R_SUCOAS_e_coli_2 + R_FRD7_e_coli_2
                       + R_ICL_e_coli_2  = 0
 M_succ_e:             - R_SUCCt2_2_e_coli_1 + R_SUCCt3_e_coli_1 - R_EX_succ_e
                       - R_SUCCt2_2_e_coli_2 + R_SUCCt3_e_coli_2  = 0
 M_succoa_c_e_coli_1:  R_AKGDH_e_coli_1 + R_SUCOAS_e_coli_1  = 0
 M_succoa_c_e_coli_2:  R_AKGDH_e_coli_2 + R_SUCOAS_e_coli_2  = 0
 M_xu5p__D_c_e_coli_1: R_RPE_e_coli_1 - R_TKT1_e_coli_1 - R_TKT2_e_coli_1  = 0
 M_xu5p__D_c_e_coli_2: R_RPE_e_coli_2 - R_TKT1_e_coli_2 - R_TKT2_e_coli_2  = 0
Bounds
 0 <= R_PFK_e_coli_1 <= 1000
 0 <= R_PFL_e_coli_1 <= 1000
-1000 <= R_PGI_e_coli_1 <= 1000
-1000 <= R_PGK_e_coli_1 <= 1000
 0 <= R_PGL_e_coli_1 <= 1000
-1000 <= R_ACALD_e_coli_1 <= 1000
-1000 <= R_AKGt2r_e_coli_1 <= 1000
-1000 <= R_PGM_e_coli_1 <= 1000
-1000 <= R_PIt2r_e_coli_1 <= 1000
-1000 <= R_ALCD2x_e_coli_1 <= 1000
-1000 <= R_ACALDt_e_coli_1 <= 1000
-1000 <= R_ACKr_e_coli_1 <= 1000
 0 <= R_PPC_e_coli_1 <= 1000
-1000 <= R_ACONTa_e_coli_1 <= 1000
-1000 <= R_ACONTb_e_coli_1 <= 1000
 8.39 <= R_ATPM_e_coli_1 <= 1000
 0 <= R_PPCK_e_coli_1 <= 1000
-1000 <= R_ACt2r_e_coli_1 <= 1000
 0 <= R_PPS_e_coli_1 <= 1000
-1000 <= R_ADK1_e_coli_1 <= 1000
 0 <= R_AKGDH_e_coli_1 <= 1000
-1000 <= R_ATPS4r_e_coli_1 <= 1000
-1000 <= R_PTAr_e_coli_1 <= 1000
 0 <= R_PYK_e_coli_1 <= 1000
 0 <= R_BIOMASS_Ecoli_core_w_GAM_e_coli_1 <= 1000
-1000 <= R_PYRt2_e_coli_1 <= 1000
-1000 <= R_CO2t_e_coli_1 <= 1000
-1000 <= R_RPE_e_coli_1 <= 1000
 0 <= R_CS_e_coli_1 <= 1000
-1000 <= R_RPI_e_coli_1 <= 1000
 0 <= R_SUCCt2_2_e_coli_1 <= 1000
 0 <= R_CYTBD_e_coli_1 <= 1000
-1000 <= R_D_LACt2_e_coli_1 <= 1000
-1000 <= R_ENO_e_coli_1 <= 1000
 0 <= R_SUCCt3_e_coli_1 <= 1000
-1000 <= R_ETOHt2r_e_coli_1 <= 1000
 0 <= R_SUCDi_e_coli_1 <= 1000
-1000 <= R_SUCOAS_e_coli_1 <= 1000
-1000 <= R_TALA_e_coli_1 <= 1000
 0 <= R_THD2_e_coli_1 <= 1000
-1000 <= R_TKT1_e_coli_1 <= 1000
-1000 <= R_TKT2_e_coli_1 <= 1000
-1000 <= R_TPI_e_coli_1 <= 1000
 0 <= R_EX_ac_e <= 1000
 0 <= R_EX_acald_e <= 1000
 0 <= R_EX_akg_e <= 1000
-1000 <= R_EX_co2_e <= 1000
 0 <= R_EX_etoh_e <= 1000
 0 <= R_EX_for_e <= 1000
 0 <= R_EX_fru_e <= 1000
 0 <= R_EX_fum_e <= 1000
-10 <= R_EX_glc__D_e <= 1000
 0 <= R_EX_gln__L_e <= 1000
 0 <= R_EX_glu__L_e <= 1000
-1000 <= R_EX_h_e <= 1000
-1000 <= R_EX_h2o_e <= 1000
 0 <= R_EX_lac__D_e <= 1000
 0 <= R_EX_mal__L_e <= 1000
-1000 <= R_EX_nh4_e <= 1000
-1000 <= R_EX_o2_e <= 1000
-1000 <= R_EX_pi_e <= 1000
 0 <= R_EX_pyr_e <= 1000
 0 <= R_EX_succ_e <= 1000
-1000 <= R_FBA_e_coli_1 <= 1000
 0 <= R_FBP_e_coli_1 <= 1000
 0 <= R_FORt2_e_coli_1 <= 1000
-1000 <= R_FORt_e_coli_1 <= 0
 0 <= R_FRD7_e_coli_1 <= 1000
 0 <= R_FRUpts2_e_coli_1 <= 1000
-1000 <= R_FUM_e_coli_1 <= 1000
 0 <= R_FUMt2_2_e_coli_1 <= 1000
-1000 <= R_G6PDH2r_e_coli_1 <= 1000
-1000 <= R_GAPD_e_coli_1 <= 1000
 0 <= R_GLCpts_e_coli_1 <= 1000
 0 <= R_GLNS_e_coli_1 <= 1000
 0 <= R_GLNabc_e_coli_1 <= 1000
-1000 <= R_GLUDy_e_coli_1 <= 1000
 0 <= R_GLUN_e_coli_1 <= 1000
 0 <= R_GLUSy_e_coli_1 <= 1000
-1000 <= R_GLUt2r_e_coli_1 <= 1000
 0 <= R_GND_e_coli_1 <= 1000
-1000 <= R_H2Ot_e_coli_1 <= 1000
-1000 <= R_ICDHyr_e_coli_1 <= 1000
 0 <= R_ICL_e_coli_1 <= 1000
-1000 <= R_LDH_D_e_coli_1 <= 1000
 0 <= R_MALS_e_coli_1 <= 1000
 0 <= R_MALt2_2_e_coli_1 <= 1000
-1000 <= R_MDH_e_coli_1 <= 1000
 0 <= R_ME1_e_coli_1 <= 1000
 0 <= R_ME2_e_coli_1 <= 1000
 0 <= R_NADH16_e_coli_1 <= 1000
 0 <= R_NADTRHD_e_coli_1 <= 1000
-1000 <= R_NH4t_e_coli_1 <= 1000
-1000 <= R_O2t_e_coli_1 <= 1000
 0 <= R_PDH_e_coli_1 <= 1000
 0 <= R_PFK_e_coli_2 <= 1000
 0 <= R_PFL_e_coli_2 <= 1000
-1000 <= R_PGI_e_coli_2 <= 1000
-1000 <= R_PGK_e_coli_2 <= 1000
 0 <= R_PGL_e_coli_2 <= 1000
-1000 <= R_ACALD_e_coli_2 <= 1000
-1000 <= R_AKGt2r_e_coli_2 <= 1000
-1000 <= R_PGM_e_coli_2 <= 1000
-1000 <= R_PIt2r_e_coli_2 <= 1000
-1000 <= R_ALCD2x_e_coli_2 <= 1000
-1000 <= R_ACALDt_e_coli_2 <= 1000
-1000 <= R_ACKr_e_coli_2 <= 1000
 0 <= R_PPC_e_coli_2 <= 1000
-1000 <= R_ACONTa_e_coli_2 <= 1000
-1000 <= R_ACONTb_e_coli_2 <= 1000
 8.39 <= R_ATPM_e_coli_2 <= 1000
 0 <= R_PPCK_e_coli_2 <= 1000
-1000 <= R_ACt2r_e_coli_2 <= 1000
 0 <= R_PPS_e_coli_2 <= 1000
-1000 <= R_ADK1_e_coli_2 <= 1000
 0 <= R_AKGDH_e_coli_2 <= 1000
-1000 <= R_ATPS4r_e_coli_2 <= 1000
-1000 <= R_PTAr_e_coli_2 <= 1000
 0 <= R_PYK_e_coli_2 <= 1000
 0 <= R_BIOMASS_Ecoli_core_w_GAM_e_coli_2 <= 1000
-1000 <= R_PYRt2_e_coli_2 <= 1000
-1000 <= R_CO2t_e_coli_2 <= 1000
-1000 <= R_RPE_e_coli_2 <= 1000
 0 <= R_CS_e_coli_2 <= 1000
-1000 <= R_RPI_e_coli_2 <= 1000
 0 <= R_SUCCt2_2_e_coli_2 <= 1000
 0 <= R_CYTBD_e_coli_2 <= 1000
-1000 <= R_D_LACt2_e_coli_2 <= 1000
-1000 <= R_ENO_e_coli_2 <= 1000
 0 <= R_SUCCt3_e_coli_2 <= 1000
-1000 <= R_ETOHt2r_e_coli_2 <= 1000
 0 <= R_SUCDi_e_coli_2 <= 1000
-1000 <= R_SUCOAS_e_coli_2 <= 1000
-1000 <= R_TALA_e_coli_2 <= 1000
 0 <= R_THD2_e_coli_2 <= 1000
-1000 <= R_TKT1_e_coli_2 <= 1000
-1000 <= R_TKT2_e_coli_2 <= 1000
-1000 <= R_TPI_e_coli_2 <= 1000
-1000 <= R_FBA_e_coli_2 <= 1000
 0 <= R_FBP_e_coli_2 <= 1000
 0 <= R_FORt2_e_coli_2 <= 1000
-1000 <= R_FORt_e_coli_2 <= 0
 0 <= R_FRD7_e_coli_2 <= 1000
 0 <= R_FRUpts2_e_coli_2 <= 1000
-1000 <= R_FUM_e_coli_2 <= 1000
 0 <= R_FUMt2_2_e_coli_2 <= 1000
-1000 <= R_G6PDH2r_e_coli_2 <= 1000
-1000 <= R_GAPD_e_coli_2 <= 1000
 10 <= R_GLCpts_e_coli_2 <= 1000
 0 <= R_GLNS_e_coli_2 <= 1000
 0 <= R_GLNabc_e_coli_2 <= 1000
-1000 <= R_GLUDy_e_coli_2 <= 1000
 0 <= R_GLUN_e_coli_2 <= 1000
 0 <= R_GLUSy_e_coli_2 <= 1000
-1000 <= R_GLUt2r_e_coli_2 <= 1000
 0 <= R_GND_e_coli_2 <= 1000
-1000 <= R_H2Ot_e_coli_2 <= 1000
-1000 <= R_ICDHyr_e_coli_2 <= 1000
 0 <= R_ICL_e_coli_2 <= 1000
-1000 <= R_LDH_D_e_coli_2 <= 1000
 0 <= R_MALS_e_coli_2 <= 1000
 0 <= R_MALt2_2_e_coli_2 <= 1000
-1000 <= R_MDH_e_coli_2 <= 1000
 0 <= R_ME1_e_coli_2 <= 1000
 0 <= R_ME2_e_coli_2 <= 1000
 0 <= R_NADH16_e_coli_2 <= 1000
 0 <= R_NADTRHD_e_coli_2 <= 1000
-1000 <= R_NH4t_e_coli_2 <= 1000
-1000 <= R_O2t_e_coli_2 <= 1000
 0 <= R_PDH_e_coli_2 <= 1000
End
