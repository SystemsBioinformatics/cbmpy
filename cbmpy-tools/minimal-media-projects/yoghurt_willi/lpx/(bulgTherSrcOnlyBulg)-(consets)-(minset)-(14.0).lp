\ENCODING=ISO-8859-1
\Problem name: LBULSTR

Minimize
 MNAF: bin_R_ASPt6_bulg + bin_R_CYSt6_bulg + bin_R_TYRt6_bulg
       + bin_R_LYSt6_bulg + bin_R_GLUt6_bulg + bin_R_GLNabc_bulg
       + bin_R_ARGt_bulg + bin_R_HISt6_bulg + bin_R_PHEt6_bulg
       + bin_R_GLYt6_bulg + bin_R_VALt6_bulg + bin_R_THRt6_bulg
       + bin_R_ALA_Lt6_bulg + bin_R_SERt6_bulg + bin_R_ILEt6_bulg
       + bin_R_ARGabc_bulg + bin_R_TRPt6_bulg + bin_R_ASNt6_bulg
       + bin_R_PROt6_bulg + bin_R_LEUt6_bulg + bin_R_gln_export_bulg
       + bin_R_METt6_bulg + bin_R_VALabc_ther + bin_R_LEUt6_ther
       + bin_R_SERt6_ther + bin_R_VALt6_ther + bin_R_METabc_ther
       + bin_R_GLNabc_ther + bin_R_PROt6_ther + bin_R_ALA_Lt6_ther
       + bin_R_TRPt6_ther + bin_R_GLYt6_ther + bin_R_GLUt6_ther
       + bin_R_ASPabc_ther + bin_R_GLUabc_ther + bin_R_LYSt6_ther
       + bin_R_PHEt6_ther + bin_R_CYSt6_ther + bin_R_HISt6_ther
       + bin_R_LEUabc_ther + bin_R_TYRt6_ther + bin_R_ARGabc_ther
       + bin_R_ILEt6_ther + bin_R_THRt6_ther + bin_R_ARGt_ther
       + bin_R_ILEabc_ther + bin_R_METt6_ther
Subject To
 M_10fthf_c_bulg:        - R_GARFT_bulg - R_AICART_bulg + R_MTHFC_bulg
                         + R_FTHFL_bulg  = 0
 M_10fthf_c_ther:        R_MTHFC_ther - R_GARFT_ther - R_AICART_ther
                         + R_FTHFL_ther  = 0
 M_12dgr_LLA_c_ther:     0.01 R_PAP_LLA_ther - 0.01 R_DAGGT_LLA_ther
                         - 0.01 R_DAGK_LLA_ther + 0.16 R_LTAS_LLA_ther
                         - 0.01 R_UGT1_LLA_ther  = 0
 M_13dpg_c_bulg:         R_GAPD_bulg - R_PGK_bulg  = 0
 M_13dpg_c_ther:         R_GAPD_ther - R_PGK_ther  = 0
 M_15dap_c_bulg:         R_LYSDC_bulg  = 0
 M_1pyr5c_c_ther:        R_G5SADs_ther - R_P5CRr_ther  = 0
 M_23dhdp_c_ther:        R_DHDPS_ther - R_DHDPRy_ther  = 0
 M_23dhmb_c_ther:        - R_KARA1_ther - R_DHAD1_ther  = 0
 M_23dhmp_c_ther:        - R_DHAD3_ther + R_KARA2_ther - R_DHAD2_ther  = 0
 M_25aics_c_bulg:        - R_ADSL2_bulg + R_PRASCS_bulg  = 0
 M_25aics_c_ther:        R_PRASCS_ther - R_ADSL2r_ther  = 0
 M_26dap_LL_c_ther:      - R_DAPE_ther + R_ADPDS_ther  = 0
 M_26dap_M_c_bulg:       - R_DAPDC_bulg - R_UAAGDS_bulg  = 0
 M_26dap_M_c_ther:       - R_DAPDC_ther + R_DAPE_ther - R_UAAGDS_ther  = 0
 M_2aeppn_c_bulg:        R_AEPabc_bulg  = 0
 M_2aeppn_e:             - R_AEPabc_bulg  = 0
 M_2ahbut_c_ther:        - R_KARA2_ther + R_ACHBS_ther  = 0
 M_2ahhmd_c_bulg:        - R_DHPS3_bulg  = 0
 M_2ahhmd_c_ther:        - R_DHPS3_ther + R_HPPK_ther  = 0
 M_2ahhmp_c_bulg:        R_DHNPA_bulg  = 0
 M_2ahhmp_c_ther:        - R_HPPK_ther + R_DHNPA_ther  = 0
 M_2beacp_c_bulg:        - R_BTMAT1_bulg + R_HBUHL1_bulg  = 0
 M_2beacp_c_ther:        R_HBUHL1_ther - R_BTMAT1_ther  = 0
 M_2chdeacp_c_bulg:      - 0.12 R_AGAT_LPL_bulg + R_FABM1_bulg
                         - 0.12 R_GAT1_LPL_bulg  = 0
 M_2chdeacp_c_ther:      - 0.03 R_AGAT_LLA_ther - 0.03 R_GAT1_LLA_ther
                         + R_FABM1_ther  = 0
 M_2cocdacp_c_bulg:      - 0.32 R_AGAT_LPL_bulg - R_CFAS180_bulg + R_FABM2_bulg
                         - 0.32 R_GAT1_LPL_bulg  = 0
 M_2cocdacp_c_ther:      R_FABM2_ther - R_CFAS180_ther - 0.44 R_AGAT_LLA_ther
                         - 0.44 R_GAT1_LLA_ther  = 0
 M_2cpr5p_c_bulg:        R_PRAI_bulg  = 0
 M_2cpr5p_c_ther:        - R_IGPS_ther + R_PRAI_ther  = 0
 M_2ctdeacp_c_ther:      - 0.005 R_AGAT_LLA_ther - 0.005 R_GAT1_LLA_ther
                         + R_FABM_ther  = 0
 M_2dda7p_c_ther:        - R_DHQS_ther + R_DAHPS_ther  = 0
 M_2dhp_c_ther:          - R_DPR_ther  = 0
 M_2dr1p_c_bulg:         R_PUNP4_bulg + R_PUNP2_bulg + R_PUNP6_bulg
                         + R_DURIPP_bulg  = 0
 M_2dr1p_c_ther:         R_PUNP2_ther + R_PUNP4_ther + R_DURIPP_ther
                         + R_PUNP6_ther - R_PPM2_ther  = 0
 M_2dr5p_c_bulg:         R_DRBK_bulg  = 0
 M_2dr5p_c_ther:         - R_DRPA_ther + R_PPM2_ther  = 0
 M_2hxic_L_c_bulg:       R_2HXICDH_bulg  = 0
 M_2hxic_L_c_ther:       R_2HXICDH_ther  = 0
 M_2ippm_c_ther:         - R_IPPMIb_ther + R_IPPMIa_ther  = 0
 M_2mahmp_c_ther:        R_PMPK_ther  = 0
 M_2mbcoa_c_bulg:        - R_KAS4_bulg - R_KAS10_bulg - R_KAS12_bulg  = 0
 M_2mp2coa_c_ther:       - R_ECOAH2c_ther  = 0
 M_2o3mpt_c_ther:        R_DHAD3_ther  = 0
 M_2obut_c_bulg:         R_CYSTGL_bulg  = 0
 M_2obut_c_ther:         R_THRD_L_ther + R_CYSTGL_ther - R_ACHBS_ther  = 0
 M_2pg_c_bulg:           - R_ENO_bulg + R_PGM_bulg  = 0
 M_2pg_c_ther:           R_PGM_ther - R_ENO_ther  = 0
 M_2pglyc_c_bulg:        - R_PGLYCP_bulg  = 0
 M_2shchc_c_bulg:        - R_SUCBZS_bulg  = 0
 M_2tddacp_c_bulg:       R_HDDHL5_bulg - R_DDMAT5_bulg  = 0
 M_2tddacp_c_ther:       R_HDDHL5_ther - R_DDMAT5_ther  = 0
 M_2tdeacp_c_bulg:       - R_DEMAT4_bulg + R_HDEHL4_bulg  = 0
 M_2tdeacp_c_ther:       R_HDEHL4_ther - R_DEMAT4_ther  = 0
 M_2thdeacp_c_bulg:      - R_HDMAT7_bulg + R_HHDHL7_bulg - R_FABM1_bulg  = 0
 M_2thdeacp_c_ther:      R_HHDHL7_ther - R_HDMAT7_ther - R_FABM1_ther  = 0
 M_2theacp_c_bulg:       - R_HEMAT2_bulg + R_HHYHL2_bulg  = 0
 M_2theacp_c_ther:       R_HHYHL2_ther - R_HEMAT2_ther  = 0
 M_2tocdacp_c_bulg:      - R_OCDMAT8_bulg - R_FABM2_bulg + R_HODHL8_bulg  = 0
 M_2tocdacp_c_ther:      - R_FABM2_ther + R_HODHL8_ther - R_OCDMAT8_ther  = 0
 M_2toceacp_c_bulg:      - R_OCMAT3_bulg + R_HOCHL3_bulg  = 0
 M_2toceacp_c_ther:      - R_OCMAT3_ther + R_HOCHL3_ther  = 0
 M_2ttdeacp_c_bulg:      - R_TDMAT6_bulg + R_HTDHL6_bulg  = 0
 M_2ttdeacp_c_ther:      - R_TDMAT6_ther + R_HTDHL6_ther - R_FABM_ther  = 0
 M_34hpp_c_bulg:         R_TYRTA_bulg  = 0
 M_34hpp_c_ther:         R_ASPTA5_ther + R_PPND_ther  = 0
 M_3c2hmp_c_ther:        - R_IPMD_ther - R_IPPMIa_ther  = 0
 M_3c3hmp_c_ther:        R_IPPS_ther + R_IPPMIb_ther  = 0
 M_3c4mop_c_ther:        - R_OMCDC_ther + R_IPMD_ther  = 0
 M_3dhq_c_ther:          R_DHQS_ther - R_DHQD_ther  = 0
 M_3dhsk_c_ther:         - R_SHK3D_ther + R_DHQD_ther  = 0
 M_3hbacp_c_bulg:        R_HBUR1_bulg - R_HBUHL1_bulg  = 0
 M_3hbacp_c_ther:        - R_HBUHL1_ther + R_HBUR1_ther  = 0
 M_3hddacp_c_bulg:       - R_HDDHL5_bulg + R_HDDR5_bulg  = 0
 M_3hddacp_c_ther:       - R_HDDHL5_ther + R_HDDR5_ther  = 0
 M_3hdeacp_c_bulg:       R_HDER4_bulg - R_HDEHL4_bulg  = 0
 M_3hdeacp_c_ther:       R_HDER4_ther - R_HDEHL4_ther  = 0
 M_3hhacp_c_bulg:        R_HHYR2_bulg - R_HHYHL2_bulg  = 0
 M_3hhacp_c_ther:        R_HHYR2_ther - R_HHYHL2_ther  = 0
 M_3hocacp_c_bulg:       R_HOCR3_bulg - R_HOCHL3_bulg  = 0
 M_3hocacp_c_ther:       R_HOCR3_ther - R_HOCHL3_ther  = 0
 M_3hocdacp_c_bulg:      - R_HODHL8_bulg + R_HODR8_bulg  = 0
 M_3hocdacp_c_ther:      - R_HODHL8_ther + R_HODR8_ther  = 0
 M_3hpaacp_c_bulg:       - R_HHDHL7_bulg + R_HHDR7_bulg  = 0
 M_3hpaacp_c_ther:       - R_HHDHL7_ther + R_HHDR7_ther  = 0
 M_3htdacp_c_bulg:       R_HTDR6_bulg - R_HTDHL6_bulg  = 0
 M_3htdacp_c_ther:       R_HTDR6_ther - R_HTDHL6_ther  = 0
 M_3ig3p_c_ther:         R_IGPS_ther - R_TRPS1_ther  = 0
 M_3mob_c_bulg:          R_VALTA_bulg  = 0
 M_3mob_c_ther:          R_VALTA_ther + R_DHAD1_ther - R_IPPS_ther  = 0
 M_3mop_c_bulg:          R_ILETA_bulg  = 0
 M_3mop_c_ther:          R_ILETA_ther + R_DHAD2_ther  = 0
 M_3oxddacp_c_bulg:      R_MCMAT5_bulg - R_HDDR5_bulg  = 0
 M_3oxddacp_c_ther:      - R_HDDR5_ther + R_MCMAT5_ther  = 0
 M_3oxdeacp_c_bulg:      R_MCMAT4_bulg - R_HDER4_bulg  = 0
 M_3oxdeacp_c_ther:      - R_HDER4_ther + R_MCMAT4_ther  = 0
 M_3oxhacp_c_bulg:       R_MCMAT2_bulg - R_HHYR2_bulg  = 0
 M_3oxhacp_c_ther:       R_MCMAT2_ther - R_HHYR2_ther  = 0
 M_3oxhdacp_c_bulg:      R_MCMAT7_bulg - R_HHDR7_bulg  = 0
 M_3oxhdacp_c_ther:      - R_HHDR7_ther + R_MCMAT7_ther  = 0
 M_3oxocacp_c_bulg:      - R_HOCR3_bulg + R_MCMAT3_bulg  = 0
 M_3oxocacp_c_ther:      - R_HOCR3_ther + R_MCMAT3_ther  = 0
 M_3oxocdacp_c_bulg:     R_MCMAT8_bulg - R_HODR8_bulg  = 0
 M_3oxocdacp_c_ther:     R_MCMAT8_ther - R_HODR8_ther  = 0
 M_3oxtdacp_c_bulg:      - R_HTDR6_bulg + R_MCMAT6_bulg  = 0
 M_3oxtdacp_c_ther:      R_MCMAT6_ther - R_HTDR6_ther  = 0
 M_3pg_c_bulg:           - R_PGCD_bulg + R_GLYCK_bulg + R_PGK_bulg - R_PGM_bulg
                          = 0
 M_3pg_c_ther:           R_GLYCK_ther - R_PGM_ther + R_PGK_ther - R_PGCD_ther
                          = 0
 M_3php_c_bulg:          R_PGCD_bulg  = 0
 M_3php_c_ther:          - R_PSERT_ther + R_PGCD_ther  = 0
 M_3psme_c_ther:         - R_CHORS_ther + R_PSCVT_ther  = 0
 M_4abz_c_bulg:          - R_DHPS3_bulg  = 0
 M_4abz_c_ther:          R_ADCL_ther - R_DHPS3_ther  = 0
 M_4adcho_c_ther:        - R_ADCL_ther + R_ADCS_ther  = 0
 M_4ahmmp_c_ther:        R_AHMMPS_ther - R_HMPK1_ther  = 0
 M_4ampm_c_ther:         - R_PMPK_ther + R_HMPK1_ther  = 0
 M_4h2oxg_c_ther:        R_ASPTA2_ther  = 0
 M_4hglu_c_ther:         - R_ASPTA2_ther  = 0
 M_4izp_c_ther:          R_URCN_ther  = 0
 M_4mop_c_bulg:          - R_2HXICDH_bulg + R_LEUTA_bulg  = 0
 M_4mop_c_ther:          R_OMCDC_ther - R_2HXICDH_ther + R_LEUTA_ther  = 0
 M_4mpetz_c_bulg:        R_THZPSN2_bulg  = 0
 M_4pasp_c_bulg:         R_ASAD_bulg + R_ASPK_bulg  = 0
 M_4pasp_c_ther:         R_ASPK_ther + R_ASAD_ther  = 0
 M_4ppan_c_bulg:         R_PNTK_bulg - R_PPNCL_bulg  = 0
 M_4ppan_c_ther:         - R_PPNCL_ther + R_PNTK_ther  = 0
 M_4ppcys_c_bulg:        - R_PPCDC_bulg + R_PPNCL_bulg  = 0
 M_4ppcys_c_ther:        R_PPNCL_ther - R_PPCDC_ther  = 0
 M_5aizc_c_bulg:         - R_PRASCS_bulg + R_AIRC_bulg  = 0
 M_5aizc_c_ther:         R_AIRC_ther - R_PRASCS_ther  = 0
 M_5dpmev_c_bulg:        R_PMEVK_bulg - R_DPMVD_bulg  = 0
 M_5dpmev_c_ther:        - R_DPMVD_ther + R_PMEVK_ther  = 0
 M_5fothf_c_bulg:        - R_FTHFCL_bulg  = 0
 M_5fothf_c_ther:        - R_FTHFCL_ther  = 0
 M_5mta_c_bulg:          - R_MTAN_bulg  = 0
 M_5mta_c_ther:          - R_MTAN_ther  = 0
 M_5mthf_c_bulg:         - R_METFR_bulg - R_METS_bulg  = 0
 M_5mthf_c_ther:         - R_METFR_ther - R_METS_ther  = 0
 M_5mthglu_c_bulg:       - R_MTHPTGHM_bulg  = 0
 M_5mthglu_c_ther:       - R_MTHPTGHM_ther  = 0
 M_5mtr_c_bulg:          R_MTAN_bulg  = 0
 M_5mtr_c_ther:          R_MTAN_ther  = 0
 M_5pmev_c_bulg:         R_MEVK_bulg - R_PMEVK_bulg  = 0
 M_5pmev_c_ther:         R_MEVK_ther - R_PMEVK_ther  = 0
 M_6pgc_c_bulg:          - R_PGDH_bulg + R_PGL_bulg  = 0
 M_6pgg_c_bulg:          - R_PGGH_bulg + R_CELBpts_bulg  = 0
 M_6pgg_c_ther:          - R_PGGH_ther + R_CELBpts_ther  = 0
 M_6pgl_c_bulg:          - R_PGL_bulg + R_G6PDHy_bulg  = 0
 M_6pthp_c_ther:         R_PTHPS_ther  = 0
 M_CPS_LLA_c_ther:       - 0.248041 R_biomass_STR_ther + R_CPSS_LLA_ther  = 0
 M_CPS_LPL_c_bulg:       - 0.1294 R_biomass_LBUL10_bulg + R_CPSS_LPL_bulg  = 0
 M_DNA_LPL_c_bulg:       0.01 R_DNAS_LPL_bulg - 0.000617 R_biomass_LBUL10_bulg
                          = 0
 M_DNA_STU_c_ther:       0.001 R_DNAS_STU_ther - 0.000032 R_biomass_STR_ther
                          = 0
 M_GTAglc_c_bulg:        - 0.01378 R_biomass_LBUL10_bulg + R_GTAPGL_bulg  = 0
 M_LTAAlaGal_LLA_c_ther: - 0.000148 R_biomass_STR_ther + 0.01 R_GALTAL_ther
                          = 0
 M_LTA_LLA_c_ther:       - 0.01 R_DALTAL_LLA_ther + 0.01 R_LTAS_LLA_ther  = 0
 M_LTA_LPL_c_bulg:       - 0.01 R_DALTAL_bulg  = 0
 M_LTA_LPL_c_ther:       - 0.01 R_GLTAL_ther  = 0
 M_LTAala_LLA_c_ther:    0.01 R_DALTAL_LLA_ther - 0.01 R_GALTAL_ther  = 0
 M_LTAala_LPL_c_bulg:    0.01 R_DALTAL_bulg  = 0
 M_LTAglc_LPL_c_ther:    0.01 R_GLTAL_ther  = 0
 M_Lglyald_c_bulg:       - R_ALCD19_L_bulg  = 0
 M_MCOOH_c_bulg:         R_THZPSN2_bulg - R_MOADCST_bulg  = 0
 M_MCOOH_c_ther:         - R_MOADCST_ther  = 0
 M_MCOSH_c_bulg:         - R_THZPSN2_bulg + R_MOADCST_bulg  = 0
 M_MCOSH_c_ther:         R_MOADCST_ther  = 0
 M_PG_c_ther:            - 0.12447 R_biomass_STR_ther + R_PGGT2_ther  = 0
 M_PGlac2_c_bulg:        R_PGGT4_bulg - 0.1462 R_biomass_LBUL10_bulg  = 0
 M_PROT_LPL_v60_c_bulg:  - 0.00243 R_biomass_LBUL10_bulg
                         + 0.001 R_PROTS_LPL_v60_bulg  = 0
 M_PROT_STR_c_ther:      0.001 R_PROTS_STR_ther - 0.002332 R_biomass_STR_ther
                          = 0
 M_PreZ_c_ther:          R_PZS_ther  = 0
 M_RNA_LLA_c_ther:       - 0.002674 R_biomass_STR_ther + 0.01 R_RNAS_LLA_ther
                          = 0
 M_RNA_LPL_c_bulg:       0.01 R_RNAS_LPL_bulg - 0.002784 R_biomass_LBUL10_bulg
                          = 0
 M_RTA_c_bulg:           - R_GRTAL_bulg - R_DARTAL_bulg  = 0
 M_RTA_c_ther:           - R_GRTAL_ther  = 0
 M_RTAala_c_bulg:        R_DARTAL_bulg  = 0
 M_RTAglc_c_bulg:        R_GRTAL_bulg  = 0
 M_RTAglc_c_ther:        R_GRTAL_ther  = 0
 M_a_gal_D_c_ther:       R_GALM_ther - R_GALK2_ther  = 0
 M_aaacp_c_bulg:         - R_HBUR1_bulg + R_kaasIII_bulg  = 0
 M_aaacp_c_ther:         - R_HBUR1_ther + R_kaasIII_ther  = 0
 M_aacoa_c_bulg:         R_ACACT1r_bulg + R_HMGCOAS_bulg  = 0
 M_aacoa_c_ther:         R_HMGCOAS_ther + R_ACACT1_ther  = 0
 M_ac_c_bulg:            R_ACt6_bulg - R_ACKr_bulg + R_CYSS_bulg + R_AGDC_bulg
                          = 0
 M_ac_c_ther:            - R_ACKr_ther + R_ACTPASE_ther + R_AHSERL3_ther
                         + R_METACH_ther + R_SLCYSS_ther + R_AGDC_ther
                         + R_ACt6_ther + R_CYSS_ther + R_AHSERL4_ther
                         + R_ADPDS_ther + R_SHSL5_ther  = 0
 M_ac_e:                 - R_ACt6_bulg - R_EX_ac_e_ - R_ACt6_ther  = 0
 M_acald_c_bulg:         - R_ACALD_bulg + R_ACALDt_bulg + R_ALCD2x_bulg  = 0
 M_acald_c_ther:         R_DRPA_ther + R_ALCD2x_ther + R_THRA_ther  = 0
 M_acald_e:              - R_ACALDt_bulg  = 0
 M_acamoxm_c_ther:       R_ADPTA_ther + R_APAT_ther  = 0
 M_accoa_c_bulg:         - R_G1PACT_bulg - R_ACCOAC_bulg + R_ACALD_bulg
                         - R_HSAT_bulg - 2 R_ACACT1r_bulg - R_KAS2_bulg
                         - R_PTAr_bulg - R_KAS7_bulg - R_KAS13_bulg
                         - R_KAS8_bulg - R_kaasIII_bulg + R_HMGCOAS_bulg
                         - R_KAS17_bulg  = 0
 M_accoa_c_ther:         - R_ACCOAC_ther + R_HMGCOAS_ther - R_G1PACT_ther
                         - R_PTAr_ther + R_PFL_ther - R_SERAT_ther
                         - 2 R_ACACT1_ther - R_HSAT_ther - R_CS_ther
                         + R_PDH_ther - R_APAT_ther - R_IPPS_ther
                         - R_kaasIII_ther  = 0
 M_acetone_c_ther:       2 R_NPPDO_ther  = 0
 M_acg5p_c_ther:         R_AGPR_ther + R_ACGK_ther  = 0
 M_acg5sa_c_ther:        - R_AGPR_ther + R_ACOTA_ther  = 0
 M_acgal6p_c_bulg:       R_ACGALpts_bulg  = 0
 M_acgal6p_c_ther:       R_ACGALpts_ther  = 0
 M_acgala_e:             - R_ACGALpts_ther - R_ACGALpts_bulg  = 0
 M_acgam1p_c_bulg:       R_G1PACT_bulg - R_UAGDP_bulg  = 0
 M_acgam1p_c_ther:       R_G1PACT_ther - R_UAGDP_ther  = 0
 M_acgam6p_c_bulg:       R_ACGApts_bulg - R_AGDC_bulg  = 0
 M_acgam6p_c_ther:       R_ACGApts_ther - R_AGDC_ther  = 0
 M_acgam_e:              - R_ACGApts_ther - R_ACGApts_bulg  = 0
 M_acglu_c_ther:         R_ORNTAC_ther - R_ACGK_ther  = 0
 M_achms_c_bulg:         R_HSAT_bulg  = 0
 M_achms_c_ther:         - R_METACH_ther + R_HSAT_ther - R_SHSL5_ther  = 0
 M_acorn_c_ther:         - R_ORNTAC_ther - R_ACOTA_ther  = 0
 M_acp_c_bulg:           R_AGAT_LPL_bulg + R_ACPS2_bulg + R_MCMAT2_bulg
                         + R_MCMAT5_bulg + R_MCMAT8_bulg + R_MCMAT7_bulg
                         - R_MACPMT_bulg + R_MCMAT3_bulg + R_MCMAT4_bulg
                         + R_MCMAT6_bulg + R_GAT1_LPL_bulg  = 0
 M_acp_c_ther:           R_MCMAT2_ther + R_MCMAT6_ther - R_MACPMT_ther
                         + R_MCMAT8_ther + R_MCMAT3_ther + R_AGAT_LLA_ther
                         + R_ACPS2_ther + R_MCMAT4_ther + R_GAT1_LLA_ther
                         + R_MCMAT5_ther + R_MCMAT7_ther  = 0
 M_acser_c_bulg:         - R_CYSS_bulg  = 0
 M_acser_c_ther:         - R_AHSERL3_ther + R_SERAT_ther - R_SLCYSS_ther
                         - R_CYSS_ther - R_AHSERL4_ther  = 0
 M_actn_R_c_ther:        R_ACLDC_ther  = 0
 M_actp_c_bulg:          R_PTAr_bulg + R_ACKr_bulg + R_PKL_bulg + R_PYROX_bulg
                          = 0
 M_actp_c_ther:          R_ACKr_ther - R_ACTPASE_ther + R_PTAr_ther
                         + R_PYROX_ther  = 0
 M_ade_c_bulg:           R_MTAN_bulg + R_PUNP1_bulg + R_PUNP2_bulg
                         - R_ADPT_bulg + R_AHCYSNS_bulg + R_PNS1_bulg  = 0
 M_ade_c_ther:           R_PUNP2_ther + R_AHCYSNS_ther - R_ADPT_ther
                         + R_MTAN_ther + R_PUNP1_ther  = 0
 M_adn_c_bulg:           - R_ADA_bulg - R_PUNP1_bulg - R_PNS1_bulg  = 0
 M_adn_c_ther:           R_NTD7_ther - R_ADA_ther - R_PUNP1_ther  = 0
 M_adp_c_bulg:           0.4 R_RNAS_LPL_bulg + R_DADNK_bulg + 25 R_DALTAL_bulg
                         + R_SHKK_bulg + R_DRBK_bulg + 1.37 R_DNAS_LPL_bulg
                         + R_MNabc_bulg + R_PRFGS_bulg + R_ACCOAC_bulg
                         + R_UAAGLS1_bulg + R_RBK_bulg + 2 R_PEPTabc_bulg
                         + R_ALALAC_bulg + R_GLYBabc_bulg + R_THFGLUS_bulg
                         + R_HSK_bulg + R_GLNS_bulg + R_NDPK1_bulg
                         + R_AEPabc_bulg + R_DGNSK_bulg + R_CTPS1_1_bulg
                         + R_UGMDDS2_bulg + R_CYTK1_bulg + R_MEVK_bulg
                         + R_GLYK_bulg + R_CHLabc_bulg + R_PRASCS_bulg
                         + R_NADK_bulg + R_NDPK3_bulg + R_DHFS_bulg
                         + R_HEX1_bulg + R_PFK_bulg + R_THMabc_bulg
                         + R_UDCPK_bulg + R_GLNabc_bulg + R_PRAGS_bulg
                         + R_GK1_bulg - R_PYK_bulg + R_PIabc_bulg
                         + R_NDPK6_bulg - R_ATPS3r_bulg
                         + 27.2 R_biomass_LBUL10_bulg + R_NDPK8_bulg
                         + R_ZN2abc2_bulg + R_UAMAS_bulg + R_UAMAGS_bulg
                         + R_CYTDK1_bulg + R_GLNTAL_bulg + R_GLYCK_bulg
                         + R_ADK2_bulg + R_PYDAMK_bulg + R_ACKr_bulg
                         + R_PNTK_bulg + R_DINSK_bulg + R_DPCOAK_bulg
                         + R_PRAIS_bulg + R_DURIK1_bulg + R_CTPS2_bulg
                         + R_NDPK9_bulg + 2 R_ADK1_bulg - R_PGK_bulg
                         + R_PMEVK_bulg + R_FTHFCL_bulg + R_NDPK5_bulg
                         + R_ASPK_bulg + R_DPMVD_bulg + R_PTRCabc_bulg
                         + R_SPMDabc_bulg + R_URIDK2_bulg + 0.111 R_CASabc_bulg
                         + R_MALTabc_bulg + R_NDPK7_bulg + R_NDPK2_bulg
                         + R_TMDK1_bulg + R_UAAGDS_bulg + R_DTMPK_bulg
                         + R_ARGabc_bulg + R_ATPM_bulg + R_URIK1_bulg
                         + 0.306 R_PROTS_LPL_v60_bulg + R_NAabcO_bulg
                         + R_PYDXK_bulg + 25 R_DARTAL_bulg + R_PPIK_bulg
                         + R_URIDK1_bulg + R_CYTK2_bulg + R_FRUK_bulg
                         + R_FTHFL_bulg + 2 R_CBPS_bulg + R_F1PK_bulg
                         + R_NDPK4_bulg + R_RBFK_bulg  = 0
 M_adp_c_ther:           R_MALTTTRabc_ther + R_VALabc_ther + R_PIabc_ther
                         + R_ACCOAC_ther + R_URIDK2_ther + R_GLYCK_ther
                         + R_PRAGS_ther + R_ACKr_ther + R_RBK_ther
                         + R_F1PK_ther + R_NAMNPP_ther + R_DPMVD_ther
                         + R_ADK2_ther + R_PYDAMK_ther + R_CTPS2_ther
                         + R_CYTK1_ther + R_SPMDabc_ther + R_PYDXK_ther
                         + 1.37 R_DNAS_STU_ther + 2 R_ADK1_ther + R_MEVK_ther
                         + R_PMPK_ther + R_URIK1_ther + 0.306 R_PROTS_STR_ther
                         - R_ATPS3r_ther + R_HSK_ther + R_ASPabc_ther
                         + R_MALTHXabc_ther + R_MALTPTabc_ther + R_PRFGS_ther
                         + 2 R_PEPTabc_ther + R_ZN2abc2_ther + R_METabc_ther
                         + R_PRAIS_ther + R_PNTK_ther + R_GLNabc_ther
                         + R_GLNS_ther + R_PTRCabc_ther + R_CYTK2_ther
                         + R_GLU5K_ther + R_THMabc_ther + R_METSRabc_ther
                         + R_MNabc_ther + R_NDPK8_ther + R_NADK_ther
                         + R_CYTDK1_ther + R_GLYBabc_ther + R_NDPK4_ther
                         + R_FRUK_ther + R_DAGK_LLA_ther + R_DHFS_ther
                         + R_UAMAGS_ther + R_NDPK2_ther + R_THFGLUS_ther
                         + R_FERRICHabc_ther + R_METSabc_ther + R_NDPK7_ther
                         + R_NDPK3_ther + R_DURIK1_ther
                         + 18.15 R_biomass_STR_ther + R_PMEVK_ther
                         + 0.4 R_RNAS_LLA_ther + 6 R_DALTAL_LLA_ther
                         + R_NDPK9_ther + R_SHKK_ther - R_PGK_ther
                         + R_UAAGDS_ther + R_RBFK_ther + R_GLUabc_ther
                         + R_UAAGLS1_ther - R_RNDR1_ther + R_NDPK1_ther
                         + R_UGLDDS1_ther + R_ALAALA_ther + R_NDPK6_ther
                         + R_PFK_ther + R_HEX1_ther - R_PYK_ther
                         + R_MALTTRabc_ther + R_TMDK1_ther + R_URIDK1_ther
                         + R_UDCPK_ther + R_ATPM_ther + R_PROLBabc_ther
                         + R_LEUabc_ther + R_ARGabc_ther + 0.111 R_CASabc_ther
                         + R_DPCOAK_ther + R_ASPK_ther + R_PRASCS_ther
                         + 2 R_CBPS_ther + R_GALK2_ther + R_METDabc_2_ther
                         + R_UAMAS_ther + R_FTHFCL_ther + R_ACGK_ther
                         + R_DTMPK_ther + R_MALTabc_ther + R_NDPK5_ther
                         + R_GK1_ther + R_GLNTAL_ther + R_HMPK1_ther
                         + R_CTPS1_1_ther + R_FTHFL_ther + R_ILEabc_ther  = 0
 M_adprib_c_bulg:        - R_ADPRDP_bulg  = 0
 M_adprib_c_ther:        - R_ADPRDP_ther  = 0
 M_agly3p_LLA_c_ther:    - 0.01 R_AGAT_LLA_ther + 0.01 R_GAT1_LLA_ther  = 0
 M_agly3p_LPL_c_bulg:    - 0.01 R_AGAT_LPL_bulg + 0.01 R_GAT1_LPL_bulg  = 0
 M_agm_c_bulg:           R_ARGDC_bulg  = 0
 M_agm_c_ther:           - R_AGMT_ther  = 0
 M_ahcys_c_bulg:         R_CFAS180_bulg + R_HCYSMT_bulg - R_AHCYSNS_bulg  = 0
 M_ahcys_c_ther:         R_CFAS180_ther - R_AHCYSNS_ther + R_HCYSMT_ther  = 0
 M_ahdt_c_bulg:          R_GTPCI_bulg - R_DNTPPA_bulg  = 0
 M_ahdt_c_ther:          - R_DNTPPA_ther - R_PTHPS_ther + R_GTPCI_ther  = 0
 M_aicar_c_bulg:         R_ADSL2_bulg - R_AICART_bulg  = 0
 M_aicar_c_ther:         - R_AICART_ther + R_ADSL2r_ther  = 0
 M_air_c_bulg:           R_PRAIS_bulg - R_AIRC_bulg  = 0
 M_air_c_ther:           R_PRAIS_ther - R_AHMMPS_ther - R_AIRC_ther  = 0
 M_akg_c_bulg:           - R_VALTA_bulg - R_AKGTA3_bulg - R_PHETA1_bulg
                         - R_LEUTA_bulg - R_TYRTA_bulg - R_ILETA_bulg  = 0
 M_akg_c_ther:           R_ICDHx_ther - R_ADPTA_ther - R_SDPTA_ther
                         - R_ASPTA5_ther + R_PSERT_ther + R_GLUDy_ther
                         - R_ASPTA6_ther - R_VALTA_ther - R_ILETA_ther
                         - R_ASPTA1_ther - R_ASPTA4_ther - R_AKGTA3_ther
                         + R_ICDHy_ther - R_ALATA_L_ther - R_ACOTA_ther
                         - R_LEUTA_ther - R_ASPTA2_ther  = 0
 M_al26da_c_ther:        - R_ADPTA_ther - R_ADPDS_ther  = 0
 M_ala_D_c_bulg:         - 25 R_DALTAL_bulg - R_ALALAC_bulg + R_ALAR_bulg
                         - 25 R_DARTAL_bulg  = 0
 M_ala_D_c_ther:         R_ALAR_ther - 6 R_DALTAL_LLA_ther - 2 R_ALAALA_ther
                          = 0
 M_ala_L_c_bulg:         - R_ALATRS_bulg - R_ALAR_bulg + M_ala_L_c_bulg_src
                         - R_UAMAS_bulg + 0.0763611709802577 R_CASL_bulg
                         + 0.032 R_PEPTL_bulg + R_ALA_Lt6_bulg
                         - R_ALATA_Lr_bulg  = 0
 M_ala_L_c_ther:         - R_ALAR_ther - R_ALATRS_ther - R_ALAD_Lr_ther
                         + R_ALA_Lt6_ther + 0.0763611709802577 R_CASL_ther
                         + 0.032 R_PEPTL_ther - R_UAMAS_ther - R_ALATA_L_ther
                          = 0
 M_ala_L_e:              - R_EX_ala_L_e_ - R_ALA_Lt6_ther - R_ALA_Lt6_bulg  = 0
 M_alaala_c_ther:        - R_UGLDDS1_ther + R_ALAALA_ther  = 0
 M_alac_S_c_ther:        - R_ACLDC_ther + R_ACLS_ther + R_KARA1_ther  = 0
 M_alalac_c_bulg:        R_ALALAC_bulg - R_UGMDDS2_bulg  = 0
 M_alatrna_c_bulg:       R_ALATRS_bulg - 0.125 R_PROTS_LPL_v60_bulg  = 0
 M_alatrna_c_ther:       - 0.082 R_PROTS_STR_ther + R_ALATRS_ther  = 0
 M_amet_c_bulg:          R_METAT_bulg - R_CFAS180_bulg - R_HCYSMT_bulg  = 0
 M_amet_c_ther:          R_METAT_ther - R_CFAS180_ther - R_HCYSMT_ther  = 0
 M_amp_c_bulg:           R_ALATRS_bulg + R_CYSTRS_bulg + R_ASNTRS_bulg
                         + R_ARGTRS_bulg + R_LEUTRS_bulg + R_SERTRS_bulg
                         + R_ADPRDP_bulg + R_ASPTRS_bulg + R_VALTRS_bulg
                         + R_METTRS_bulg + R_ASNS2_bulg + R_PROTRS_bulg
                         + R_PRPPS_bulg + R_HISTRS_bulg + R_MOADCST_bulg
                         + R_PHETRS_bulg + R_GLYTRS_bulg - R_ADK2_bulg
                         + R_GTPDPK_bulg + R_GMPS2_bulg - R_ADK1_bulg
                         + R_THRTRS_bulg + R_NADS1_bulg + R_GLUTRS_bulg
                         + R_ADPT_bulg + R_LYSTRS_bulg + R_TYRTRS_bulg
                         + R_TMDPK_bulg + R_ILETRS_bulg + R_ASNS1_bulg
                         + R_TRPTRS_bulg + R_ADSL1_bulg  = 0
 M_amp_c_ther:           R_GLUTRS_ther + R_HISTRS_ther + R_NADS1_ther
                         - R_NTD7_ther + R_ASNS2_ther + R_GTPDPK_ther
                         + R_ARGSS_ther + R_TMDPK_ther - R_ADK2_ther
                         - R_ADK1_ther + R_VALTRS_ther + R_CYSTRS_ther
                         + R_THRTRS_ther + R_ADPRDP_ther + R_ADSL1r_ther
                         + R_PHETRS_ther + R_ALATRS_ther + R_GMPS2_ther
                         + R_GLYTRS_ther + R_PRPPS_ther + R_MOADCST_ther
                         + R_ADPT_ther + R_ARGTRS_ther + R_ASNTRS_ther
                         + R_ILETRS_ther + R_SERTRS_ther + R_TYRTRS_ther
                         + R_HPPK_ther + R_LEUTRS_ther + R_LYSTRS_ther
                         + R_METTRS_ther + R_TRPTRS_ther + R_ASPTRS_ther
                         + R_PROTRS_ther  = 0
 M_anth_c_ther:          R_ANS1_ther - R_ANPRT_ther  = 0
 M_apoACP_c_bulg:        - R_ACPS2_bulg  = 0
 M_apoACP_c_ther:        - R_ACPS2_ther  = 0
 M_arg_L_c_bulg:         - R_ARGTRS_bulg + R_ARGt_bulg + M_arg_L_c_bulg_src
                         + 0.0109838387295756 R_CASL_bulg + 0.031 R_PEPTL_bulg
                         + R_ARGabc_bulg - R_ARGDC_bulg  = 0
 M_arg_L_c_ther:         0.0109838387295756 R_CASL_ther + R_ARGSL_ther
                         - R_ARGTRS_ther + R_ARGabc_ther + 0.031 R_PEPTL_ther
                         + R_ARGt_ther  = 0
 M_arg_L_e:              - R_ARGt_bulg - R_EX_arg_L_e_ - R_ARGabc_bulg
                         - R_ARGabc_ther - R_ARGt_ther  = 0
 M_argsuc_c_ther:        R_ARGSS_ther - R_ARGSL_ther  = 0
 M_argtrna_c_bulg:       R_ARGTRS_bulg - 0.04 R_PROTS_LPL_v60_bulg  = 0
 M_argtrna_c_ther:       - 0.023 R_PROTS_STR_ther + R_ARGTRS_ther  = 0
 M_asn_L_c_bulg:         - R_ASNTRS_bulg + R_ASNt6_bulg + R_ASNS2_bulg
                         + M_asn_L_c_bulg_src - R_ASNN_bulg + R_ASNS1_bulg  = 0
 M_asn_L_c_ther:         R_ASNS2_ther - R_ASNN_ther - R_ASNTRS_ther  = 0
 M_asn_L_e:              - R_EX_asn_L_e_ - R_ASNt6_bulg  = 0
 M_asntrna_c_bulg:       R_ASNTRS_bulg - 0.06 R_PROTS_LPL_v60_bulg  = 0
 M_asntrna_c_ther:       - 0.08 R_PROTS_STR_ther + R_ASNTRS_ther  = 0
 M_asp_D_c_bulg:         - R_ASPR_bulg  = 0
 M_asp_L_c_bulg:         R_ASPt6_bulg + M_asp_L_c_bulg_src - R_PRASCS_bulg
                         - R_ASPT_bulg - R_ASPTRS_bulg - R_ASNS2_bulg
                         - R_ADSS_bulg - R_ASPK_bulg
                         + 0.137523384267553 R_CASL_bulg + 0.076 R_PEPTL_bulg
                         - R_ASPCT_bulg + R_ASNN_bulg - R_ASNS1_bulg
                         + R_ASPR_bulg  = 0
 M_asp_L_c_ther:         - R_ASNS2_ther - R_ARGSS_ther - R_ASPCT_ther
                         + R_ASPabc_ther - R_ASPTA1_ther + R_ASNN_ther
                         + 0.137523384267553 R_CASL_ther - R_ADSS_ther
                         + 0.076 R_PEPTL_ther - R_ASPK_ther - R_PRASCS_ther
                         - R_ASPTRS_ther  = 0
 M_asp_L_e:              - R_ASPt6_bulg - R_EX_asp_L_e_ - R_ASPabc_ther  = 0
 M_aspsa_c_bulg:         R_HSDy_bulg - R_ASAD_bulg  = 0
 M_aspsa_c_ther:         - R_DHDPS_ther + R_HSDy_ther - R_ASAD_ther  = 0
 M_asptrna_c_bulg:       R_ASPTRS_bulg - 0.06 R_PROTS_LPL_v60_bulg  = 0
 M_asptrna_c_ther:       - 0.054 R_PROTS_STR_ther + R_ASPTRS_ther  = 0
 M_atp_c_bulg:           - 0.66 R_RNAS_LPL_bulg - R_DADNK_bulg
                         - 25 R_DALTAL_bulg - R_SHKK_bulg - R_METAT_bulg
                         - R_DRBK_bulg - 1.37 R_DNAS_LPL_bulg - R_MNabc_bulg
                         - R_ALATRS_bulg - R_CYSTRS_bulg - R_PRFGS_bulg
                         - R_ACCOAC_bulg - R_ASNTRS_bulg - R_UAAGLS1_bulg
                         - R_RBK_bulg - 2 R_PEPTabc_bulg - R_ARGTRS_bulg
                         - R_ALALAC_bulg - R_GLYBabc_bulg - R_THFGLUS_bulg
                         - R_HSK_bulg - R_NNAT_bulg - R_GLNS_bulg
                         - R_NDPK1_bulg - R_AEPabc_bulg - R_LEUTRS_bulg
                         - R_DGNSK_bulg - R_CTPS1_1_bulg - R_UGMDDS2_bulg
                         - R_CYTK1_bulg - R_MEVK_bulg - R_GLYK_bulg
                         - R_CHLabc_bulg - R_PRASCS_bulg - R_NADK_bulg
                         - R_NDPK3_bulg - R_DHFS_bulg - R_SERTRS_bulg
                         - R_HEX1_bulg - R_PFK_bulg - R_THMabc_bulg
                         - R_UDCPK_bulg - R_GLNabc_bulg - R_PRAGS_bulg
                         - R_GK1_bulg - R_FMNAT_bulg + R_PYK_bulg
                         - R_PIabc_bulg - R_ASPTRS_bulg - R_VALTRS_bulg
                         - R_NDPK6_bulg - R_METTRS_bulg + R_ATPS3r_bulg
                         - R_ASNS2_bulg - R_PROTRS_bulg
                         - 27.2 R_biomass_LBUL10_bulg - R_NDPK8_bulg
                         - R_ZN2abc2_bulg - R_UAMAS_bulg - R_PRPPS_bulg
                         - R_UAMAGS_bulg - R_HISTRS_bulg - R_CYTDK1_bulg
                         - R_MOADCST_bulg - R_PHETRS_bulg - R_GLYTRS_bulg
                         - R_GLNTAL_bulg - R_GLYCK_bulg - R_PYDAMK_bulg
                         - R_ACKr_bulg - R_PNTK_bulg - R_GTPDPK_bulg
                         - R_DINSK_bulg - R_GMPS2_bulg - R_DPCOAK_bulg
                         - R_PRAIS_bulg - R_DURIK1_bulg - R_CTPS2_bulg
                         - R_NDPK9_bulg - R_ADK1_bulg - R_PTPAT_bulg
                         + R_PGK_bulg - R_THRTRS_bulg - R_PMEVK_bulg
                         - R_FTHFCL_bulg - R_NDPK5_bulg - R_ASPK_bulg
                         - R_DPMVD_bulg - R_PTRCabc_bulg - R_NADS1_bulg
                         - R_GLUTRS_bulg - R_SPMDabc_bulg - R_URIDK2_bulg
                         - 0.111 R_CASabc_bulg - R_MALTabc_bulg - R_NDPK7_bulg
                         - R_NDPK2_bulg - R_LYSTRS_bulg - R_TYRTRS_bulg
                         - R_TMDK1_bulg - R_UAAGDS_bulg - R_DTMPK_bulg
                         - R_ARGabc_bulg - R_TMDPK_bulg - R_ATPM_bulg
                         - R_URIK1_bulg - 0.306 R_PROTS_LPL_v60_bulg
                         - R_NAabcO_bulg - R_PYDXK_bulg - R_ILETRS_bulg
                         - 25 R_DARTAL_bulg - R_ASNS1_bulg - R_TRPTRS_bulg
                         - R_PPIK_bulg - R_URIDK1_bulg - R_CYTK2_bulg
                         - R_FRUK_bulg - R_FTHFL_bulg - 2 R_CBPS_bulg
                         - R_F1PK_bulg - R_RNTR1_bulg - R_NDPK4_bulg
                         - R_RBFK_bulg  = 0
 M_atp_c_ther:           - R_MALTTTRabc_ther - R_GLUTRS_ther - R_HISTRS_ther
                         - R_VALabc_ther - R_METAT_ther - R_PIabc_ther
                         - R_ACCOAC_ther - R_URIDK2_ther - R_NADS1_ther
                         - R_PTPAT_ther - R_GLYCK_ther - R_PRAGS_ther
                         - R_ACKr_ther - R_ASNS2_ther - R_RBK_ther
                         - R_GTPDPK_ther - R_F1PK_ther - R_ARGSS_ther
                         - R_NAMNPP_ther - R_DPMVD_ther - R_TMDPK_ther
                         - R_PYDAMK_ther - R_CTPS2_ther - R_CYTK1_ther
                         - R_SPMDabc_ther - R_PYDXK_ther - 1.37 R_DNAS_STU_ther
                         - R_ADK1_ther - R_MEVK_ther - R_PMPK_ther
                         - R_URIK1_ther - R_VALTRS_ther
                         - 0.306 R_PROTS_STR_ther + R_ATPS3r_ther
                         - R_FMNAT_ther - R_HSK_ther - R_ASPabc_ther
                         - R_MALTHXabc_ther - R_MALTPTabc_ther - R_PRFGS_ther
                         - 2 R_PEPTabc_ther - R_ZN2abc2_ther - R_CYSTRS_ther
                         - R_METabc_ther - R_PRAIS_ther - R_PNTK_ther
                         - R_GLNabc_ther - R_GLNS_ther - R_PTRCabc_ther
                         - R_CYTK2_ther - R_THRTRS_ther - R_GLU5K_ther
                         - R_THMabc_ther - R_METSRabc_ther - R_RNTR1_ther
                         - R_MNabc_ther - R_NDPK8_ther - R_PHETRS_ther
                         - R_ALATRS_ther - R_NADK_ther - R_CYTDK1_ther
                         - R_GLYBabc_ther - R_NDPK4_ther - R_FRUK_ther
                         - R_DAGK_LLA_ther - R_DHFS_ther - R_UAMAGS_ther
                         - R_NDPK2_ther - R_THFGLUS_ther - R_GMPS2_ther
                         - R_FERRICHabc_ther - R_METSabc_ther - R_GLYTRS_ther
                         - R_NDPK7_ther - R_NDPK3_ther - R_PRPPS_ther
                         - R_DURIK1_ther - 18.15 R_biomass_STR_ther
                         - R_PMEVK_ther - 0.66 R_RNAS_LLA_ther
                         - 6 R_DALTAL_LLA_ther - R_NDPK9_ther - R_SHKK_ther
                         - R_MOADCST_ther + R_PGK_ther - R_UAAGDS_ther
                         - R_RBFK_ther - R_GLUabc_ther - R_UAAGLS1_ther
                         - R_NDPK1_ther - R_UGLDDS1_ther - R_ALAALA_ther
                         - R_NDPK6_ther - R_PFK_ther - R_HEX1_ther + R_PYK_ther
                         - R_MALTTRabc_ther - R_TMDK1_ther - R_ARGTRS_ther
                         - R_URIDK1_ther - R_UDCPK_ther - R_ATPM_ther
                         - R_ASNTRS_ther - R_ILETRS_ther - R_PROLBabc_ther
                         - R_LEUabc_ther - R_SERTRS_ther - R_ARGabc_ther
                         - R_TYRTRS_ther - R_HPPK_ther - 0.111 R_CASabc_ther
                         - R_LEUTRS_ther - R_DPCOAK_ther - R_ASPK_ther
                         - R_PRASCS_ther - 2 R_CBPS_ther - R_LYSTRS_ther
                         - R_METTRS_ther - R_GALK2_ther - R_METDabc_2_ther
                         - R_UAMAS_ther - R_FTHFCL_ther - R_NNAT_ther
                         - R_ACGK_ther - R_DTMPK_ther - R_TRPTRS_ther
                         - R_MALTabc_ther - R_NDPK5_ther - R_ASPTRS_ther
                         - R_GK1_ther - R_PROTRS_ther - R_GLNTAL_ther
                         - R_HMPK1_ther - R_CTPS1_1_ther - R_FTHFL_ther
                         - R_ILEabc_ther  = 0
 M_bilea_c_ther:         R_BILEt4_ther  = 0
 M_bilea_e:              - R_BILEt4_ther  = 0
 M_btn_c_bulg:           R_BTNt2i_bulg - 0.00001 R_biomass_LBUL10_bulg  = 0
 M_btn_e:                - R_BTNt2i_bulg - R_EX_btn_e_  = 0
 M_butacp_c_bulg:        - R_MCMAT2_bulg + R_BTMAT1_bulg  = 0
 M_butacp_c_ther:        - R_MCMAT2_ther + R_BTMAT1_ther  = 0
 M_cas_e:                - R_CASPROTL_bulg - R_EX_cas  = 0
 M_caspep_c_bulg:        - 0.001 R_CASL_bulg + 0.001 R_CASabc_bulg  = 0
 M_caspep_c_ther:        - 0.001 R_CASL_ther + 0.001 R_CASabc_ther  = 0
 M_caspep_e:             R_CASPROTL_bulg - R_EX_caspep_e_ - 0.001 R_CASabc_bulg
                         - 0.001 R_CASabc_ther  = 0
 M_cbasp_c_bulg:         R_DHORTS_bulg + R_ASPCT_bulg  = 0
 M_cbasp_c_ther:         R_ASPCT_ther + R_DHORTS_ther  = 0
 M_cbp_c_bulg:           - R_ASPCT_bulg + R_CBPS_bulg  = 0
 M_cbp_c_ther:           - R_OCBT_ther - R_ASPCT_ther + R_CBPS_ther  = 0
 M_cdp_c_bulg:           R_CYTK1_bulg - R_NDPK3_bulg + R_PPNCL_bulg  = 0
 M_cdp_c_ther:           - R_RNDR3_ther + R_CYTK1_ther + R_PPNCL_ther
                         - R_NDPK3_ther  = 0
 M_cdpdag_LLA_c_ther:    0.01 R_DASYN_LLA_ther - 0.01 R_PGSA_LLA_ther  = 0
 M_cdpdag_LPL_c_bulg:    0.01 R_DASYN_LPL_bulg - 0.01 R_PGSA_LPL_bulg  = 0
 M_cdpglyc_c_bulg:       R_G3PCT_bulg - 25 R_GTAglcS_bulg - R_CGPT_bulg  = 0
 M_cellb_e:              - R_CELBpts_bulg - R_CELBpts_ther  = 0
 M_chol_c_bulg:          R_GPDDA1_bulg + R_CHLabc_bulg  = 0
 M_chol_e:               - R_CHLabc_bulg  = 0
 M_chor_c_bulg:          - R_CHORM_bulg  = 0
 M_chor_c_ther:          - R_ANS1_ther + R_CHORS_ther - R_CHORM_ther
                         - R_ADCS_ther  = 0
 M_cit_c_ther:           - R_ACONT_ther + R_CITt6_ther + R_CS_ther  = 0
 M_cit_e:                - R_CITt6_ther - R_EX_cit_e_  = 0
 M_citr_L_c_ther:        R_OCBT_ther - R_ARGSS_ther  = 0
 M_clpn_LLA_c_ther:      - 0.000241 R_biomass_STR_ther + 0.01 R_CLPNS_LLA_ther
                          = 0
 M_clpn_LPL_c_bulg:      0.01 R_CLPNS_LPL_bulg - 0.000018 R_biomass_LBUL10_bulg
                          = 0
 M_cmp_c_bulg:           - R_CYTK1_bulg + 25 R_GTAglcS_bulg + R_CYTDK3_bulg
                         + R_CYTDK1_bulg + R_CYTDK2_bulg + R_CGPT_bulg
                         + R_PGSA_LPL_bulg  = 0
 M_cmp_c_ther:           - R_CYTK1_ther - R_NTD4_ther + R_CYTDK3_ther
                         + R_CYTDK1_ther + R_CYTDK2_ther + R_PGSA_LLA_ther  = 0
 M_co2_c_bulg:           6 R_KAS5_bulg + R_MCMAT2_bulg + R_MCMAT5_bulg
                         + R_OMPDC_bulg + R_DAPDC_bulg + R_THZPSN2_bulg
                         + 6 R_KAS2_bulg + 5 R_KAS1_bulg + R_MCMAT8_bulg
                         + R_PGDH_bulg - R_HCO3E_bulg + 5 R_KAS3_bulg
                         + 6 R_KAS9_bulg + 5 R_KAS4_bulg + R_MCMAT7_bulg
                         + 6 R_KAS10_bulg + R_LYSDC_bulg + R_PPCDC_bulg
                         + 7 R_KAS7_bulg + R_DPMVD_bulg + 8 R_KAS13_bulg
                         + 7 R_KAS8_bulg + R_kaasIII_bulg + R_MCMAT3_bulg
                         + R_MCMAT4_bulg + R_CO2t_bulg - R_AIRC_bulg
                         - R_PPC_bulg + R_MCMAT6_bulg + 6 R_KAS12_bulg
                         + 6 R_KAS11_bulg + R_ARGDC_bulg + R_PYROX_bulg
                         + 6 R_KAS6_bulg + 8 R_KAS17_bulg  = 0
 M_co2_c_ther:           R_ICDHx_ther + R_MCMAT2_ther + R_DAPDC_ther
                         + R_IGPS_ther + R_ACLDC_ther + R_DPMVD_ther
                         + R_OMCDC_ther + R_ACLS_ther + R_MCMAT6_ther
                         + R_MCMAT8_ther + R_OMPDC_ther + R_MCMAT3_ther
                         + R_CO2t_ther + R_PPNDH_ther + R_PPCDC_ther
                         + R_MCMAT4_ther + R_ICDHy_ther - R_PPC_ther
                         + R_PDH_ther - R_AIRC_ther + R_ACHBS_ther
                         + R_UREA2_ther + R_MCMAT5_ther - R_HCO3E_ther
                         + R_PPND_ther + R_kaasIII_ther + R_MCMAT7_ther
                         + R_PYROX_ther + 2 R_CPPPGO_ther  = 0
 M_co2_e:                - R_CO2t_ther - R_CO2t_bulg - R_EX_co2_e_  = 0
 M_coa_c_bulg:           7 R_KAS5_bulg + R_G1PACT_bulg - R_ACPS2_bulg
                         - R_HMGCOAR_bulg - R_ACALD_bulg + R_HSAT_bulg
                         + R_ACACT1r_bulg + 7 R_KAS2_bulg + 6 R_KAS1_bulg
                         - 0.0002 R_biomass_LBUL10_bulg + 6 R_KAS3_bulg
                         + R_PTAr_bulg + 7 R_KAS9_bulg + 6 R_KAS4_bulg
                         + 7 R_KAS10_bulg + R_DPCOAK_bulg + 8 R_KAS7_bulg
                         + 9 R_KAS13_bulg + 8 R_KAS8_bulg + R_kaasIII_bulg
                         + R_MACPMT_bulg + 7 R_KAS12_bulg - R_HMGCOAS_bulg
                         + 7 R_KAS11_bulg + 7 R_KAS6_bulg + 9 R_KAS17_bulg  = 0
 M_coa_c_ther:           R_THDPS_ther - R_HMGCOAS_ther + R_MACPMT_ther
                         + R_G1PACT_ther + R_PTAr_ther - R_PFL_ther
                         + R_SERAT_ther + R_ACACT1_ther
                         - 0.00018 R_biomass_STR_ther - R_ACPS2_ther
                         + R_HSAT_ther + R_CS_ther - R_PDH_ther
                         - R_HMGCOAR_ther + R_APAT_ther + R_DPCOAK_ther
                         + R_IPPS_ther + R_kaasIII_ther  = 0
 M_cpocdacp_c_bulg:      - 0.25 R_AGAT_LPL_bulg + R_CFAS180_bulg
                         - 0.25 R_GAT1_LPL_bulg  = 0
 M_cpocdacp_c_ther:      R_CFAS180_ther - 0.13 R_AGAT_LLA_ther
                         - 0.13 R_GAT1_LLA_ther  = 0
 M_cpppg3_c_ther:        - R_CPPPGO_ther  = 0
 M_csn_c_bulg:           - R_CSND_bulg  = 0
 M_csn_c_ther:           - R_CSND_ther  = 0
 M_ctp_c_bulg:           - 0.2 R_RNAS_LPL_bulg + R_CTPS1_1_bulg - R_G3PCT_bulg
                         + R_NDPK3_bulg + R_CTPS2_bulg - R_RNTR3_bulg
                         - R_DASYN_LPL_bulg - R_PPNCL_bulg  = 0
 M_ctp_c_ther:           R_CTPS2_ther - R_PPNCL_ther - R_RNTR3_ther
                         - R_DASYN_LLA_ther + R_NDPK3_ther
                         - 0.2 R_RNAS_LLA_ther + R_CTPS1_1_ther  = 0
 M_cys_L_c_bulg:         M_cys_L_c_bulg_src - R_CYSTRS_bulg + R_CYSTGL_bulg
                         + R_CYSt6_bulg - R_MOADCST_bulg
                         + 0.00227831258377267 R_CASL_bulg - R_PPNCL_bulg
                         + 0.014 R_PEPTL_bulg + R_CYSS_bulg  = 0
 M_cys_L_c_ther:         - R_PPNCL_ther - R_CYSTRS_ther + R_CYSt6_ther
                         - R_ASPTA4_ther + 0.00227831258377267 R_CASL_ther
                         - R_MOADCST_ther + R_CYSS_ther + R_AHSERL4_ther
                         + R_CYSTGL_ther + 0.014 R_PEPTL_ther - R_SHSL5_ther
                          = 0
 M_cys_L_e:              - R_EX_cys_L_e_ - R_CYSt6_ther - R_CYSt6_bulg  = 0
 M_cysth_L_c_bulg:       - R_CYSTGL_bulg - R_CYSTL_bulg  = 0
 M_cysth_L_c_ther:       - R_CYSTL_ther - R_CYSTGL_ther + R_SHSL5_ther  = 0
 M_cystrna_c_bulg:       R_CYSTRS_bulg - 0.011 R_PROTS_LPL_v60_bulg  = 0
 M_cystrna_c_ther:       - 0.006 R_PROTS_STR_ther + R_CYSTRS_ther  = 0
 M_cytd_c_bulg:          - R_CYTD_bulg - R_CYTDK3_bulg - R_CYTDK1_bulg
                         - R_CYTDK2_bulg  = 0
 M_cytd_c_ther:          R_NTD4_ther - R_CYTDK3_ther - R_CYTDK1_ther
                         - R_CYTDK2_ther - R_CYTD_ther  = 0
 M_d12dg_LLA_c_ther:     0.01 R_DAGGT_LLA_ther - 0.000172 R_biomass_STR_ther
                         - 0.01 R_LTAS_LLA_ther + 0.01 R_UGT2_LLA_ther  = 0
 M_dad_2_c_bulg:         - R_DADNK_bulg - R_DADA_bulg - R_PUNP2_bulg  = 0
 M_dad_2_c_ther:         - R_PUNP2_ther + R_NTD6_ther - R_DADA_ther  = 0
 M_dadp_c_bulg:          R_GK2_bulg - R_NDPK8_bulg  = 0
 M_dadp_c_ther:          R_GK2_ther - R_NDPK8_ther + R_RNDR1_ther  = 0
 M_dag_LPL_c_bulg:       - 0.01 R_DAGGT_LPL_bulg  = 0
 M_damp_c_bulg:          R_DADNK_bulg  = 0
 M_damp_c_ther:          - R_NTD6_ther  = 0
 M_datp_c_bulg:          - 0.28 R_DNAS_LPL_bulg - R_GK2_bulg + R_NDPK8_bulg
                         + R_RNTR1_bulg  = 0
 M_datp_c_ther:          - 0.305 R_DNAS_STU_ther - R_GK2_ther + R_RNTR1_ther
                         + R_NDPK8_ther  = 0
 M_dcamp_c_bulg:         R_ADSS_bulg - R_ADSL1_bulg  = 0
 M_dcamp_c_ther:         - R_ADSL1r_ther + R_ADSS_ther  = 0
 M_dcdp_c_bulg:          - R_NDPK7_bulg + R_CYTK2_bulg  = 0
 M_dcdp_c_ther:          R_RNDR3_ther + R_CYTK2_ther - R_NDPK7_ther  = 0
 M_dcmp_c_bulg:          - R_CYTK2_bulg  = 0
 M_dcmp_c_ther:          - R_CYTK2_ther - R_NTD3_ther  = 0
 M_dctp_c_bulg:          - 0.22 R_DNAS_LPL_bulg + R_RNTR3_bulg + R_NDPK7_bulg
                          = 0
 M_dctp_c_ther:          - 0.195 R_DNAS_STU_ther + R_RNTR3_ther + R_NDPK7_ther
                          = 0
 M_dcyt_c_bulg:          - R_DCYTD_bulg  = 0
 M_dcyt_c_ther:          - R_DCYTD_ther + R_NTD3_ther  = 0
 M_ddeacp_c_bulg:        R_DDMAT5_bulg - R_MCMAT6_bulg  = 0
 M_ddeacp_c_ther:        R_DDMAT5_ther - R_MCMAT6_ther  = 0
 M_decacp_c_bulg:        - R_MCMAT5_bulg + R_DEMAT4_bulg  = 0
 M_decacp_c_ther:        R_DEMAT4_ther - R_MCMAT5_ther  = 0
 M_dgdag_LPL_c_bulg:     0.01 R_DAGGT_LPL_bulg  = 0
 M_dgdp_c_bulg:          - R_NDPK5_bulg  = 0
 M_dgdp_c_ther:          R_RNDR2_ther - R_NDPK5_ther  = 0
 M_dgmp_c_bulg:          R_DGNSK_bulg  = 0
 M_dgmp_c_ther:          - R_NTD8_ther  = 0
 M_dgsn_c_bulg:          - R_DGNSK_bulg - R_PUNP4_bulg  = 0
 M_dgsn_c_ther:          R_NTD8_ther - R_PUNP4_ther  = 0
 M_dgtp_c_bulg:          - 0.22 R_DNAS_LPL_bulg + R_RNTR2_bulg + R_NDPK5_bulg
                          = 0
 M_dgtp_c_ther:          - 0.195 R_DNAS_STU_ther + R_NDPK5_ther + R_RNTR2_ther
                          = 0
 M_dha_c_bulg:           - R_DHAPT_bulg  = 0
 M_dha_c_ther:           R_DHAt_ther - R_DHAPT_ther  = 0
 M_dha_e:                - R_DHAt_ther - R_EX_dha_e_  = 0
 M_dhap_c_bulg:          R_DHAPT_bulg - R_TPI_bulg + R_G3PD1_bulg + R_FBA_bulg
                          = 0
 M_dhap_c_ther:          - R_TPI_ther + R_G3PD1_ther + R_DHAPT_ther
                         + R_FBA_ther  = 0
 M_dhf_c_bulg:           R_TMDS_bulg + R_DHFS_bulg - R_DHFR_bulg
                         - R_DHFOR2_bulg  = 0
 M_dhf_c_ther:           - R_DHFOR2_ther + R_TMDS_ther + R_DHFS_ther
                         - R_DHFR_ther  = 0
 M_dhnpt_c_bulg:         - R_DHNPA_bulg + R_DNMPPA_bulg  = 0
 M_dhnpt_c_ther:         R_DNMPPA_ther - R_DHNPA_ther  = 0
 M_dhor_S_c_bulg:        - R_DHORD6_bulg - R_DHORTS_bulg  = 0
 M_dhor_S_c_ther:        - R_DHORD6_ther - R_DHORTS_ther - R_DHORD1_ther  = 0
 M_dhpmp_c_bulg:         R_DNTPPA_bulg - R_DNMPPA_bulg  = 0
 M_dhpmp_c_ther:         R_DNTPPA_ther - R_DNMPPA_ther  = 0
 M_dhpt_c_bulg:          - R_DHFS_bulg + R_DHPS3_bulg  = 0
 M_dhpt_c_ther:          - R_DHFS_ther + R_DHPS3_ther  = 0
 M_diact_c_ther:         - R_DIACTt_ther  = 0
 M_diact_e:              R_DIACTt_ther - R_EX_diact_e_  = 0
 M_dimp_c_bulg:          R_DINSK_bulg  = 0
 M_din_c_bulg:           R_DADA_bulg - R_DINSK_bulg - R_PUNP6_bulg  = 0
 M_din_c_ther:           - R_PUNP6_ther + R_DADA_ther  = 0
 M_dmpp_c_bulg:          R_IPDDI_bulg - R_DMATT_bulg  = 0
 M_dmpp_c_ther:          R_IPDDI_ther - R_DMATT_ther  = 0
 M_dnad_c_bulg:          R_NNAT_bulg - R_NADS1_bulg  = 0
 M_dnad_c_ther:          - R_NADS1_ther + R_NNAT_ther  = 0
 M_dpcoa_c_bulg:         - R_DPCOAK_bulg + R_PTPAT_bulg  = 0
 M_dpcoa_c_ther:         R_PTPAT_ther - R_DPCOAK_ther  = 0
 M_drib_c_bulg:          - R_DRBK_bulg  = 0
 M_dtdp6dm_c_bulg:       - R_TDPDRR_bulg - R_CPSS_LPL_bulg  = 0
 M_dtdp6dm_c_ther:       - R_TDPDRR_ther - R_CPSS_LLA_ther  = 0
 M_dtdp_c_bulg:          R_CPSS_LPL_bulg + R_DTMPK_bulg - R_NDPK4_bulg  = 0
 M_dtdp_c_ther:          - R_NDPK4_ther + R_DTMPK_ther + R_CPSS_LLA_ther  = 0
 M_dtdpddg_c_bulg:       - R_TDPDRE_bulg + R_TDPGDH_bulg  = 0
 M_dtdpddg_c_ther:       - R_TDPDRE_ther + R_TDPGDH_ther  = 0
 M_dtdpddm_c_bulg:       R_TDPDRR_bulg + R_TDPDRE_bulg  = 0
 M_dtdpddm_c_ther:       R_TDPDRR_ther + R_TDPDRE_ther  = 0
 M_dtdpglc_c_bulg:       R_G1PTMT_bulg - R_TDPGDH_bulg  = 0
 M_dtdpglc_c_ther:       R_G1PTMT_ther - R_TDPGDH_ther  = 0
 M_dtmp_c_bulg:          R_TMDS_bulg + R_TMDK1_bulg - R_DTMPK_bulg
                         + R_TMDK2_bulg  = 0
 M_dtmp_c_ther:          R_TMDS_ther + R_TMDK2_ther + R_TMDK1_ther
                         - R_NTD5_ther - R_DTMPK_ther  = 0
 M_dttp_c_bulg:          - 0.28 R_DNAS_LPL_bulg - R_G1PTMT_bulg + R_NDPK4_bulg
                          = 0
 M_dttp_c_ther:          - 0.305 R_DNAS_STU_ther - R_G1PTMT_ther + R_NDPK4_ther
                          = 0
 M_dudp_c_bulg:          - R_NDPK6_bulg + R_URIDK2_bulg  = 0
 M_dudp_c_ther:          R_URIDK2_ther + R_RNDR4_ther - R_NDPK6_ther  = 0
 M_dump_c_bulg:          - R_TMDS_bulg + R_DURIK1_bulg + R_DUTPDP_bulg
                         - R_URIDK2_bulg  = 0
 M_dump_c_ther:          - R_URIDK2_ther - R_TMDS_ther - R_NTD1_ther
                         + R_DUTPDP_ther + R_DURIK1_ther  = 0
 M_duri_c_bulg:          - R_DURIK1_bulg - R_DURIPP_bulg + R_DCYTD_bulg  = 0
 M_duri_c_ther:          R_NTD1_ther - R_DURIPP_ther - R_DURIK1_ther
                         + R_DCYTD_ther  = 0
 M_dutp_c_bulg:          R_NDPK6_bulg - R_DUTPDP_bulg + R_RNTR4_bulg  = 0
 M_dutp_c_ther:          R_RNTR4_ther - R_DUTPDP_ther + R_NDPK6_ther  = 0
 M_dxyl5p_c_bulg:        - R_THZPSN2_bulg  = 0
 M_e4p_c_ther:           - R_TKT2_ther - R_DAHPS_ther  = 0
 M_etha_c_bulg:          R_GPDDA2_bulg  = 0
 M_etoh_c_bulg:          - R_ALCD2x_bulg + R_ETOHt1_bulg  = 0
 M_etoh_c_ther:          - R_ALCD2x_ther + R_ETOHt1_ther  = 0
 M_etoh_e:               - R_EX_etoh_e_ - R_ETOHt1_ther - R_ETOHt1_bulg  = 0
 M_f1p_c_bulg:           R_FRUpts_bulg - R_F1PK_bulg  = 0
 M_f1p_c_ther:           - R_F1PK_ther  = 0
 M_f6p_c_bulg:           R_G6PDA_bulg - R_GF6PTA_bulg - R_PFK_bulg + R_PGI_bulg
                         + R_MAN6PI_bulg + R_FRUK_bulg  = 0
 M_f6p_c_ther:           R_G6PDA_ther + R_TKT2_ther - R_GF6PTA_ther
                         + R_FRUK_ther - R_PFK_ther + R_PGI_ther
                         + R_MAN6PI_ther  = 0
 M_fa10_c_bulg:          R_KAS10_bulg  = 0
 M_fa11_c_bulg:          R_KAS11_bulg  = 0
 M_fa12_c_bulg:          R_KAS12_bulg  = 0
 M_fa1_c_bulg:           R_KAS1_bulg  = 0
 M_fa3_c_bulg:           R_KAS3_bulg  = 0
 M_fa4_c_bulg:           R_KAS4_bulg  = 0
 M_fa5_c_bulg:           R_KAS5_bulg  = 0
 M_fa6_c_bulg:           R_KAS6_bulg  = 0
 M_fa9_c_bulg:           R_KAS9_bulg  = 0
 M_fad_c_bulg:           R_FMNAT_bulg  = 0
 M_fad_c_ther:           R_FMNAT_ther - 0.00001 R_biomass_STR_ther  = 0
 M_fdp_c_bulg:           R_PFK_bulg - R_FBA_bulg + R_F1PK_bulg  = 0
 M_fdp_c_ther:           R_F1PK_ther - R_FBA_ther + R_PFK_ther  = 0
 M_fe3_c_ther:           R_FE3t4_ther  = 0
 M_fe3_e:                - R_FE3t4_ther  = 0
 M_ferrich_c_ther:       R_FERRICHabc_ther  = 0
 M_ferrich_e:            - R_FERRICHabc_ther  = 0
 M_fgam_c_bulg:          - R_PRFGS_bulg + R_GARFT_bulg  = 0
 M_fgam_c_ther:          - R_PRFGS_ther + R_GARFT_ther  = 0
 M_fmet_c_bulg:          - R_FMETDF_bulg  = 0
 M_fmet_c_ther:          - R_FMETDF_ther  = 0
 M_fmn_c_bulg:           - R_FMNAT_bulg + R_RBFK_bulg  = 0
 M_fmn_c_ther:           - R_FMNAT_ther + R_RBFK_ther  = 0
 M_fol_c_bulg:           R_FOLt_bulg + R_DHFOR2_bulg  = 0
 M_fol_c_ther:           R_DHFOR2_ther + R_FOLt_ther  = 0
 M_fol_e:                - R_FOLt_bulg - R_FOLt_ther - R_EX_fol_e_  = 0
 M_for_c_bulg:           R_GTPCI_bulg + R_FMETDF_bulg + R_FORt_bulg
                         - R_FTHFL_bulg  = 0
 M_for_c_ther:           R_PFL_ther - R_OXFOt4_ther + R_FORt2_ther
                         + R_GTPCI_ther + R_FMETDF_ther - R_FTHFL_ther  = 0
 M_for_e:                - R_EX_for_e_ - R_FORt_bulg + R_OXFOt4_ther
                         - R_FORt2_ther  = 0
 M_fpram_c_bulg:         R_PRFGS_bulg - R_PRAIS_bulg  = 0
 M_fpram_c_ther:         R_PRFGS_ther - R_PRAIS_ther  = 0
 M_fprica_c_bulg:        R_AICART_bulg + R_IMPC_bulg  = 0
 M_fprica_c_ther:        R_AICART_ther + R_IMPC_ther  = 0
 M_frdp_c_bulg:          - R_UDCPDPS_bulg + R_GRTT_bulg  = 0
 M_frdp_c_ther:          - R_UDCPDPS_ther + R_GRTT_ther  = 0
 M_fru_c_bulg:           - R_FRUK_bulg  = 0
 M_fru_c_ther:           - R_FRUK_ther  = 0
 M_fru_e:                - R_EX_fru_e_ - R_FRUpts_bulg  = 0
 M_fum_c_bulg:           - R_FRDx_bulg + R_ADSL2_bulg + R_ASPT_bulg
                         + R_SUCD3BC_bulg + R_ADSL1_bulg - R_FUM_bulg  = 0
 M_fum_c_ther:           - R_FRDx_ther + R_ADSL1r_ther + R_ARGSL_ther
                         + R_ADSL2r_ther  = 0
 M_g1p_B_c_ther:         R_PGMT_B_ther  = 0
 M_g1p_c_bulg:           - R_G1PTMT_bulg - R_PGMT_bulg - R_GALU_bulg  = 0
 M_g1p_c_ther:           R_GLCP_ther - R_G1PTMT_ther - R_GALU_ther
                         - R_PGMT_ther  = 0
 M_g3p_c_bulg:           R_TPI_bulg - R_GAPD_bulg + R_FBA_bulg + R_PKL_bulg
                          = 0
 M_g3p_c_ther:           R_TPI_ther + R_TKT2_ther + R_DRPA_ther + R_TRPS1_ther
                         + R_TKT1_ther - R_GAPD_ther + R_FBA_ther  = 0
 M_g3pc_c_bulg:          - R_GPDDA1_bulg  = 0
 M_g3pe_c_bulg:          - R_GPDDA2_bulg  = 0
 M_g3pg_c_bulg:          - R_GPDDA4_bulg  = 0
 M_g3pi_c_bulg:          - R_GPDDA5_bulg  = 0
 M_g3ps_c_bulg:          - R_GPDDA3_bulg  = 0
 M_g6p_B_c_bulg:         R_G6PI_bulg  = 0
 M_g6p_B_c_ther:         R_G6PI_ther - R_PGMT_B_ther  = 0
 M_g6p_c_bulg:           R_PGGH_bulg + R_HEX1_bulg + R_PGMT_bulg
                         + R_GLCpts_bulg - R_PGI_bulg - R_G6PI_bulg
                         - R_G6PDHy_bulg  = 0
 M_g6p_c_ther:           R_PGGH_ther - R_G6PI_ther + R_TRE6PH_ther
                         + R_HEX1_ther + R_PGMT_ther - R_PGI_ther
                         + R_GLCpts_ther  = 0
 M_gal1p_c_ther:         - R_GALT_ther + R_GALK2_ther  = 0
 M_gal_c_bulg:           R_LACZ_bulg + R_GALt6_bulg  = 0
 M_gal_c_ther:           - R_LCTSGALex_ther + R_LACZ_ther - R_GALM_ther  = 0
 M_gal_e:                - R_GALt6_bulg + R_LCTSGALex_ther - R_EX_gal_e_  = 0
 M_galt1p_c_bulg:        R_GALTpts_bulg  = 0
 M_galt1p_c_ther:        R_GALTpts_ther  = 0
 M_galt_e:               - R_GALTpts_bulg - R_GALTpts_ther  = 0
 M_gam1p_c_bulg:         - R_G1PACT_bulg - R_PGAMT_bulg  = 0
 M_gam1p_c_ther:         - R_G1PACT_ther - R_PGAMT_ther  = 0
 M_gam6p_c_bulg:         - R_G6PDA_bulg + R_GF6PTA_bulg + R_AGDC_bulg
                         + R_PGAMT_bulg  = 0
 M_gam6p_c_ther:         - R_G6PDA_ther + R_GF6PTA_ther + R_PGAMT_ther
                         + R_AGDC_ther  = 0
 M_gar_c_bulg:           - R_GARFT_bulg + R_PRAGS_bulg  = 0
 M_gar_c_ther:           R_PRAGS_ther - R_GARFT_ther  = 0
 M_gcald_c_bulg:         R_DHNPA_bulg  = 0
 M_gcald_c_ther:         R_GCALDt_ther + R_AHMMPS_ther + R_DHNPA_ther  = 0
 M_gcald_e:              - R_GCALDt_ther - R_EX_gcald_e_  = 0
 M_gdp_c_bulg:           R_GK2_bulg - R_NDPK1_bulg + R_GK1_bulg + R_ADSS_bulg
                         + R_CYTDK2_bulg + R_URIK2_bulg
                         + 2 R_PROTS_LPL_v60_bulg + R_TMDK2_bulg  = 0
 M_gdp_c_ther:           2 R_PROTS_STR_ther + R_GK2_ther + R_TMDK2_ther
                         + R_CYTDK2_ther + R_URIK2_ther - R_RNDR2_ther
                         - R_NDPK1_ther + R_ADSS_ther + R_GK1_ther  = 0
 M_gdptp_c_bulg:         R_GTPDPK_bulg  = 0
 M_gdptp_c_ther:         R_GTPDPK_ther  = 0
 M_glc_D_c_bulg:         R_PGGH_bulg - R_HEX1_bulg + R_LACZ_bulg
                         + 2 R_MALT_bulg  = 0
 M_glc_D_c_ther:         R_PGGH_ther + R_AMALT4_ther + R_TRE6PH_ther
                         + R_AMALT2_ther + R_AMALT3_ther + R_LACZ_ther
                         - R_HEX1_ther + R_AMALT1_ther  = 0
 M_glc_D_e:              - R_GLCpts_bulg - R_EX_glc_e_ - R_GLCpts_ther  = 0
 M_gln_L_c_bulg:         - R_PRFGS_bulg - R_GF6PTA_bulg + R_GLNS_bulg
                         + R_GLNabc_bulg - R_GLNTAL_bulg - R_GMPS2_bulg
                         - R_CTPS2_bulg - R_ASNS1_bulg + M_gln_L_c_bulg_src
                         - R_CBPS_bulg - R_GLUPRT_bulg + R_gln_export_bulg  = 0
 M_gln_L_c_ther:         - R_ANS1_ther - R_CTPS2_ther - R_PRFGS_ther
                         - R_GF6PTA_ther + R_GLNabc_ther + R_GLNS_ther
                         - R_GMPS2_ther - R_ADCS_ther - R_CBPS_ther
                         - R_GLUPRT_ther - R_GLNTAL_ther  = 0
 M_gln_L_e:              - R_GLNabc_bulg - R_GLNabc_ther - R_EX_gln_L_e_
                         - R_gln_export_bulg  = 0
 M_glntrna_c_bulg:       R_GLNTAL_bulg - 0.083 R_PROTS_LPL_v60_bulg  = 0
 M_glntrna_c_ther:       - 0.093 R_PROTS_STR_ther + R_GLNTAL_ther  = 0
 M_glu5p_c_ther:         R_GLU5K_ther - R_G5SD_ther  = 0
 M_glu5sa_c_ther:        R_G5SD_ther - R_G5SADs_ther  = 0
 M_glu_D_c_bulg:         - R_GLUR_bulg - R_UAMAGS_bulg  = 0
 M_glu_D_c_ther:         - R_UAMAGS_ther - R_GLUR_ther  = 0
 M_glu_L_c_bulg:         R_PRFGS_bulg + R_GF6PTA_bulg + R_GLUR_bulg
                         - R_THFGLUS_bulg - R_GLNS_bulg - R_DHFS_bulg
                         + R_GLUt6_bulg + M_glu_L_c_bulg_src + R_VALTA_bulg
                         + R_AKGTA3_bulg + R_GLNTAL_bulg + R_PHETA1_bulg
                         + R_LEUTA_bulg + R_GMPS2_bulg + R_TYRTA_bulg
                         + R_CTPS2_bulg + 0.132294586048683 R_CASL_bulg
                         - R_GLUTRS_bulg + 0.224 R_PEPTL_bulg + R_ILETA_bulg
                         + R_ASNS1_bulg + R_CBPS_bulg + R_GLUPRT_bulg  = 0
 M_glu_L_c_ther:         - R_GLUTRS_ther + R_ADPTA_ther + R_SDPTA_ther
                         + R_ANS1_ther + R_ASPTA5_ther - R_ORNTAC_ther
                         - R_PSERT_ther + R_CTPS2_ther - R_GLUDy_ther
                         + R_ASPTA6_ther + R_PRFGS_ther + R_VALTA_ther
                         + R_GF6PTA_ther - R_GLNS_ther - R_GLU5K_ther
                         + R_ILETA_ther - R_DHFS_ther + R_ASPTA1_ther
                         - R_THFGLUS_ther + R_GMPS2_ther + R_ASPTA4_ther
                         + R_GLUt6_ther + R_ADCS_ther
                         + 0.132294586048683 R_CASL_ther + R_GLUR_ther
                         + R_GLUabc_ther + R_AKGTA3_ther + 0.224 R_PEPTL_ther
                         + R_CBPS_ther + R_GLUPRT_ther + R_ALATA_L_ther
                         + R_ACOTA_ther + R_LEUTA_ther + R_GLNTAL_ther
                         + R_ASPTA2_ther  = 0
 M_glu_L_e:              - R_GLUt6_bulg - R_EX_glu_L_e_ - R_GLUt6_ther
                         - R_GLUabc_ther  = 0
 M_glutrna_c_bulg:       - R_GLNTAL_bulg + R_GLUTRS_bulg
                         - 0.023 R_PROTS_LPL_v60_bulg  = 0
 M_glutrna_c_ther:       R_GLUTRS_ther - 0.025 R_PROTS_STR_ther - R_GLNTAL_ther
                          = 0
 M_gly_c_bulg:           - R_THZPSN2_bulg - R_PRAGS_bulg - R_GLYTRS_bulg
                         + R_GLYt6_bulg + M_gly_c_bulg_src
                         + 0.10708544031126 R_CASL_bulg + 0.017 R_PEPTL_bulg
                          = 0
 M_gly_c_ther:           - R_PRAGS_ther - R_GLYTRS_ther
                         + 0.10708544031126 R_CASL_ther + 0.017 R_PEPTL_ther
                         + R_THRA_ther + R_GLYt6_ther + R_GHMT_ther  = 0
 M_gly_e:                - R_EX_gly_e_ - R_GLYt6_bulg - R_GLYt6_ther  = 0
 M_glyald_c_bulg:        - R_ALCD19_D_bulg  = 0
 M_glyald_c_ther:        - R_ALCD19_ther  = 0
 M_glyb_c_bulg:          R_GLYBabc_bulg  = 0
 M_glyb_c_ther:          R_GLYBabc_ther  = 0
 M_glyb_e:               - R_GLYBabc_bulg - R_GLYBabc_ther  = 0
 M_glyc3p_c_bulg:        R_GPDDA1_bulg - R_G3PCT_bulg + R_GLYK_bulg
                         + R_GPDDA2_bulg + R_GPDDA3_bulg - R_G3PD1_bulg
                         + R_GPDDA5_bulg - R_PGSA_LPL_bulg - R_GAT1_LPL_bulg
                         + R_GPDDA4_bulg  = 0
 M_glyc3p_c_ther:        - R_G3PD1_ther - R_PGSA_LLA_ther - R_GAT1_LLA_ther
                          = 0
 M_glyc_R_c_bulg:        - R_GLYCK_bulg  = 0
 M_glyc_R_c_ther:        - R_GLYCK_ther  = 0
 M_glyc_c_bulg:          R_ALCD19_D_bulg - R_GLYK_bulg + R_ALCD19_L_bulg
                         + R_GLYCt1_bulg + R_CLPNS_LPL_bulg + R_GPDDA4_bulg
                          = 0
 M_glyc_c_ther:          R_ALCD19_ther + R_CLPNS_LLA_ther - R_GLYCt5_ther  = 0
 M_glyc_e:               - R_EX_glyc_e_ - R_GLYCt1_bulg + R_GLYCt5_ther  = 0
 M_glyclt_c_bulg:        R_PGLYCP_bulg  = 0
 M_glycogen_c_ther:      - R_GLCP_ther  = 0
 M_glytrna_c_bulg:       R_GLYTRS_bulg - 0.084 R_PROTS_LPL_v60_bulg  = 0
 M_glytrna_c_ther:       - 0.054 R_PROTS_STR_ther + R_GLYTRS_ther  = 0
 M_gmp_c_bulg:           - R_GK2_bulg + R_GUAPRT_bulg - R_GK1_bulg
                         + R_GMPS2_bulg - R_GMPR_bulg  = 0
 M_gmp_c_ther:           - R_NTD9_ther - R_GMPR_ther - R_GK2_ther
                         + R_GUAPRT_ther + R_GMPS2_ther - R_GK1_ther  = 0
 M_grdp_c_bulg:          R_DMATT_bulg - R_GRTT_bulg  = 0
 M_grdp_c_ther:          - R_GRTT_ther + R_DMATT_ther  = 0
 M_gsn_c_bulg:           - R_PUNP3_bulg - R_PNS2_bulg  = 0
 M_gsn_c_ther:           R_NTD9_ther - R_PUNP3_ther  = 0
 M_gthox_c_ther:         R_GTHRD_ther  = 0
 M_gthrd_c_ther:         - R_LGTHL_ther - 2 R_GTHRD_ther  = 0
 M_gtp_c_bulg:           - 0.32 R_RNAS_LPL_bulg - R_GTPCI_bulg + R_NDPK1_bulg
                         - R_RNTR2_bulg - R_ADSS_bulg - R_GTPDPK_bulg
                         - R_CYTDK2_bulg - R_URIK2_bulg
                         - 2 R_PROTS_LPL_v60_bulg - R_TMDK2_bulg  = 0
 M_gtp_c_ther:           - R_GTPDPK_ther - 2 R_PROTS_STR_ther - R_TMDK2_ther
                         - R_PZS_ther - R_CYTDK2_ther - 0.32 R_RNAS_LLA_ther
                         - R_URIK2_ther + R_NDPK1_ther - R_ADSS_ther
                         - R_GTPCI_ther - R_RNTR2_ther  = 0
 M_gua_c_bulg:           R_GUAt2_bulg + R_PUNP3_bulg + R_PNS2_bulg
                         - R_GUAPRT_bulg + R_PUNP4_bulg  = 0
 M_gua_c_ther:           - R_GUAPRT_ther + R_PUNP3_ther + R_PUNP4_ther  = 0
 M_gua_e:                - R_GUAt2_bulg - R_EX_gua_e_  = 0
 M_h2o2_c_bulg:          R_PYROX_bulg + R_LOXD_bulg  = 0
 M_h2o2_c_ther:          R_SOD_ther + R_DHORD1_ther + R_NOX1_ther
                         + R_PYROX_ther  = 0
 M_h2o_c_bulg:           - 0.4 R_RNAS_LPL_bulg - R_METAT_bulg
                         - 1.37 R_DNAS_LPL_bulg - R_MNabc_bulg - R_SSALx_bulg
                         + 5 R_KAS5_bulg - R_G6PDA_bulg + R_ENO_bulg
                         - R_PGGH_bulg - R_PRFGS_bulg - R_GTPCI_bulg
                         + R_H2Ot5_bulg - 2 R_PEPTabc_bulg - R_UDCPDP_bulg
                         - R_GLYBabc_bulg - R_GPDDA1_bulg - R_RHC_bulg
                         - R_AEPabc_bulg - R_PPA_bulg - R_CHLabc_bulg
                         - R_CYSTGL_bulg + R_HHDHL7_bulg - R_CYSTL_bulg
                         - R_THMabc_bulg - R_LACZ_bulg - R_PNS4_bulg
                         + 2 R_THZPSN2_bulg + R_RNTR2_bulg - R_THRS_bulg
                         + 5 R_KAS2_bulg - R_PNS2_bulg - R_ADPRDP_bulg
                         - R_GLNabc_bulg - R_FMETDF_bulg + 4 R_KAS1_bulg
                         + R_YUMPS_bulg - R_PIabc_bulg - R_CSND_bulg
                         - R_DNTPPA_bulg - R_CYTD_bulg - R_ADA_bulg
                         - R_DNMPPA_bulg + R_ATPS3r_bulg - R_GPDDA2_bulg
                         + R_HHYHL2_bulg - 27.2 R_biomass_LBUL10_bulg
                         - R_ZN2abc2_bulg - R_DADA_bulg - R_HCO3E_bulg
                         + R_TDPGDH_bulg + 4 R_KAS3_bulg - R_GPDDA3_bulg
                         - R_MOADCST_bulg - R_PYRNS1_bulg - R_MTAN_bulg
                         - R_GLNTAL_bulg - R_PGL_bulg + 5 R_KAS9_bulg
                         - R_PGLYCP_bulg - R_EPPP_bulg - R_GMPS2_bulg
                         + 4 R_KAS4_bulg + 5 R_KAS10_bulg + R_HDDHL5_bulg
                         - R_CTPS2_bulg - 5 R_CPSS_LPL_bulg - R_MALT_bulg
                         + 6 R_KAS7_bulg + R_RNTR3_bulg - R_DUTPDP_bulg
                         - R_PTRCabc_bulg - R_NNAMr_bulg + 7 R_KAS13_bulg
                         - R_CASL_bulg - R_SPMDabc_bulg + R_HODHL8_bulg
                         + 6 R_KAS8_bulg - 0.111 R_CASabc_bulg - R_MALTabc_bulg
                         - R_MTHFC_bulg + R_PKL_bulg - R_PEPTL_bulg
                         - R_DHORTS_bulg - R_ARGabc_bulg - R_IMPD_bulg
                         + R_HTDHL6_bulg + R_RNTR4_bulg - R_PNS3_bulg
                         - R_ATPM_bulg - R_PGPP_LPL_bulg
                         - 2.306 R_PROTS_LPL_v60_bulg - R_NAabcO_bulg
                         - R_AHCYSNS_bulg - R_PPC_bulg - R_GPDDA5_bulg
                         - R_AGDC_bulg + 5 R_KAS12_bulg - R_ASNN_bulg
                         + R_HMGCOAS_bulg - R_ASNS1_bulg + 5 R_KAS11_bulg
                         + R_HBUHL1_bulg - R_SSALy_bulg - R_CBPS_bulg
                         - R_PNS1_bulg + R_SUCBZS_bulg - R_GPDDA4_bulg
                         - R_DCYTD_bulg + R_HDEHL4_bulg - R_GLUPRT_bulg
                         + 5 R_KAS6_bulg - R_FUM_bulg + R_RNTR1_bulg
                         + 7 R_KAS17_bulg - R_IMPC_bulg + R_HOCHL3_bulg  = 0
 M_h2o_c_ther:           - R_NTD11_ther + R_RNDR3_ther - R_MALTTTRabc_ther
                         - R_PSP_L_ther - R_CYSTL_ther - R_PGGH_ther
                         - R_VALabc_ther - R_G6PDA_ther - R_METAT_ther
                         - R_PIabc_ther + R_IGPS_ther - R_MTHFC_ther
                         - R_MI1PP_ther + R_HHDHL7_ther - R_PGPP_LLA_ther
                         - R_THDPS_ther + 2 R_DHDPS_ther - R_NTD7_ther
                         - R_NTD9_ther - R_AGMT_ther - R_PAP_LLA_ther
                         - R_NAMNPP_ther - R_CTPS2_ther + R_HBUHL1_ther
                         - R_SPMDabc_ther - R_ACTPASE_ther - R_GLUDy_ther
                         + R_HDEHL4_ther - 1.37 R_DNAS_STU_ther - R_NTD10_ther
                         - R_NTD4_ther - R_NNAMr_ther + R_RNDR4_ther
                         + R_HDDHL5_ther - 2.306 R_PROTS_STR_ther
                         + R_ATPS3r_ther + R_HMGCOAS_ther - R_DNTPPA_ther
                         - R_ASPabc_ther + R_ENO_ther - R_MALTHXabc_ther
                         - R_MALTPTabc_ther - R_PRFGS_ther - 2 R_PEPTabc_ther
                         - R_ZN2abc2_ther - R_METabc_ther + R_RNTR3_ther
                         - R_GLNabc_ther - R_UDCPDP_ther + R_TRPS1_ther
                         - R_PTRCabc_ther - R_NTD8_ther - R_USHD_ther
                         - R_THMabc_ther - R_METSRabc_ther - R_DHORTS_ther
                         + R_DHAD3_ther - R_CSND_ther - R_ADPRDP_ther
                         + R_RNTR1_ther - R_NTD2_ther - R_SUC6PH_ther
                         - R_MNabc_ther - R_DNMPPA_ther - R_NTD1_ther
                         + 2 R_NOX2_ther - R_GLYBabc_ther - R_PPA_ther
                         + R_RNTR4_ther - R_AHCYSNS_ther - R_CYTD_ther
                         - R_TRE6PH_ther + R_HHYHL2_ther - R_SSALx_ther
                         - R_GMPS2_ther - R_ALAD_Lr_ther - R_FERRICHabc_ther
                         - R_DUTPDP_ther - R_METSabc_ther
                         - 18.15 R_biomass_STR_ther + R_HODHL8_ther
                         + R_G5SADs_ther - 0.4 R_RNAS_LLA_ther + R_PPNDH_ther
                         - R_AGDC_ther + R_DHQD_ther - R_ASNN_ther
                         - R_CASL_ther - R_MOADCST_ther - R_THRS_ther
                         + R_RNDR2_ther - R_GLUabc_ther + R_RNDR1_ther
                         - R_ECOAH2c_ther - R_LACZ_ther - R_MALTTRabc_ther
                         + R_YUMPS_ther - R_AHMMPS_ther + R_HTDHL6_ther
                         - R_PPC_ther - R_ATPM_ther - R_CS_ther - R_CYSTGL_ther
                         - R_GTPCI_ther - R_PROLBabc_ther - R_NTD6_ther
                         - R_LEUabc_ther - R_ARGabc_ther - R_PEPTL_ther
                         - R_APAT_ther + R_TDPGDH_ther - R_DCYTD_ther
                         - 0.111 R_CASabc_ther - R_UREA2_ther + R_DHAD1_ther
                         - R_DAHPS_ther - R_RHC_ther - R_IMPD_ther
                         - R_CBPS_ther - R_GLUPRT_ther - R_MTAN_ther
                         - R_ADPDS_ther - R_NTD5_ther - R_METDabc_2_ther
                         - R_IPPS_ther - R_IMPC_ther - R_HCO3E_ther
                         - R_FMETDF_ther - R_ADA_ther - R_IPPMIb_ther
                         + R_DHAD2_ther - R_MALTabc_ther + R_H2Ot5_ther
                         - R_URCN_ther + R_GHMT_ther - R_NTD3_ther
                         - R_SSALy_ther + R_IPPMIa_ther - 4 R_CPSS_LLA_ther
                         - R_GLNTAL_ther - R_DADA_ther + R_RNTR2_ther
                         + R_HOCHL3_ther + 2 R_CPPPGO_ther - R_ILEabc_ther  = 0
 M_h2o_e:                - R_H2Ot5_bulg - R_EX_h2o_e_ - R_H2Ot5_ther  = 0
 M_h2s_c_bulg:           - R_CYSS_bulg  = 0
 M_h2s_c_ther:           - R_METACH_ther - R_CYSS_ther  = 0
 M_h_c_bulg:             0.4 R_RNAS_LPL_bulg + R_DADNK_bulg + R_SHKK_bulg
                         + R_DRBK_bulg + 1.37 R_DNAS_LPL_bulg + R_MNabc_bulg
                         + R_ALATRS_bulg + 2 R_SSALx_bulg + R_GUAt2_bulg
                         - 16 R_KAS5_bulg + R_CYSTRS_bulg + R_ASPt6_bulg
                         - R_TRDR_bulg + 2 R_METFR_bulg + R_MELIBt2_bulg
                         + 2 R_PRFGS_bulg + R_G1PACT_bulg + R_ACCOAC_bulg
                         + R_ASNTRS_bulg + R_L_LACt2_bulg + R_GTPCI_bulg
                         + 2 R_UAAGLS1_bulg + R_RBK_bulg + 2 R_PEPTabc_bulg
                         - R_G1PTMT_bulg + R_UDCPDP_bulg + R_PGCD_bulg
                         + R_ARGTRS_bulg + R_GARFT_bulg - R_HOCR3_bulg
                         + R_URAt2_bulg + R_ACPS2_bulg - R_MCMAT2_bulg
                         + R_GLYBabc_bulg + R_THFGLUS_bulg + R_HSK_bulg
                         + R_UAGPT5_bulg + R_GPDDA1_bulg - R_NNAT_bulg
                         + R_GLNS_bulg - R_MCMAT5_bulg + R_AEPabc_bulg
                         + R_LEUTRS_bulg - R_HDMAT7_bulg - R_OMPDC_bulg
                         + R_DGNSK_bulg - R_FRDx_bulg + 2 R_HMGCOAR_bulg
                         + R_TYRt6_bulg + 2 R_CTPS1_1_bulg + R_CFAS180_bulg
                         + R_ACALD_bulg - R_G3PCT_bulg - R_DAPDC_bulg
                         + R_PPA_bulg + R_LYSt6_bulg + R_UGMDDS2_bulg
                         - R_ALCD19_D_bulg + R_MEVK_bulg + R_GAPD_bulg
                         + R_LCTSt6_bulg + R_GLYK_bulg + R_CHLabc_bulg
                         + 50 R_GTAglcS_bulg + R_PRASCS_bulg + R_NADK_bulg
                         + R_DHFS_bulg + R_SERTRS_bulg + R_HEX1_bulg
                         - R_UAPGR_bulg - R_HHYR2_bulg + R_PFK_bulg
                         + R_THMabc_bulg - R_ALCD19_L_bulg + 2 R_THZPSN2_bulg
                         + R_UDCPK_bulg + R_BTNt2i_bulg - 17 R_KAS2_bulg
                         + R_HSDy_bulg + 2 R_ADPRDP_bulg - R_OCDMAT8_bulg
                         + R_GLUt6_bulg - R_2HXICDH_bulg + R_GLNabc_bulg
                         - R_HEMAT2_bulg + R_MALt6_bulg + R_PRAGS_bulg
                         - 14 R_KAS1_bulg - R_FMNAT_bulg - R_MCMAT8_bulg
                         + R_ARGt_bulg - R_PYK_bulg + R_HXANt2_bulg
                         + R_HISt6_bulg + R_ASAD_bulg + R_PIabc_bulg
                         + R_ASNt6_bulg + R_TDPDRR_bulg + R_FOLt_bulg
                         - R_CSND_bulg + R_ASPTRS_bulg + R_NH4DIS_bulg
                         - R_CYTD_bulg + R_ACt6_bulg + R_VALTRS_bulg
                         - R_ADA_bulg + R_CYTDK3_bulg - R_DEMAT4_bulg
                         + 25 R_GRTAL_bulg + R_METTRS_bulg + R_DNMPPA_bulg
                         + R_GALt6_bulg + 2 R_ATPS3r_bulg + R_CYSt6_bulg
                         + R_GPDDA2_bulg + R_ASNS2_bulg + R_PROTRS_bulg
                         + 27.2 R_biomass_LBUL10_bulg + R_ZN2abc2_bulg
                         + R_UAMAS_bulg + R_PRPPS_bulg - R_DADA_bulg
                         + R_HCO3E_bulg + R_UAMAGS_bulg + R_HISTRS_bulg
                         - 14 R_KAS3_bulg + R_CYTDK1_bulg + R_GPDDA3_bulg
                         + R_PROt6_bulg - R_BTMAT1_bulg + R_PHEt6_bulg
                         + R_PHETRS_bulg - R_HBUR1_bulg + 2 R_ADSS_bulg
                         - R_TDMAT6_bulg + R_GLYTRS_bulg + R_GLNTAL_bulg
                         - R_UAGDP_bulg + R_GLYCK_bulg + R_PYDAMK_bulg
                         + R_GLYt6_bulg + R_PNTK_bulg + R_GTPDPK_bulg
                         + R_PGL_bulg - 16 R_KAS9_bulg + R_SERt6_bulg
                         + R_G3PD1_bulg + R_VALt6_bulg + R_XANt2_bulg
                         - R_DHFR_bulg + 2 R_EPPP_bulg + 2 R_GMPS2_bulg
                         - 14 R_KAS4_bulg - R_MCMAT7_bulg + R_PYDAMt_bulg
                         + 2 R_DAGGT_LPL_bulg - 16 R_KAS10_bulg + R_DPCOAK_bulg
                         + R_SUCCt6_bulg + R_PRAIS_bulg + R_DURIK1_bulg
                         - R_LYSDC_bulg + 2 R_CTPS2_bulg + 6 R_CPSS_LPL_bulg
                         - R_PPCDC_bulg + R_METS_bulg - R_PTPAT_bulg
                         + R_THRTRS_bulg - 19 R_KAS7_bulg + R_DUTPDP_bulg
                         - R_DASYN_LPL_bulg + R_DHORD6_bulg + R_PTRCabc_bulg
                         - 23 R_KAS13_bulg + R_NADS1_bulg + R_RIBFLVt2_bulg
                         + R_G6PDHy_bulg + R_GLUTRS_bulg + R_SPMDabc_bulg
                         - R_HTDR6_bulg + R_PPNCL_bulg - 20 R_KAS8_bulg
                         - R_kaasIII_bulg + 0.111 R_CASabc_bulg
                         + R_MALTabc_bulg + R_LDH_D_bulg - R_HODR8_bulg
                         + R_CYTDK2_bulg + R_MTHFC_bulg + R_THRt6_bulg
                         - R_DDMAT5_bulg + R_HCYSMT_bulg + R_LYSTRS_bulg
                         + R_TYRTRS_bulg + R_CYSS_bulg + R_TMDK1_bulg
                         + R_DHORTS_bulg + R_CGPT_bulg + R_ALA_Lt6_bulg
                         + R_ASPCT_bulg + R_UAAGDS_bulg - R_HHDR7_bulg
                         - R_MCMAT3_bulg - R_MCMAT4_bulg + R_ARGabc_bulg
                         + R_URIK2_bulg + R_IMPD_bulg - R_GALU_bulg
                         + R_TMDPK_bulg + R_TRPt6_bulg + R_ATPM_bulg
                         - 2 R_GMPR_bulg + R_URIK1_bulg + R_AIRC_bulg
                         + 2.306 R_PROTS_LPL_v60_bulg + R_NAabcO_bulg
                         + R_PPC_bulg - R_HDER4_bulg - R_MCMAT6_bulg
                         + 2 R_PYDXK_bulg + R_GPDDA5_bulg - 17 R_KAS12_bulg
                         + R_ILETRS_bulg - R_HMGCOAS_bulg + R_LEUt6_bulg
                         + R_UNAGAMAMT_bulg + R_ASNS1_bulg - 17 R_KAS11_bulg
                         + R_TRPTRS_bulg - R_OCMAT3_bulg + R_PGSA_LPL_bulg
                         + R_URIK3_bulg - R_ARGDC_bulg + R_LDH_L_bulg
                         - R_PYROX_bulg + R_ALCD2x_bulg - R_HDDR5_bulg
                         + R_METt6_bulg + 2 R_SSALy_bulg + R_FRUK_bulg
                         + 2 R_CBPS_bulg + R_TMDK2_bulg - R_NAPRT_bulg
                         + R_GPDDA4_bulg + R_PNTOt2_bulg - R_DCYTD_bulg
                         + R_F1PK_bulg + R_gln_export_bulg + R_ILEt6_bulg
                         - 17 R_KAS6_bulg + R_ALATA_Lr_bulg - 22 R_KAS17_bulg
                         + R_RBFK_bulg  = 0
 M_h_c_ther:             R_MALTTTRabc_ther + R_GLUTRS_ther + R_OCBT_ther
                         - R_ADPTA_ther - R_MCMAT2_ther + R_HISTRS_ther
                         + R_ANS1_ther - 2 R_SOD_ther - R_DAPDC_ther
                         + R_VALabc_ther + R_PIabc_ther + R_ACCOAC_ther
                         - R_GALT_ther - R_IGPS_ther + R_MTHFC_ther
                         + R_NADS1_ther - R_PTPAT_ther + R_GLYCK_ther
                         + R_PRAGS_ther + R_DHDPS_ther + R_ASNS2_ther
                         + R_RBK_ther + R_GTPDPK_ther + R_UAGPT1_ther
                         - R_ACLDC_ther - R_HHDR7_ther - R_HDER4_ther
                         + R_TDPDRR_ther + R_F1PK_ther + R_ARGSS_ther
                         - R_SHK3D_ther + R_TMDPK_ther + R_PYDAMK_ther
                         + 2 R_CTPS2_ther + R_VALt6_ther - R_OMCDC_ther
                         - R_UAPGR_ther + R_SPMDabc_ther + R_ACTPASE_ther
                         + R_SERt6_ther - 2 R_GMPR_ther + R_GLUDy_ther
                         + R_LCTSt6_ther - R_NMNS_ther + 2 R_PYDXK_ther
                         - R_FRDx_ther + R_ASPCT_ther + 1.37 R_DNAS_STU_ther
                         + R_MEVK_ther + R_URIK1_ther + R_HSDy_ther
                         + R_ILEt6_ther + R_VALTRS_ther
                         + 2.306 R_PROTS_STR_ther + 2 R_ATPS3r_ther
                         - R_HMGCOAS_ther - R_FMNAT_ther + R_CYTDK3_ther
                         - R_ACLS_ther - R_HOCR3_ther + 2 R_METFR_ther
                         + R_HSK_ther + R_URIK3_ther + R_DHORD6_ther
                         + R_ASPabc_ther - R_DDMAT5_ther + R_LDH_L_ther
                         + R_MALTHXabc_ther + R_MALTPTabc_ther - R_MCMAT6_ther
                         + 2 R_PRFGS_ther + R_PPNCL_ther + 2 R_PEPTabc_ther
                         + R_ZN2abc2_ther + R_CYSTRS_ther - R_DEMAT4_ther
                         + R_METabc_ther + R_PRAIS_ther + R_PNTK_ther
                         + R_GLNabc_ther + R_UDCPDP_ther + R_MALt6_ther
                         + R_GLNS_ther - R_HDDR5_ther + R_PTRCabc_ther
                         + R_G1PACT_ther + R_THRTRS_ther - R_MCMAT8_ther
                         - R_DPR_ther + R_G3PD1_ther + R_ADCL_ther
                         + 2 R_USHD_ther + R_CFAS180_ther + R_THMabc_ther
                         + R_METSRabc_ther + R_DHORTS_ther - R_OMPDC_ther
                         + R_CYSt6_ther - R_CSND_ther + 2 R_ADPRDP_ther
                         + R_MNabc_ther + R_TMDK2_ther - R_HHYR2_ther
                         + R_TRPt6_ther + R_AGPR_ther - R_TDMAT6_ther
                         + R_DNMPPA_ther - R_MCMAT3_ther + R_GAPD_ther
                         - R_KARA2_ther - 2 R_NOX2_ther + R_PHETRS_ther
                         - R_G1PTMT_ther + R_ALATRS_ther + R_THRt6_ther
                         + R_NADK_ther + R_CYTDK1_ther + R_GLYBabc_ther
                         + R_PPA_ther - R_GALU_ther + R_METACH_ther
                         + R_METS_ther + R_CYTDK2_ther + R_FRUK_ther
                         - R_G5SD_ther + R_CITt6_ther - R_CYTD_ther
                         - R_DASYN_LLA_ther + 2 R_DAGGT_LLA_ther + R_HISt6_ther
                         + R_DAGK_LLA_ther + 25 R_GRTAL_ther + R_DHFS_ther
                         + R_UAMAGS_ther - R_2HXICDH_ther - R_OCMAT3_ther
                         + 2 R_SSALx_ther + R_THFGLUS_ther + 2 R_GMPS2_ther
                         + R_ALAD_Lr_ther + R_FERRICHabc_ther + R_DUTPDP_ther
                         + R_METSabc_ther - R_TRDR_ther + R_GLYTRS_ther
                         - R_ALCD19_ther + R_ALA_Lt6_ther + R_PRPPS_ther
                         + R_DURIK1_ther + 18.15 R_biomass_STR_ther
                         + R_ACPS2_ther - R_HEMAT2_ther - R_ASPTA4_ther
                         - R_DHFR_ther + R_ALCD2x_ther + R_G5SADs_ther
                         - R_DHDPRy_ther + R_SUCCt6_ther - R_OCDMAT8_ther
                         + R_GLUt6_ther + 0.4 R_RNAS_LLA_ther - R_HTDR6_ther
                         - R_PPNDH_ther + R_FOLt_ther + R_SHKK_ther
                         - R_PPCDC_ther - R_HDMAT7_ther + R_GTHRD_ther
                         + R_UAAGDS_ther + R_URIK2_ther - R_MCMAT4_ther
                         + R_RBFK_ther + 2 R_NPPDO_ther + R_GLUabc_ther
                         + 2 R_UAAGLS1_ther + R_LYSt6_ther - R_NAPRT_ther
                         + R_KARA1_ther + 9.8 R_GALTAL_ther + R_ALAALA_ther
                         + R_GARFT_ther + R_PHEt6_ther - R_ECOAH2c_ther
                         + R_IPMD_ther + R_PFK_ther + R_HEX1_ther + R_ACt6_ther
                         - R_PYK_ther + R_MALTTRabc_ther + R_CYSS_ther
                         - R_FLVR_ther + 2 R_ADSS_ther + R_PGSA_LLA_ther
                         + R_HCYSMT_ther + R_TMDK1_ther + R_ARGTRS_ther
                         + R_AHSERL4_ther + R_URAt2_ther - R_UAGDP_ther
                         + R_PPC_ther + R_UDCPK_ther + R_ATPM_ther
                         + R_FORt2_ther + R_ASNTRS_ther + R_CS_ther
                         + R_ILETRS_ther + R_L_LACt2_ther + R_GTPCI_ther
                         + R_PROLBabc_ther + R_LEUabc_ther + R_SERTRS_ther
                         + R_ARGabc_ther + R_AIRC_ther + 2 R_HMGCOAR_ther
                         + R_TYRt6_ther + R_TYRTRS_ther + R_HPPK_ther
                         - R_BTMAT1_ther + R_PNTOt2_ther - R_ACHBS_ther
                         + R_LEUt6_ther - R_DCYTD_ther + R_PYDAMt_ther
                         + 0.111 R_CASabc_ther - 2 R_UREA2_ther - R_HBUR1_ther
                         + R_LEUTRS_ther - R_PUNP8_ther - R_MCMAT5_ther
                         + R_IMPD_ther + R_DPCOAK_ther - 2 R_P5CRr_ther
                         + R_PRASCS_ther + 2 R_CBPS_ther + R_LYSTRS_ther
                         + R_METTRS_ther - R_ADPDS_ther - R_NOX1_ther
                         + R_GALK2_ther + R_ASAD_ther + R_UGT1_LLA_ther
                         + R_PROt6_ther + R_METDabc_2_ther + R_UAMAS_ther
                         + R_Kt2_ther + R_IPPS_ther + 25 R_GLTAL_ther
                         + R_HCO3E_ther - R_NNAT_ther - R_ADA_ther
                         + R_GLYt6_ther + R_TRPTRS_ther - R_kaasIII_ther
                         + R_SHSL5_ther + R_MALTabc_ther + R_RIBFLVt2_ther
                         + R_ASPTRS_ther + R_PROTRS_ther + 2 R_SSALy_ther
                         + R_ARGt_ther + 5 R_CPSS_LLA_ther + R_GLNTAL_ther
                         - R_DADA_ther - R_MCMAT7_ther + R_UGT2_LLA_ther
                         - R_PYROX_ther + R_HMPK1_ther - R_HODR8_ther
                         + 2 R_CTPS1_1_ther - 2 R_CPPPGO_ther + R_PGCD_ther
                         + R_ILEabc_ther + R_METt6_ther  = 0
 M_h_e:                  - R_GUAt2_bulg - R_ASPt6_bulg - R_MELIBt2_bulg
                         - R_L_LACt2_bulg - R_URAt2_bulg - R_TYRt6_bulg
                         - R_VALt6_ther - R_LYSt6_bulg - R_LCTSt6_bulg
                         - R_SERt6_ther - R_LCTSt6_ther - R_ILEt6_ther
                         - 3 R_ATPS3r_ther - R_BTNt2i_bulg - R_GLUt6_bulg
                         - R_MALt6_bulg - R_ARGt_bulg - R_HXANt2_bulg
                         - R_MALt6_ther - R_HISt6_bulg - R_ASNt6_bulg
                         - R_FOLt_bulg - R_ACt6_bulg - R_GALt6_bulg
                         - 3 R_ATPS3r_bulg - R_CYSt6_ther - R_CYSt6_bulg
                         - R_TRPt6_ther - R_PROt6_bulg - R_PHEt6_bulg
                         + R_NH4DISex_bulg - R_THRt6_ther - R_GLYt6_bulg
                         - R_CITt6_ther - R_HISt6_ther - R_SERt6_bulg
                         - R_VALt6_bulg - R_XANt2_bulg - R_PYDAMt_bulg
                         - R_ALA_Lt6_ther - R_SUCCt6_bulg - R_SUCCt6_ther
                         - R_GLUt6_ther - R_FOLt_ther - R_RIBFLVt2_bulg
                         - R_LYSt6_ther - R_PHEt6_ther - R_THRt6_bulg
                         - R_ACt6_ther - R_ALA_Lt6_bulg - R_EX_h_e_
                         - R_URAt2_ther - R_FORt2_ther - R_L_LACt2_ther
                         - R_TRPt6_bulg - R_TYRt6_ther - R_PNTOt2_ther
                         - R_LEUt6_ther - R_PYDAMt_ther - R_LEUt6_bulg
                         - R_PROt6_ther - R_Kt2_ther - R_METt6_bulg
                         - R_GLYt6_ther - R_RIBFLVt2_ther - R_PNTOt2_bulg
                         - R_ARGt_ther - R_gln_export_bulg - R_ILEt6_bulg
                         - R_METt6_ther  = 0
 M_hco3_c_bulg:          - R_ACCOAC_bulg + R_HCO3E_bulg - R_CBPS_bulg  = 0
 M_hco3_c_ther:          - R_ACCOAC_ther - R_CBPS_ther + R_HCO3E_ther  = 0
 M_hcys_L_c_bulg:        R_RHC_bulg + R_CYSTL_bulg - R_METS_bulg
                         - R_MTHPTGHM_bulg - R_HCYSMT_bulg  = 0
 M_hcys_L_c_ther:        R_CYSTL_ther + R_METACH_ther - R_METS_ther
                         - R_MTHPTGHM_ther - R_HCYSMT_ther + R_RHC_ther  = 0
 M_hdca_c_bulg:          R_KAS8_bulg  = 0
 M_hdcea_c_bulg:         R_KAS7_bulg  = 0
 M_hdeacp_c_bulg:        - 0.26 R_AGAT_LPL_bulg + R_HDMAT7_bulg - R_MCMAT8_bulg
                         - 0.26 R_GAT1_LPL_bulg  = 0
 M_hdeacp_c_ther:        - R_MCMAT8_ther - 0.295 R_AGAT_LLA_ther
                         + R_HDMAT7_ther - 0.295 R_GAT1_LLA_ther  = 0
 M_hexacp_c_bulg:        R_HEMAT2_bulg - R_MCMAT3_bulg  = 0
 M_hexacp_c_ther:        - R_MCMAT3_ther + R_HEMAT2_ther  = 0
 M_hibcoa_c_ther:        R_ECOAH2c_ther  = 0
 M_his_L_c_bulg:         R_HISt6_bulg + M_his_L_c_bulg_src - R_HISTRS_bulg
                         + 0.0646784210833032 R_CASL_bulg + 0.024 R_PEPTL_bulg
                          = 0
 M_his_L_c_ther:         - R_HISTRS_ther + R_HISt6_ther
                         + 0.0646784210833032 R_CASL_ther + 0.024 R_PEPTL_ther
                         - R_HISD1_ther  = 0
 M_his_L_e:              - R_HISt6_bulg - R_HISt6_ther - R_EX_his_L_e_  = 0
 M_histrna_c_bulg:       R_HISTRS_bulg - 0.017 R_PROTS_LPL_v60_bulg  = 0
 M_histrna_c_ther:       R_HISTRS_ther - 0.017 R_PROTS_STR_ther  = 0
 M_hmgcoa_c_bulg:        R_HMGCOAR_bulg - R_HMGCOAS_bulg  = 0
 M_hmgcoa_c_ther:        - R_HMGCOAS_ther + R_HMGCOAR_ther  = 0
 M_hom_L_c_bulg:         - R_HSK_bulg - R_HSAT_bulg - R_HSDy_bulg  = 0
 M_hom_L_c_ther:         - R_HSDy_ther - R_HSK_ther - R_HSAT_ther  = 0
 M_hxan_c_bulg:          R_PUNP5_bulg + R_HXANt2_bulg - R_HXPRT_bulg
                         + R_PUNP6_bulg + R_PNS3_bulg  = 0
 M_hxan_c_ther:          - R_HXPRT_ther + R_PUNP5_ther + R_PUNP6_ther  = 0
 M_hxan_e:               - R_EX_hxan_e_ - R_HXANt2_bulg  = 0
 M_ibcoa_c_bulg:         - R_KAS5_bulg - R_KAS1_bulg - R_KAS6_bulg  = 0
 M_icit_c_ther:          - R_ICDHx_ther + R_ACONT_ther - R_ICDHy_ther  = 0
 M_idp_c_bulg:           R_CYTDK3_bulg - R_NDPK9_bulg + R_URIK3_bulg  = 0
 M_idp_c_ther:           R_CYTDK3_ther + R_URIK3_ther - R_NDPK9_ther  = 0
 M_ile_L_c_bulg:         M_ile_L_c_bulg_src + 0.0206098522514103 R_CASL_bulg
                         + 0.041 R_PEPTL_bulg - R_ILETA_bulg - R_ILETRS_bulg
                         + R_ILEt6_bulg  = 0
 M_ile_L_c_ther:         R_ILEt6_ther - R_ILETA_ther
                         + 0.0206098522514103 R_CASL_ther - R_ILETRS_ther
                         + 0.041 R_PEPTL_ther + R_ILEabc_ther  = 0
 M_ile_L_e:              - R_EX_ile_L_e_ - R_ILEt6_ther - R_ILEt6_bulg
                         - R_ILEabc_ther  = 0
 M_iletrna_c_bulg:       - 0.043 R_PROTS_LPL_v60_bulg + R_ILETRS_bulg  = 0
 M_iletrna_c_ther:       - 0.036 R_PROTS_STR_ther + R_ILETRS_ther  = 0
 M_imp_c_bulg:           R_HXPRT_bulg - R_ADSS_bulg - R_IMPD_bulg + R_GMPR_bulg
                         - R_IMPC_bulg  = 0
 M_imp_c_ther:           - R_NTD11_ther + R_HXPRT_ther + R_GMPR_ther
                         - R_ADSS_ther - R_IMPD_ther - R_IMPC_ther  = 0
 M_indpyr_c_bulg:        R_AKGTA3_bulg  = 0
 M_indpyr_c_ther:        R_AKGTA3_ther  = 0
 M_inost_c_bulg:         R_GPDDA5_bulg  = 0
 M_inost_c_ther:         R_MI1PP_ther  = 0
 M_ins_c_bulg:           - R_PUNP5_bulg + R_ADA_bulg - R_PNS3_bulg  = 0
 M_ins_c_ther:           R_NTD11_ther - R_PUNP5_ther + R_ADA_ther  = 0
 M_ipdp_c_bulg:          - 8 R_UDCPDPS_bulg - R_IPDDI_bulg - R_DMATT_bulg
                         + R_DPMVD_bulg - R_GRTT_bulg  = 0
 M_ipdp_c_ther:          R_DPMVD_ther - 8 R_UDCPDPS_ther - R_IPDDI_ther
                         - R_GRTT_ther - R_DMATT_ther  = 0
 M_itp_c_bulg:           - R_CYTDK3_bulg + R_NDPK9_bulg - R_URIK3_bulg  = 0
 M_itp_c_ther:           - R_CYTDK3_ther - R_URIK3_ther + R_NDPK9_ther  = 0
 M_ivcoa_c_bulg:         - R_KAS3_bulg - R_KAS9_bulg - R_KAS11_bulg  = 0
 M_k_c_ther:             R_Kt2_ther  = 0
 M_k_e:                  - R_Kt2_ther  = 0
 M_lac_D_c_bulg:         - R_ALALAC_bulg - R_LDH_D_bulg - R_LOXD_bulg  = 0
 M_lac_L_c_bulg:         R_L_LACt2_bulg - R_LDH_L_bulg  = 0
 M_lac_L_c_ther:         - R_LDH_L_ther + R_L_LACt2_ther - R_MALLACap_ther  = 0
 M_lac_L_e:              - R_L_LACt2_bulg - R_L_LACt2_ther - R_EX_lac_L_e_
                         + R_MALLACap_ther  = 0
 M_lcts_c_bulg:          R_LCTSt6_bulg - R_LACZ_bulg  = 0
 M_lcts_c_ther:          R_LCTSt6_ther + R_LCTSGALex_ther - R_LACZ_ther  = 0
 M_lcts_e:               - R_LCTSt6_bulg - R_LCTSt6_ther - R_LCTSGALex_ther
                         - R_EX_lcts_e_  = 0
 M_leu_L_c_bulg:         - R_LEUTRS_bulg - R_LEUTA_bulg + M_leu_L_c_bulg_src
                         + 0.0751011189004616 R_CASL_bulg + 0.09 R_PEPTL_bulg
                         + R_LEUt6_bulg  = 0
 M_leu_L_c_ther:         0.0751011189004616 R_CASL_ther + R_LEUabc_ther
                         + 0.09 R_PEPTL_ther + R_LEUt6_ther - R_LEUTRS_ther
                         - R_LEUTA_ther  = 0
 M_leu_L_e:              - R_EX_leu_L_e_ - R_LEUabc_ther - R_LEUt6_ther
                         - R_LEUt6_bulg  = 0
 M_leutrna_c_bulg:       R_LEUTRS_bulg - 0.078 R_PROTS_LPL_v60_bulg  = 0
 M_leutrna_c_ther:       - 0.072 R_PROTS_STR_ther + R_LEUTRS_ther  = 0
 M_lgt_S_c_ther:         R_LGTHL_ther  = 0
 M_lipidX_c_ther:        R_USHD_ther  = 0
 M_lys_L_c_bulg:         - R_UAAGLS1_bulg + R_DAPDC_bulg + R_LYSt6_bulg
                         + M_lys_L_c_bulg_src - R_LYSDC_bulg
                         + 0.0273245227741142 R_CASL_bulg + 0.076 R_PEPTL_bulg
                         - R_LYSTRS_bulg  = 0
 M_lys_L_c_ther:         R_DAPDC_ther + 0.0273245227741142 R_CASL_ther
                         - R_UAAGLS1_ther + R_LYSt6_ther + 0.076 R_PEPTL_ther
                         - R_LYSTRS_ther  = 0
 M_lys_L_e:              - R_EX_lys_L_e_ - R_LYSt6_bulg - R_LYSt6_ther  = 0
 M_lyspg_LLA_c_ther:     - 0.000024 R_biomass_STR_ther + 0.01 R_LPGS_LLA_ther
                          = 0
 M_lyspg_LPL_c_bulg:     - 0.000164 R_biomass_LBUL10_bulg
                         + 0.01 R_LPGS_LPL_bulg  = 0
 M_lystrna_c_bulg:       R_LYSTRS_bulg - 0.066 R_PROTS_LPL_v60_bulg
                         - R_LPGS_LPL_bulg  = 0
 M_lystrna_c_ther:       - 0.051 R_PROTS_STR_ther - R_LPGS_LLA_ther
                         + R_LYSTRS_ther  = 0
 M_m12dg_LLA_c_ther:     - 0.000023 R_biomass_STR_ther + 0.01 R_UGT1_LLA_ther
                         - 0.01 R_UGT2_LLA_ther  = 0
 M_mal_L_c_bulg:         R_MALt6_bulg + R_FUM_bulg  = 0
 M_mal_L_c_ther:         R_MALt6_ther + R_MALLACap_ther  = 0
 M_mal_L_e:              - R_MALt6_bulg - R_MALt6_ther - R_EX_mal_L_e_
                         - R_MALLACap_ther  = 0
 M_malacp_c_bulg:        - R_MCMAT2_bulg - R_MCMAT5_bulg - R_MCMAT8_bulg
                         - R_MCMAT7_bulg - R_kaasIII_bulg + R_MACPMT_bulg
                         - R_MCMAT3_bulg - R_MCMAT4_bulg - R_MCMAT6_bulg  = 0
 M_malacp_c_ther:        - R_MCMAT2_ther - R_MCMAT6_ther + R_MACPMT_ther
                         - R_MCMAT8_ther - R_MCMAT3_ther - R_MCMAT4_ther
                         - R_MCMAT5_ther - R_kaasIII_ther - R_MCMAT7_ther  = 0
 M_malcoa_c_bulg:        - 6 R_KAS5_bulg + R_ACCOAC_bulg - 6 R_KAS2_bulg
                         - 5 R_KAS1_bulg - 5 R_KAS3_bulg - 6 R_KAS9_bulg
                         - 5 R_KAS4_bulg - 6 R_KAS10_bulg - 7 R_KAS7_bulg
                         - 8 R_KAS13_bulg - 7 R_KAS8_bulg - R_MACPMT_bulg
                         - 6 R_KAS12_bulg - 6 R_KAS11_bulg - 6 R_KAS6_bulg
                         - 8 R_KAS17_bulg  = 0
 M_malcoa_c_ther:        R_ACCOAC_ther - R_MACPMT_ther  = 0
 M_malt_c_bulg:          - R_MALT_bulg + R_MALTabc_bulg  = 0
 M_malt_c_ther:          - R_AMALT4_ther - R_AMALT2_ther - R_AMALT3_ther
                         - R_AMALT1_ther + R_MALTabc_ther  = 0
 M_malt_e:               - R_MALTabc_bulg - R_MALTabc_ther  = 0
 M_malthp_c_ther:        R_AMALT4_ther  = 0
 M_malthx_c_ther:        R_MALTHXabc_ther - R_AMALT4_ther + R_AMALT3_ther  = 0
 M_malthx_e:             - R_MALTHXabc_ther  = 0
 M_maltpt_c_ther:        R_MALTPTabc_ther + R_AMALT2_ther - R_AMALT3_ther  = 0
 M_maltpt_e:             - R_MALTPTabc_ther  = 0
 M_malttr_c_ther:        R_MALTTRabc_ther - R_AMALT1_ther  = 0
 M_malttr_e:             - R_MALTTRabc_ther  = 0
 M_maltttr_c_ther:       R_MALTTTRabc_ther - R_AMALT2_ther + R_AMALT1_ther  = 0
 M_maltttr_e:            - R_MALTTTRabc_ther  = 0
 M_man1p_c_bulg:         - R_PMANM_bulg  = 0
 M_man1p_c_ther:         - R_PMANM_ther  = 0
 M_man6p_c_bulg:         R_PMANM_bulg + R_MANpts_bulg - R_MAN6PI_bulg  = 0
 M_man6p_c_ther:         R_PMANM_ther + R_MANpts_ther - R_MAN6PI_ther  = 0
 M_man_e:                - R_MANpts_ther - R_MANpts_bulg  = 0
 M_melib_c_bulg:         R_MELIBt2_bulg  = 0
 M_melib_e:              - R_MELIBt2_bulg  = 0
 M_mercppyr_c_ther:      R_ASPTA4_ther  = 0
 M_met_D_2_c_ther:       R_METDabc_2_ther  = 0
 M_met_D_2_e:            - R_METDabc_2_ther  = 0
 M_met_L_c_bulg:         - R_METAT_bulg + R_FMETDF_bulg + M_met_L_c_bulg_src
                         - R_METTRS_bulg + R_METS_bulg + R_MTHPTGHM_bulg
                         + 0.0737678055718706 R_CASL_bulg + R_HCYSMT_bulg
                         + 0.023 R_PEPTL_bulg + R_METt6_bulg  = 0
 M_met_L_c_ther:         - R_METAT_ther + R_METabc_ther + R_METS_ther
                         + R_MTHPTGHM_ther + 0.0737678055718706 R_CASL_ther
                         + R_HCYSMT_ther + 0.023 R_PEPTL_ther - R_METTRS_ther
                         + R_FMETDF_ther + R_METt6_ther  = 0
 M_met_L_e:              - R_METabc_ther - R_EX_met_L_e_ - R_METt6_bulg
                         - R_METt6_ther  = 0
 M_methf_c_bulg:         R_MTHFD_bulg + R_FTHFCL_bulg - R_MTHFC_bulg  = 0
 M_methf_c_ther:         - R_MTHFC_ther + R_MTHFD_ther + R_FTHFCL_ther  = 0
 M_metox_R_c_ther:       R_METSRabc_ther  = 0
 M_metox_R_e:            - R_METSRabc_ther  = 0
 M_metox_c_ther:         R_METSabc_ther  = 0
 M_metox_e:              - R_METSabc_ther  = 0
 M_mettrna_c_bulg:       R_METTRS_bulg - 0.022 R_PROTS_LPL_v60_bulg  = 0
 M_mettrna_c_ther:       - 0.015 R_PROTS_STR_ther + R_METTRS_ther  = 0
 M_mev_R_c_bulg:         - R_HMGCOAR_bulg - R_MEVK_bulg  = 0
 M_mev_R_c_ther:         - R_MEVK_ther - R_HMGCOAR_ther  = 0
 M_mi1p_D_c_ther:        - R_MI1PP_ther  = 0
 M_mlthf_c_bulg:         - R_TMDS_bulg + R_METFR_bulg - R_MTHFD_bulg  = 0
 M_mlthf_c_ther:         R_METFR_ther - R_TMDS_ther - R_MTHFD_ther
                         + R_GHMT_ther  = 0
 M_mn2_c_bulg:           R_MNabc_bulg  = 0
 M_mn2_c_ther:           R_MNabc_ther  = 0
 M_mn2_e:                - R_MNabc_bulg - R_MNabc_ther  = 0
 M_mnl1p_c_bulg:         R_MNLpts_bulg  = 0
 M_mnl1p_c_ther:         R_MNLpts_ther  = 0
 M_mnl_e:                - R_MNLpts_ther - R_MNLpts_bulg  = 0
 M_mpept_c_bulg:         0.001 R_PEPTabc_bulg - 0.001 R_PEPTL_bulg  = 0
 M_mpept_c_ther:         0.001 R_PEPTabc_ther - 0.001 R_PEPTL_ther  = 0
 M_mpept_e:              - 0.001 R_PEPTabc_bulg - 0.001 R_PEPTabc_ther
                         + R_MPROTL_bulg - R_EX_mpept_e_  = 0
 M_mprot_e:              - R_MPROTL_bulg - R_EX_mprot_e_  = 0
 M_mql7_c_bulg:          R_SUCD3BC_bulg  = 0
 M_mqn7_c_bulg:          - R_SUCD3BC_bulg  = 0
 M_mthgxl_c_ther:        - R_LGTHL_ther  = 0
 M_myrt_c_bulg:          R_KAS2_bulg  = 0
 M_n2ppn_c_ther:         - 2 R_NPPDO_ther  = 0
 M_na1_c_bulg:           - R_NAabcO_bulg  = 0
 M_na1_c_ther:           R_BILEt4_ther  = 0
 M_na1_e:                R_NAabcO_bulg - R_BILEt4_ther  = 0
 M_nac_c_bulg:           R_NACUP_bulg + R_NNAMr_bulg - R_NAPRT_bulg  = 0
 M_nac_c_ther:           R_NACUP_ther - R_NAMNPP_ther + R_NNAMr_ther
                         - R_NAPRT_ther  = 0
 M_nac_e:                - R_NACUP_ther - R_NACUP_bulg - R_EX_nac_e_  = 0
 M_nad_c_bulg:           - R_SSALx_bulg - R_PGCD_bulg + R_HDMAT7_bulg
                         + R_FRDx_bulg - R_ACALD_bulg - R_GAPD_bulg
                         - R_NADK_bulg + R_OCDMAT8_bulg + R_2HXICDH_bulg
                         + R_HEMAT2_bulg + R_DEMAT4_bulg
                         - 0.002 R_biomass_LBUL10_bulg + R_BTMAT1_bulg
                         + R_TDMAT6_bulg - R_G3PD1_bulg - R_DHORD6_bulg
                         + R_NADS1_bulg - R_LDH_D_bulg + R_DDMAT5_bulg
                         - R_IMPD_bulg + R_OCMAT3_bulg - R_LDH_L_bulg
                         - R_ALCD2x_bulg  = 0
 M_nad_c_ther:           - R_ICDHx_ther + R_NADS1_ther + R_FRDx_ther
                         - R_DHORD6_ther + R_DDMAT5_ther - R_LDH_L_ther
                         + R_DEMAT4_ther - R_G3PD1_ther + R_TDMAT6_ther
                         - R_GAPD_ther + 2 R_NOX2_ther - R_NADK_ther
                         + R_2HXICDH_ther + R_OCMAT3_ther - R_SSALx_ther
                         - R_ALAD_Lr_ther + R_ALCD19_ther
                         - 0.0018 R_biomass_STR_ther + R_HEMAT2_ther
                         - R_ALCD2x_ther + R_OCDMAT8_ther + R_HDMAT7_ther
                         - R_IPMD_ther - R_PDH_ther + R_BTMAT1_ther
                         - R_IMPD_ther + R_NOX1_ther - R_PPND_ther
                         - R_PGCD_ther  = 0
 M_nadh_c_bulg:          R_SSALx_bulg + R_PGCD_bulg - R_HDMAT7_bulg
                         - R_FRDx_bulg + R_ACALD_bulg + R_GAPD_bulg
                         - R_OCDMAT8_bulg - R_2HXICDH_bulg - R_HEMAT2_bulg
                         - R_DEMAT4_bulg - R_BTMAT1_bulg - R_TDMAT6_bulg
                         + R_G3PD1_bulg + R_DHORD6_bulg + R_LDH_D_bulg
                         - R_DDMAT5_bulg + R_IMPD_bulg - R_OCMAT3_bulg
                         + R_LDH_L_bulg + R_ALCD2x_bulg  = 0
 M_nadh_c_ther:          R_ICDHx_ther - R_FRDx_ther + R_DHORD6_ther
                         - R_DDMAT5_ther + R_LDH_L_ther - R_DEMAT4_ther
                         + R_G3PD1_ther - R_TDMAT6_ther + R_GAPD_ther
                         - 2 R_NOX2_ther - R_2HXICDH_ther - R_OCMAT3_ther
                         + R_SSALx_ther + R_ALAD_Lr_ther - R_ALCD19_ther
                         - R_HEMAT2_ther + R_ALCD2x_ther - R_OCDMAT8_ther
                         - R_HDMAT7_ther + R_IPMD_ther + R_PDH_ther
                         - R_BTMAT1_ther + R_IMPD_ther - R_NOX1_ther
                         + R_PPND_ther + R_PGCD_ther  = 0
 M_nadp_c_bulg:          11 R_KAS5_bulg + R_TRDR_bulg - R_METFR_bulg
                         + R_HOCR3_bulg - 2 R_HMGCOAR_bulg - R_MTHFD_bulg
                         + R_ALCD19_D_bulg + R_NADK_bulg + R_UAPGR_bulg
                         + R_HHYR2_bulg + R_ALCD19_L_bulg - R_THZPSN2_bulg
                         + 12 R_KAS2_bulg - R_HSDy_bulg + 10 R_KAS1_bulg
                         - R_PGDH_bulg - R_ASAD_bulg - R_TDPDRR_bulg
                         + 10 R_KAS3_bulg + R_HBUR1_bulg + 11 R_KAS9_bulg
                         + R_DHFR_bulg + 10 R_KAS4_bulg + 11 R_KAS10_bulg
                         + 13 R_KAS7_bulg + 16 R_KAS13_bulg - R_G6PDHy_bulg
                         + R_HTDR6_bulg + 14 R_KAS8_bulg + R_HODR8_bulg
                         + R_HHDR7_bulg + R_GMPR_bulg + R_HDER4_bulg
                         + 12 R_KAS12_bulg + 12 R_KAS11_bulg - R_DHFOR2_bulg
                         + R_HDDR5_bulg - R_SSALy_bulg + 12 R_KAS6_bulg
                         + 15 R_KAS17_bulg  = 0
 M_nadp_c_ther:          R_HHDR7_ther + R_HDER4_ther - R_TDPDRR_ther
                         + R_SHK3D_ther + R_UAPGR_ther + R_GMPR_ther
                         - R_GLUDy_ther - R_HSDy_ther + R_HOCR3_ther
                         - R_METFR_ther + R_HDDR5_ther + R_DPR_ther
                         - R_DHFOR2_ther + R_HHYR2_ther - R_AGPR_ther
                         + R_KARA2_ther + R_NADK_ther + R_G5SD_ther
                         + R_TRDR_ther + R_DHFR_ther + R_DHDPRy_ther
                         + R_HTDR6_ther - R_GTHRD_ther - R_KARA1_ther
                         + R_FLVR_ther - R_ICDHy_ther - 2 R_HMGCOAR_ther
                         - R_MTHFD_ther + R_HBUR1_ther + R_P5CRr_ther
                         - R_ASAD_ther - R_SSALy_ther + R_HODR8_ther  = 0
 M_nadph_c_bulg:         - 11 R_KAS5_bulg - R_TRDR_bulg + R_METFR_bulg
                         - R_HOCR3_bulg + 2 R_HMGCOAR_bulg + R_MTHFD_bulg
                         - R_ALCD19_D_bulg - R_UAPGR_bulg - R_HHYR2_bulg
                         - R_ALCD19_L_bulg + R_THZPSN2_bulg - 12 R_KAS2_bulg
                         + R_HSDy_bulg - 10 R_KAS1_bulg + R_PGDH_bulg
                         + R_ASAD_bulg + R_TDPDRR_bulg - 10 R_KAS3_bulg
                         - R_HBUR1_bulg - 11 R_KAS9_bulg - R_DHFR_bulg
                         - 10 R_KAS4_bulg - 11 R_KAS10_bulg - 13 R_KAS7_bulg
                         - 16 R_KAS13_bulg + R_G6PDHy_bulg - R_HTDR6_bulg
                         - 14 R_KAS8_bulg - R_HODR8_bulg - R_HHDR7_bulg
                         - R_GMPR_bulg - R_HDER4_bulg - 12 R_KAS12_bulg
                         - 12 R_KAS11_bulg + R_DHFOR2_bulg - R_HDDR5_bulg
                         + R_SSALy_bulg - 12 R_KAS6_bulg - 15 R_KAS17_bulg  = 0
 M_nadph_c_ther:         - R_HHDR7_ther - R_HDER4_ther + R_TDPDRR_ther
                         - R_SHK3D_ther - R_UAPGR_ther - R_GMPR_ther
                         + R_GLUDy_ther + R_HSDy_ther - R_HOCR3_ther
                         + R_METFR_ther - R_HDDR5_ther - R_DPR_ther
                         + R_DHFOR2_ther - R_HHYR2_ther + R_AGPR_ther
                         - R_KARA2_ther - R_G5SD_ther - R_TRDR_ther
                         - R_DHFR_ther - R_DHDPRy_ther - R_HTDR6_ther
                         + R_GTHRD_ther + R_KARA1_ther - R_FLVR_ther
                         + R_ICDHy_ther + 2 R_HMGCOAR_ther + R_MTHFD_ther
                         - R_HBUR1_ther - R_P5CRr_ther + R_ASAD_ther
                         + R_SSALy_ther - R_HODR8_ther  = 0
 M_ncam_c_bulg:          - R_NNAMr_bulg  = 0
 M_ncam_c_ther:          - R_NMNS_ther - R_NNAMr_ther  = 0
 M_nebari_c_ther:        - R_PUNP8_ther  = 0
 M_nh3_c_bulg:           R_NH4DIS_bulg + R_NH3t_bulg  = 0
 M_nh3_e:                R_NH4DISex_bulg - R_NH3t_bulg  = 0
 M_nh4_c_bulg:           R_G6PDA_bulg - R_GLNS_bulg - R_CTPS1_1_bulg
                         + R_CYSTGL_bulg + R_CYSTL_bulg + R_ASPT_bulg
                         + R_CSND_bulg - R_NH4DIS_bulg + R_CYTD_bulg
                         + R_ADA_bulg - R_ASNS2_bulg + R_DADA_bulg
                         + R_NNAMr_bulg - R_NADS1_bulg + R_GMPR_bulg
                         + R_ASNN_bulg + R_DCYTD_bulg  = 0
 M_nh4_c_ther:           R_CYSTL_ther + R_G6PDA_ther - R_NADS1_ther
                         + R_NH4t_ther - R_ASNS2_ther + R_GMPR_ther
                         + R_GLUDy_ther + R_NNAMr_ther + R_THRD_L_ther
                         - R_GLNS_ther + R_SERD_L_ther + R_CSND_ther
                         + R_CYTD_ther + R_ALAD_Lr_ther + R_ASNN_ther
                         + R_SERD_D_ther + R_CYSTGL_ther + R_DCYTD_ther
                         + 2 R_UREA2_ther + R_ADA_ther + R_HISD1_ther
                         + R_DADA_ther - R_CTPS1_1_ther  = 0
 M_nh4_e:                - R_NH4t_ther - R_NH4DISex_bulg - R_EX_nh4_e_  = 0
 M_nicrnt_c_bulg:        - R_NNAT_bulg + R_NAPRT_bulg  = 0
 M_nicrnt_c_ther:        R_NAMNPP_ther + R_NAPRT_ther - R_NNAT_ther  = 0
 M_nmn_c_ther:           R_NMNS_ther  = 0
 M_no2_c_ther:           2 R_NPPDO_ther  = 0
 M_o2__c_ther:           - 2 R_SOD_ther  = 0
 M_o2_c_bulg:            - R_PYROX_bulg - R_LOXD_bulg  = 0
 M_o2_c_ther:            R_SOD_ther - R_NOX2_ther - R_DHORD1_ther + R_O2t_ther
                         - R_NPPDO_ther + 0.5 R_AHMMPS_ther - R_NOX1_ther
                         - R_PYROX_ther - R_CPPPGO_ther  = 0
 M_o2_e:                 - R_O2t_ther - R_EX_o2_e_  = 0
 M_oaa_c_bulg:           R_PPC_bulg  = 0
 M_oaa_c_ther:           R_ASPTA1_ther + R_PPC_ther - R_CS_ther  = 0
 M_ocdacp_c_bulg:        - 0.02 R_AGAT_LPL_bulg + R_OCDMAT8_bulg
                         - 0.02 R_GAT1_LPL_bulg  = 0
 M_ocdacp_c_ther:        - 0.01 R_AGAT_LLA_ther + R_OCDMAT8_ther
                         - 0.01 R_GAT1_LLA_ther  = 0
 M_ocdca_c_bulg:         R_KAS13_bulg  = 0
 M_ocdcea_c_bulg:        R_KAS17_bulg  = 0
 M_octacp_c_bulg:        - R_MCMAT4_bulg + R_OCMAT3_bulg  = 0
 M_octacp_c_ther:        R_OCMAT3_ther - R_MCMAT4_ther  = 0
 M_orn_L_c_ther:         - R_OCBT_ther + R_ORNTAC_ther  = 0
 M_orot5p_c_bulg:        - R_OMPDC_bulg - R_ORPT_bulg  = 0
 M_orot5p_c_ther:        - R_ORPT_ther - R_OMPDC_ther  = 0
 M_orot_c_bulg:          R_ORPT_bulg + R_DHORD6_bulg  = 0
 M_orot_c_ther:          R_ORPT_ther + R_DHORD6_ther + R_DHORD1_ther  = 0
 M_oxa_c_ther:           R_OXFOt4_ther  = 0
 M_oxa_e:                - R_OXFOt4_ther - R_EX_oxa_e_  = 0
 M_pa_LLA_c_ther:        - 0.01 R_PAP_LLA_ther - 0.01 R_DASYN_LLA_ther
                         + 0.01 R_AGAT_LLA_ther + 0.01 R_DAGK_LLA_ther  = 0
 M_pa_LPL_c_bulg:        0.01 R_AGAT_LPL_bulg - 0.01 R_DASYN_LPL_bulg  = 0
 M_pan4p_c_bulg:         R_PPCDC_bulg - R_PTPAT_bulg  = 0
 M_pan4p_c_ther:         - R_PTPAT_ther + R_PPCDC_ther  = 0
 M_pant_R_c_ther:        R_DPR_ther  = 0
 M_pap_c_bulg:           R_ACPS2_bulg  = 0
 M_pap_c_ther:           R_ACPS2_ther  = 0
 M_pep_c_bulg:           - R_ACGApts_bulg + R_ENO_bulg - R_GALTpts_bulg
                         - R_DHAPT_bulg - R_SBTpts_bulg - R_PYK_bulg
                         - R_SUCpts_bulg - R_MANpts_bulg - R_GLCpts_bulg
                         - R_TREpts_bulg - R_CELBpts_bulg - R_MNLpts_bulg
                         - R_UAGCVT_bulg - R_ACGALpts_bulg - R_PPC_bulg
                         - R_FRUpts_bulg  = 0
 M_pep_c_ther:           - R_ACGApts_ther - R_MNLpts_ther - R_MANpts_ther
                         - R_GALTpts_ther + R_ENO_ther - R_UAGCVT_ther
                         - R_SUCpts_ther - R_SBTpts_ther - R_PSCVT_ther
                         - R_TREpts_ther - R_DHAPT_ther - R_ACGALpts_ther
                         - R_PYK_ther - R_PPC_ther - R_DAHPS_ther
                         - R_CELBpts_ther - R_GLCpts_ther  = 0
 M_pg_LLA_c_ther:        0.01 R_PGPP_LLA_ther - 0.000107 R_biomass_STR_ther
                         - 0.16 R_LTAS_LLA_ther - 0.02 R_CLPNS_LLA_ther
                         - 0.01 R_LPGS_LLA_ther  = 0
 M_pg_LPL_c_bulg:        - 0.02 R_CLPNS_LPL_bulg
                         - 0.000625 R_biomass_LBUL10_bulg
                         + 0.01 R_PGPP_LPL_bulg - 0.01 R_LPGS_LPL_bulg  = 0
 M_pgp_LLA_c_ther:       - 0.01 R_PGPP_LLA_ther + 0.01 R_PGSA_LLA_ther  = 0
 M_pgp_LPL_c_bulg:       - 0.01 R_PGPP_LPL_bulg + 0.01 R_PGSA_LPL_bulg  = 0
 M_phe_L_c_bulg:         R_PHEt6_bulg - R_PHETRS_bulg - R_PHETA1_bulg
                         + 0.0309081111820145 R_CASL_bulg + 0.054 R_PEPTL_bulg
                         + M_phe_L_c_bulg_src  = 0
 M_phe_L_c_ther:         - R_ASPTA6_ther - R_PHETRS_ther
                         + 0.0309081111820145 R_CASL_ther + R_PHEt6_ther
                         + 0.054 R_PEPTL_ther  = 0
 M_phe_L_e:              - R_EX_phe_L_e_ - R_PHEt6_bulg - R_PHEt6_ther  = 0
 M_phetrna_c_bulg:       R_PHETRS_bulg - 0.034 R_PROTS_LPL_v60_bulg  = 0
 M_phetrna_c_ther:       - 0.031 R_PROTS_STR_ther + R_PHETRS_ther  = 0
 M_phom_c_bulg:          R_HSK_bulg - R_THRS_bulg  = 0
 M_phom_c_ther:          R_HSK_ther - R_THRS_ther  = 0
 M_phpyr_c_bulg:         R_PHETA1_bulg  = 0
 M_phpyr_c_ther:         R_ASPTA6_ther + R_PPNDH_ther  = 0
 M_pi_c_bulg:            0.4 R_RNAS_LPL_bulg + 25 R_DALTAL_bulg + R_METAT_bulg
                         + 1.37 R_DNAS_LPL_bulg + R_MNabc_bulg + R_PRFGS_bulg
                         + R_ACCOAC_bulg + R_UAAGLS1_bulg + 2 R_PEPTabc_bulg
                         + R_UDCPDP_bulg + R_ALALAC_bulg + R_GLYBabc_bulg
                         + R_THFGLUS_bulg + R_GLNS_bulg + R_AEPabc_bulg
                         + R_CTPS1_1_bulg + 2 R_PPA_bulg + R_UGMDDS2_bulg
                         - R_GAPD_bulg + R_CHLabc_bulg - R_PUNP3_bulg
                         + R_PRASCS_bulg + R_DHFS_bulg + R_THMabc_bulg
                         + R_THRS_bulg + R_GLNabc_bulg + R_PRAGS_bulg
                         - R_PUNP5_bulg - R_PUNP4_bulg - R_ASAD_bulg
                         + 2 R_PIabc_bulg + R_DNMPPA_bulg - R_ATPS3r_bulg
                         + 27.2 R_biomass_LBUL10_bulg + R_ZN2abc2_bulg
                         + R_UAMAS_bulg + R_UAMAGS_bulg - R_PTAr_bulg
                         + R_ADSS_bulg + R_GLNTAL_bulg - R_PUNP1_bulg
                         + R_PGLYCP_bulg + R_EPPP_bulg + R_PRAIS_bulg
                         + R_CTPS2_bulg - R_PUNP2_bulg + R_FTHFCL_bulg
                         + R_DPMVD_bulg + R_PTRCabc_bulg + R_SPMDabc_bulg
                         + R_PPNCL_bulg + R_UAGCVT_bulg + 0.111 R_CASabc_bulg
                         + R_MALTabc_bulg - R_PKL_bulg + R_ASPCT_bulg
                         + R_UAAGDS_bulg + R_ARGabc_bulg - R_PUNP6_bulg
                         + R_ATPM_bulg + R_PGPP_LPL_bulg
                         + 2.306 R_PROTS_LPL_v60_bulg + R_NAabcO_bulg
                         + R_PPC_bulg - R_PUNP7_bulg + 25 R_DARTAL_bulg
                         - R_PYROX_bulg - R_DURIPP_bulg + R_FTHFL_bulg
                         + R_CBPS_bulg  = 0
 M_pi_c_ther:            R_NTD11_ther + R_MALTTTRabc_ther + R_OCBT_ther
                         + R_PSP_L_ther - R_PUNP2_ther + R_VALabc_ther
                         + R_METAT_ther + 2 R_PIabc_ther + R_CHORS_ther
                         + R_ACCOAC_ther + R_MI1PP_ther + R_DHQS_ther
                         + R_PGPP_LLA_ther - R_GLCP_ther + R_PRAGS_ther
                         + R_NTD7_ther + R_NTD9_ther + R_PAP_LLA_ther
                         + R_NAMNPP_ther + R_DPMVD_ther + R_CTPS2_ther
                         + R_SPMDabc_ther + R_ACTPASE_ther + R_ASPCT_ther
                         + 1.37 R_DNAS_STU_ther + R_NTD10_ther + R_NTD4_ther
                         - R_PUNP7_ther + 2.306 R_PROTS_STR_ther
                         - R_ATPS3r_ther + R_ASPabc_ther + R_MALTHXabc_ther
                         + R_MALTPTabc_ther + R_PRFGS_ther + R_PPNCL_ther
                         + 2 R_PEPTabc_ther + R_ZN2abc2_ther + R_METabc_ther
                         + R_PRAIS_ther + R_GLNabc_ther + R_UDCPDP_ther
                         + R_UAGCVT_ther + R_GLNS_ther + R_PTRCabc_ther
                         + R_NTD8_ther - R_PUNP5_ther + R_THMabc_ther
                         + R_METSRabc_ther - R_PUNP3_ther + R_NTD2_ther
                         + R_SUC6PH_ther + R_MNabc_ther - R_AGPR_ther
                         + R_DNMPPA_ther + R_NTD1_ther - R_PTAr_ther
                         - R_GAPD_ther + R_PSCVT_ther + R_GLYBabc_ther
                         + 2 R_PPA_ther + R_G5SD_ther + R_DHFS_ther
                         + R_UAMAGS_ther - R_PUNP4_ther + R_THFGLUS_ther
                         + R_FERRICHabc_ther - R_DURIPP_ther + R_METSabc_ther
                         + 18.15 R_biomass_STR_ther + 0.4 R_RNAS_LLA_ther
                         + 6 R_DALTAL_LLA_ther + R_THRS_ther + R_UAAGDS_ther
                         + R_GLUabc_ther + R_UAAGLS1_ther + R_UGLDDS1_ther
                         + R_ALAALA_ther - R_PUNP6_ther + R_MALTTRabc_ther
                         + R_ADSS_ther + R_AHMMPS_ther + R_PPC_ther
                         + R_ATPM_ther + R_PROLBabc_ther + R_NTD6_ther
                         + R_LEUabc_ther + R_ARGabc_ther + 0.111 R_CASabc_ther
                         + R_DAHPS_ther - R_PUNP8_ther + R_PRASCS_ther
                         + R_CBPS_ther + R_NTD5_ther - R_ASAD_ther
                         + R_METDabc_2_ther + R_UAMAS_ther + R_FTHFCL_ther
                         + R_MALTabc_ther - R_PUNP1_ther + R_NTD3_ther
                         + R_GLNTAL_ther - R_PYROX_ther + R_CTPS1_1_ther
                         + R_FTHFL_ther + R_ILEabc_ther  = 0
 M_pi_e:                 - R_PIabc_ther - R_PIabc_bulg - R_EX_pi_e_  = 0
 M_pnto_R_c_bulg:        - R_PNTK_bulg + R_PNTOt2_bulg  = 0
 M_pnto_R_c_ther:        - R_PNTK_ther + R_PNTOt2_ther  = 0
 M_pnto_R_e:             - R_EX_pnto_R_e_ - R_PNTOt2_ther - R_PNTOt2_bulg  = 0
 M_polypi_c_bulg:        - R_EPPP_bulg  = 0
 M_pphn_c_bulg:          R_CHORM_bulg  = 0
 M_pphn_c_ther:          R_CHORM_ther - R_PPNDH_ther - R_PPND_ther  = 0
 M_ppi_c_bulg:           R_RNAS_LPL_bulg + R_METAT_bulg + R_DNAS_LPL_bulg
                         + R_ALATRS_bulg + R_CYSTRS_bulg + R_ASNTRS_bulg
                         + R_G1PTMT_bulg + R_ARGTRS_bulg + 8 R_UDCPDPS_bulg
                         + R_UPPRT_bulg + R_NNAT_bulg + R_LEUTRS_bulg
                         + R_G3PCT_bulg - R_PPA_bulg + R_SERTRS_bulg
                         + R_DHPS3_bulg + R_GUAPRT_bulg + R_FMNAT_bulg
                         + R_HXPRT_bulg + R_ASPTRS_bulg + R_DNTPPA_bulg
                         + R_VALTRS_bulg + R_METTRS_bulg + R_ASNS2_bulg
                         + R_PROTRS_bulg + R_HISTRS_bulg + R_MOADCST_bulg
                         + R_PHETRS_bulg + R_GLYTRS_bulg + R_UAGDP_bulg
                         + R_ADK2_bulg + R_DMATT_bulg - R_ORPT_bulg
                         + R_GMPS2_bulg + R_PTPAT_bulg + R_THRTRS_bulg
                         + R_DUTPDP_bulg + R_DASYN_LPL_bulg + R_NADS1_bulg
                         + R_GLUTRS_bulg + R_ADPT_bulg + R_XPPT_bulg
                         + R_LYSTRS_bulg + R_TYRTRS_bulg + R_GRTT_bulg
                         + R_GALU_bulg + R_ILETRS_bulg + R_ASNS1_bulg
                         + R_TRPTRS_bulg - R_PPIK_bulg + R_NAPRT_bulg
                         + R_GLUPRT_bulg  = 0
 M_ppi_c_ther:           R_GLUTRS_ther - R_ORPT_ther + R_HISTRS_ther
                         + R_METAT_ther + R_GALT_ther + R_NADS1_ther
                         + R_PTPAT_ther + R_ASNS2_ther + R_HXPRT_ther
                         + R_ARGSS_ther + R_NAMNPP_ther + R_ADK2_ther
                         + 8 R_UDCPDPS_ther + R_NMNS_ther + R_DNAS_STU_ther
                         + R_VALTRS_ther + R_FMNAT_ther + R_DNTPPA_ther
                         + R_CYSTRS_ther + R_GUAPRT_ther + R_THRTRS_ther
                         + R_PZS_ther + R_PHETRS_ther + R_G1PTMT_ther
                         + R_ALATRS_ther - R_PPA_ther + R_GALU_ther
                         + R_DASYN_LLA_ther + R_XPPT_ther + R_GRTT_ther
                         + R_GMPS2_ther + R_DUTPDP_ther + R_GLYTRS_ther
                         + R_RNAS_LLA_ther + R_MOADCST_ther + R_DHPS3_ther
                         + R_NAPRT_ther + R_DMATT_ther + R_ADPT_ther
                         + R_ARGTRS_ther + R_UAGDP_ther + R_ASNTRS_ther
                         + R_ILETRS_ther + R_SERTRS_ther + R_TYRTRS_ther
                         + R_LEUTRS_ther + R_GLUPRT_ther + R_LYSTRS_ther
                         + R_METTRS_ther + R_ANPRT_ther + R_NNAT_ther
                         + R_TRPTRS_ther + R_ASPTRS_ther + R_PROTRS_ther
                         + R_UPPRT_ther  = 0
 M_pppg9_c_ther:         R_CPPPGO_ther  = 0
 M_pppi_c_bulg:          - R_ADK2_bulg + R_PPIK_bulg  = 0
 M_pppi_c_ther:          - R_ADK2_ther + R_PTHPS_ther  = 0
 M_pram_c_bulg:          - R_PRAGS_bulg + R_GLUPRT_bulg  = 0
 M_pram_c_ther:          - R_PRAGS_ther + R_GLUPRT_ther  = 0
 M_pran_c_bulg:          - R_PRAI_bulg  = 0
 M_pran_c_ther:          - R_PRAI_ther + R_ANPRT_ther  = 0
 M_pro_L_c_bulg:         - R_PROTRS_bulg + R_PROt6_bulg
                         + 0.0150459290356785 R_CASL_bulg + 0.088 R_PEPTL_bulg
                         + M_pro_L_c_bulg_src  = 0
 M_pro_L_c_ther:         0.0150459290356785 R_CASL_ther + 0.088 R_PEPTL_ther
                         + R_P5CRr_ther + R_PROt6_ther - R_PROTRS_ther  = 0
 M_pro_L_e:              - R_PROt6_bulg - R_EX_pro_L_e_ - R_PROt6_ther  = 0
 M_prolb_c_ther:         R_PROLBabc_ther  = 0
 M_prolb_e:              - R_PROLBabc_ther  = 0
 M_protrna_c_bulg:       R_PROTRS_bulg - 0.04 R_PROTS_LPL_v60_bulg  = 0
 M_protrna_c_ther:       - 0.128 R_PROTS_STR_ther + R_PROTRS_ther  = 0
 M_prpp_c_bulg:          - R_UPPRT_bulg - R_GUAPRT_bulg - R_HXPRT_bulg
                         + R_PRPPS_bulg + R_ORPT_bulg - R_ADPT_bulg
                         - R_XPPT_bulg - R_NAPRT_bulg - R_GLUPRT_bulg  = 0
 M_prpp_c_ther:          R_ORPT_ther - R_HXPRT_ther - R_NAMNPP_ther
                         - R_NMNS_ther - R_GUAPRT_ther - R_XPPT_ther
                         + R_PRPPS_ther - R_NAPRT_ther - R_ADPT_ther
                         - R_GLUPRT_ther - R_ANPRT_ther - R_UPPRT_ther  = 0
 M_psd5p_c_bulg:         R_YUMPS_bulg  = 0
 M_psd5p_c_ther:         R_YUMPS_ther  = 0
 M_pser_L_c_ther:        - R_PSP_L_ther + R_PSERT_ther  = 0
 M_ptrc_c_bulg:          R_PTRCabc_bulg  = 0
 M_ptrc_c_ther:          R_AGMT_ther + R_PTRCabc_ther  = 0
 M_ptrc_e:               - R_PTRCabc_ther - R_PTRCabc_bulg  = 0
 M_purine_1_c_ther:      R_PUNP8_ther  = 0
 M_pyam5p_c_bulg:        R_PYDAMK_bulg  = 0
 M_pyam5p_c_ther:        R_PYDAMK_ther  = 0
 M_pydam_c_bulg:         - R_PYDAMK_bulg + R_PYDAMt_bulg + R_ALATA_Lr_bulg  = 0
 M_pydam_c_ther:         - R_PYDAMK_ther + R_PYDAMt_ther  = 0
 M_pydam_e:              - R_PYDAMt_bulg - R_EX_pydam_e_ - R_PYDAMt_ther  = 0
 M_pydx5p_c_bulg:        - 0.000001 R_biomass_LBUL10_bulg + R_PYDXK_bulg  = 0
 M_pydx5p_c_ther:        R_PYDXK_ther  = 0
 M_pydx_c_bulg:          - R_PYDXK_bulg - R_ALATA_Lr_bulg  = 0
 M_pydx_c_ther:          - R_PYDXK_ther  = 0
 M_pyr_c_bulg:           R_ACGApts_bulg + R_GALTpts_bulg + R_DHAPT_bulg
                         + R_CYSTL_bulg + R_SBTpts_bulg + R_PYK_bulg
                         + R_SUCpts_bulg + R_MANpts_bulg + R_GLCpts_bulg
                         + R_TREpts_bulg + R_CELBpts_bulg + R_MNLpts_bulg
                         + R_LDH_D_bulg + R_ACGALpts_bulg + R_FRUpts_bulg
                         + R_LDH_L_bulg - R_PYROX_bulg + R_LOXD_bulg
                         + R_ALATA_Lr_bulg  = 0
 M_pyr_c_ther:           R_ACGApts_ther + R_CYSTL_ther + R_ANS1_ther
                         - R_DHDPS_ther + R_MNLpts_ther + R_MANpts_ther
                         - 2 R_ACLS_ther + R_GALTpts_ther + R_LDH_L_ther
                         + R_SUCpts_ther + R_SBTpts_ther + R_SERD_L_ther
                         + R_ADCL_ther - R_PFL_ther + R_TREpts_ther
                         + R_ALAD_Lr_ther + R_SERD_D_ther + R_DHAPT_ther
                         + R_ACGALpts_ther + R_PYK_ther - R_PDH_ther
                         - R_ACHBS_ther + R_CELBpts_ther + R_GLCpts_ther
                         + R_ALATA_L_ther - R_PYROX_ther  = 0
 M_r1p_c_bulg:           R_PUNP3_bulg + R_PUNP5_bulg + R_PUNP1_bulg
                         + R_PUNP7_bulg  = 0
 M_r1p_c_ther:           - R_PPM_ther + R_PUNP7_ther + R_PUNP5_ther
                         + R_PUNP3_ther + R_PUNP8_ther + R_PUNP1_ther  = 0
 M_r5p_c_bulg:           R_RBK_bulg + R_ADPRDP_bulg - R_YUMPS_bulg
                         - R_PRPPS_bulg - R_RPI_bulg  = 0
 M_r5p_c_ther:           R_PPM_ther + R_RBK_ther - R_TKT1_ther + R_ADPRDP_ther
                         - R_PRPPS_ther - R_RPI_ther - R_YUMPS_ther  = 0
 M_rhcys_c_bulg:         - R_RHC_bulg + R_AHCYSNS_bulg  = 0
 M_rhcys_c_ther:         R_AHCYSNS_ther - R_RHC_ther  = 0
 M_rib_D_c_bulg:         - R_RBK_bulg + R_RHC_bulg + R_PNS4_bulg + R_PNS2_bulg
                         + R_PYRNS1_bulg + R_PNS3_bulg + R_PNS1_bulg  = 0
 M_rib_D_c_ther:         - R_RBK_ther + R_RHC_ther  = 0
 M_ribflvRD_c_ther:      R_FLVR_ther  = 0
 M_ribflv_c_bulg:        R_RIBFLVt2_bulg - R_RBFK_bulg  = 0
 M_ribflv_c_ther:        - R_RBFK_ther - R_FLVR_ther + R_RIBFLVt2_ther  = 0
 M_ribflv_e:             - R_RIBFLVt2_bulg - R_EX_ribflv_e_ - R_RIBFLVt2_ther
                          = 0
 M_ru5p_D_c_bulg:        R_PGDH_bulg + R_RPI_bulg - R_RPE_bulg  = 0
 M_ru5p_D_c_ther:        - R_RPE_ther + R_RPI_ther  = 0
 M_s7p_c_ther:           R_TKT1_ther  = 0
 M_sbt6p_c_bulg:         R_SBTpts_bulg  = 0
 M_sbt6p_c_ther:         R_SBTpts_ther  = 0
 M_sbt_D_e:              - R_SBTpts_bulg - R_SBTpts_ther  = 0
 M_scys_L_c_ther:        R_AHSERL3_ther  = 0
 M_seln_c_ther:          - R_AHSERL3_ther  = 0
 M_ser_D_c_ther:         - R_SERD_D_ther  = 0
 M_ser_L_c_bulg:         - R_SERTRS_bulg + R_GPDDA3_bulg + R_MOADCST_bulg
                         + R_SERt6_bulg + 0.0414598460716848 R_CASL_bulg
                         + 0.058 R_PEPTL_bulg + M_ser_L_c_bulg_src  = 0
 M_ser_L_c_ther:         R_PSP_L_ther + R_SERt6_ther - R_TRPS1_ther
                         - R_SERD_L_ther - R_SERAT_ther
                         + 0.0414598460716848 R_CASL_ther + R_MOADCST_ther
                         - R_SERTRS_ther + 0.058 R_PEPTL_ther - R_GHMT_ther
                          = 0
 M_ser_L_e:              - R_SERt6_ther - R_SERt6_bulg - R_EX_ser_L_e_  = 0
 M_sertrna_c_bulg:       R_SERTRS_bulg - 0.056 R_PROTS_LPL_v60_bulg  = 0
 M_sertrna_c_ther:       - 0.058 R_PROTS_STR_ther + R_SERTRS_ther  = 0
 M_skm5p_c_bulg:         R_SHKK_bulg  = 0
 M_skm5p_c_ther:         - R_PSCVT_ther + R_SHKK_ther  = 0
 M_skm_c_bulg:           - R_SHKK_bulg  = 0
 M_skm_c_ther:           R_SHK3D_ther - R_SHKK_ther  = 0
 M_sl26da_c_ther:        - R_SDPTA_ther  = 0
 M_sl2a6o_c_ther:        R_SDPTA_ther + R_THDPS_ther  = 0
 M_slcys_c_ther:         R_SLCYSS_ther  = 0
 M_so3_c_ther:           R_AHSERL4_ther  = 0
 M_spmd_c_bulg:          R_SPMDabc_bulg  = 0
 M_spmd_c_ther:          R_SPMDabc_ther  = 0
 M_spmd_e:               - R_SPMDabc_ther - R_SPMDabc_bulg  = 0
 M_suc6p_c_bulg:         R_SUCpts_bulg  = 0
 M_suc6p_c_ther:         R_SUCpts_ther - R_SUC6PH_ther  = 0
 M_sucbz_c_bulg:         R_SUCBZS_bulg  = 0
 M_succ_c_bulg:          R_SSALx_bulg + R_FRDx_bulg + R_SUCCt6_bulg
                         - R_SUCD3BC_bulg + R_SSALy_bulg  = 0
 M_succ_c_ther:          R_FRDx_ther + R_SSALx_ther + R_SUCCt6_ther
                         + R_SSALy_ther  = 0
 M_succ_e:               - R_EX_succ_e_ - R_SUCCt6_bulg - R_SUCCt6_ther  = 0
 M_succoa_c_ther:        - R_THDPS_ther  = 0
 M_sucr_c_ther:          R_SUC6PH_ther  = 0
 M_sucr_e:               - R_SUCpts_bulg - R_SUCpts_ther  = 0
 M_sucsal_c_bulg:        - R_SSALx_bulg - R_SSALy_bulg  = 0
 M_sucsal_c_ther:        - R_SSALx_ther - R_SSALy_ther  = 0
 M_tdeacp_c_bulg:        - 0.03 R_AGAT_LPL_bulg + R_TDMAT6_bulg - R_MCMAT7_bulg
                         - 0.03 R_GAT1_LPL_bulg  = 0
 M_tdeacp_c_ther:        R_TDMAT6_ther - 0.09 R_AGAT_LLA_ther
                         - 0.09 R_GAT1_LLA_ther - R_MCMAT7_ther  = 0
 M_thdp_c_ther:          - R_THDPS_ther + R_DHDPRy_ther - R_APAT_ther  = 0
 M_thf_c_bulg:           R_GARFT_bulg - R_THFGLUS_bulg + R_AICART_bulg
                         - 0.00001 R_biomass_LBUL10_bulg + R_DHFR_bulg
                         + R_METS_bulg - R_FTHFL_bulg  = 0
 M_thf_c_ther:           R_METS_ther - R_THFGLUS_ther
                         - 0.00001 R_biomass_STR_ther + R_DHFR_ther
                         + R_GARFT_ther + R_AICART_ther - R_GHMT_ther
                         - R_FTHFL_ther  = 0
 M_thfglu_c_bulg:        R_THFGLUS_bulg  = 0
 M_thfglu_c_ther:        R_THFGLUS_ther  = 0
 M_thglu_c_bulg:         R_MTHPTGHM_bulg  = 0
 M_thglu_c_ther:         R_MTHPTGHM_ther  = 0
 M_thm_c_bulg:           R_THMabc_bulg - R_TMDPK_bulg  = 0
 M_thm_c_ther:           - R_TMDPK_ther + R_THMabc_ther  = 0
 M_thm_e:                - R_THMabc_bulg - R_THMabc_ther - R_EX_thm_e_  = 0
 M_thmpp_c_bulg:         - 0.00001 R_biomass_LBUL10_bulg + R_TMDPK_bulg  = 0
 M_thmpp_c_ther:         R_TMDPK_ther - 0.00001 R_biomass_STR_ther  = 0
 M_thr_L_c_bulg:         R_THRS_bulg + M_thr_L_c_bulg_src - R_THRTRS_bulg
                         + R_THRt6_bulg + 0.043 R_PEPTL_bulg  = 0
 M_thr_L_c_ther:         - R_THRD_L_ther - R_THRTRS_ther + R_THRt6_ther
                         + R_THRS_ther + 0.043 R_PEPTL_ther - R_THRA_ther  = 0
 M_thr_L_e:              - R_THRt6_ther - R_EX_thr_L_e_ - R_THRt6_bulg  = 0
 M_thrtrna_c_bulg:       R_THRTRS_bulg - 0.064 R_PROTS_LPL_v60_bulg  = 0
 M_thrtrna_c_ther:       - 0.044 R_PROTS_STR_ther + R_THRTRS_ther  = 0
 M_thymd_c_bulg:         - R_TMDK1_bulg - R_TMDK2_bulg  = 0
 M_thymd_c_ther:         - R_TMDK2_ther - R_TMDK1_ther + R_NTD5_ther  = 0
 M_trdox_c_bulg:         - R_TRDR_bulg + R_RNTR2_bulg + R_RNTR3_bulg
                         + R_RNTR4_bulg + R_RNTR1_bulg  = 0
 M_trdox_c_ther:         R_RNDR3_ther + R_RNDR4_ther + R_RNTR3_ther
                         + R_RNTR1_ther + R_RNTR4_ther - R_TRDR_ther
                         + R_RNDR2_ther + R_RNDR1_ther + R_AHSERL4_ther
                         + R_RNTR2_ther  = 0
 M_trdrd_c_bulg:         R_TRDR_bulg - R_RNTR2_bulg - R_RNTR3_bulg
                         - R_RNTR4_bulg - R_RNTR1_bulg  = 0
 M_trdrd_c_ther:         - R_RNDR3_ther - R_RNDR4_ther - R_RNTR3_ther
                         - R_RNTR1_ther - R_RNTR4_ther + R_TRDR_ther
                         - R_RNDR2_ther - R_RNDR1_ther - R_AHSERL4_ther
                         - R_RNTR2_ther  = 0
 M_tre6p_c_bulg:         R_TREpts_bulg  = 0
 M_tre6p_c_ther:         - R_TRE6PH_ther + R_TREpts_ther  = 0
 M_tre_e:                - R_TREpts_bulg - R_TREpts_ther  = 0
 M_trnaala_c_bulg:       - R_ALATRS_bulg + 0.125 R_PROTS_LPL_v60_bulg  = 0
 M_trnaala_c_ther:       0.082 R_PROTS_STR_ther - R_ALATRS_ther  = 0
 M_trnaarg_c_bulg:       - R_ARGTRS_bulg + 0.04 R_PROTS_LPL_v60_bulg  = 0
 M_trnaarg_c_ther:       0.023 R_PROTS_STR_ther - R_ARGTRS_ther  = 0
 M_trnaasn_c_bulg:       - R_ASNTRS_bulg + 0.06 R_PROTS_LPL_v60_bulg  = 0
 M_trnaasn_c_ther:       0.08 R_PROTS_STR_ther - R_ASNTRS_ther  = 0
 M_trnaasp_c_bulg:       - R_ASPTRS_bulg + 0.06 R_PROTS_LPL_v60_bulg  = 0
 M_trnaasp_c_ther:       0.054 R_PROTS_STR_ther - R_ASPTRS_ther  = 0
 M_trnacys_c_bulg:       - R_CYSTRS_bulg + 0.011 R_PROTS_LPL_v60_bulg  = 0
 M_trnacys_c_ther:       0.006 R_PROTS_STR_ther - R_CYSTRS_ther  = 0
 M_trnaglu_c_bulg:       - R_GLUTRS_bulg + 0.106 R_PROTS_LPL_v60_bulg  = 0
 M_trnaglu_c_ther:       - R_GLUTRS_ther + 0.118 R_PROTS_STR_ther  = 0
 M_trnagly_c_bulg:       - R_GLYTRS_bulg + 0.084 R_PROTS_LPL_v60_bulg  = 0
 M_trnagly_c_ther:       0.054 R_PROTS_STR_ther - R_GLYTRS_ther  = 0
 M_trnahis_c_bulg:       - R_HISTRS_bulg + 0.017 R_PROTS_LPL_v60_bulg  = 0
 M_trnahis_c_ther:       - R_HISTRS_ther + 0.017 R_PROTS_STR_ther  = 0
 M_trnaile_c_bulg:       0.043 R_PROTS_LPL_v60_bulg - R_ILETRS_bulg  = 0
 M_trnaile_c_ther:       0.036 R_PROTS_STR_ther - R_ILETRS_ther  = 0
 M_trnaleu_c_bulg:       - R_LEUTRS_bulg + 0.078 R_PROTS_LPL_v60_bulg  = 0
 M_trnaleu_c_ther:       0.072 R_PROTS_STR_ther - R_LEUTRS_ther  = 0
 M_trnalys_c_bulg:       - R_LYSTRS_bulg + 0.066 R_PROTS_LPL_v60_bulg
                         + R_LPGS_LPL_bulg  = 0
 M_trnalys_c_ther:       0.051 R_PROTS_STR_ther + R_LPGS_LLA_ther
                         - R_LYSTRS_ther  = 0
 M_trnamet_c_bulg:       - R_METTRS_bulg + 0.022 R_PROTS_LPL_v60_bulg  = 0
 M_trnamet_c_ther:       0.015 R_PROTS_STR_ther - R_METTRS_ther  = 0
 M_trnaphe_c_bulg:       - R_PHETRS_bulg + 0.034 R_PROTS_LPL_v60_bulg  = 0
 M_trnaphe_c_ther:       0.031 R_PROTS_STR_ther - R_PHETRS_ther  = 0
 M_trnapro_c_bulg:       - R_PROTRS_bulg + 0.04 R_PROTS_LPL_v60_bulg  = 0
 M_trnapro_c_ther:       0.128 R_PROTS_STR_ther - R_PROTRS_ther  = 0
 M_trnaser_c_bulg:       - R_SERTRS_bulg + 0.056 R_PROTS_LPL_v60_bulg  = 0
 M_trnaser_c_ther:       0.058 R_PROTS_STR_ther - R_SERTRS_ther  = 0
 M_trnathr_c_bulg:       - R_THRTRS_bulg + 0.064 R_PROTS_LPL_v60_bulg  = 0
 M_trnathr_c_ther:       0.044 R_PROTS_STR_ther - R_THRTRS_ther  = 0
 M_trnatrp_c_bulg:       0.006 R_PROTS_LPL_v60_bulg - R_TRPTRS_bulg  = 0
 M_trnatrp_c_ther:       0.047 R_PROTS_STR_ther - R_TRPTRS_ther  = 0
 M_trnatyr_c_bulg:       - R_TYRTRS_bulg + 0.028 R_PROTS_LPL_v60_bulg  = 0
 M_trnatyr_c_ther:       0.032 R_PROTS_STR_ther - R_TYRTRS_ther  = 0
 M_trnaval_c_bulg:       - R_VALTRS_bulg + 0.06 R_PROTS_LPL_v60_bulg  = 0
 M_trnaval_c_ther:       - R_VALTRS_ther + 0.052 R_PROTS_STR_ther  = 0
 M_trp_L_c_bulg:         - R_AKGTA3_bulg + 0.0177309955194894 R_CASL_bulg
                         + M_trp_L_c_bulg_src + 0.008 R_PEPTL_bulg
                         + R_TRPt6_bulg - R_TRPTRS_bulg  = 0
 M_trp_L_c_ther:         R_TRPS1_ther + R_TRPt6_ther
                         + 0.0177309955194894 R_CASL_ther - R_AKGTA3_ther
                         + 0.008 R_PEPTL_ther - R_TRPTRS_ther  = 0
 M_trp_L_e:              - R_EX_trp_L_e_ - R_TRPt6_ther - R_TRPt6_bulg  = 0
 M_trptrna_c_bulg:       - 0.006 R_PROTS_LPL_v60_bulg + R_TRPTRS_bulg  = 0
 M_trptrna_c_ther:       - 0.047 R_PROTS_STR_ther + R_TRPTRS_ther  = 0
 M_tsul_c_ther:          - R_SLCYSS_ther - R_AHSERL4_ther  = 0
 M_tyr_L_c_bulg:         R_TYRt6_bulg + M_tyr_L_c_bulg_src - R_TYRTA_bulg
                         + 0.0563104091865341 R_CASL_bulg + 0.046 R_PEPTL_bulg
                         - R_TYRTRS_bulg  = 0
 M_tyr_L_c_ther:         - R_ASPTA5_ther + 0.0563104091865341 R_CASL_ther
                         + 0.046 R_PEPTL_ther + R_TYRt6_ther - R_TYRTRS_ther
                          = 0
 M_tyr_L_e:              - R_TYRt6_bulg - R_EX_tyr_L_e_ - R_TYRt6_ther  = 0
 M_tyrtrna_c_bulg:       R_TYRTRS_bulg - 0.028 R_PROTS_LPL_v60_bulg  = 0
 M_tyrtrna_c_ther:       - 0.032 R_PROTS_STR_ther + R_TYRTRS_ther  = 0
 M_u23ga_c_ther:         - R_USHD_ther  = 0
 M_uAgl_c_bulg:          R_UAAGLS1_bulg  = 0
 M_uAgl_c_ther:          R_UAAGLS1_ther - R_UGLDDS1_ther  = 0
 M_uAgla_c_ther:         - R_PAPPT1_ther + R_UGLDDS1_ther  = 0
 M_uaAgla_c_ther:        - R_UAGPT1_ther + R_PAPPT1_ther  = 0
 M_uaaAgla_c_ther:       R_UAGPT1_ther - R_PGGT2_ther  = 0
 M_uaagmdalac_c_bulg:    R_UAGPT5_bulg - R_PGGT4_bulg  = 0
 M_uaccg_c_bulg:         - R_UAPGR_bulg + R_UAGCVT_bulg  = 0
 M_uaccg_c_ther:         - R_UAPGR_ther + R_UAGCVT_ther  = 0
 M_uacgam_c_bulg:        - R_UAGPT5_bulg + R_UAGDP_bulg - R_ACGAMT_bulg
                         - R_UAGCVT_bulg - R_UAG2E_bulg  = 0
 M_uacgam_c_ther:        - R_UAGPT1_ther - R_UAGCVT_ther - R_ACGAMT_ther
                         + R_UAGDP_ther  = 0
 M_uacmam_c_bulg:        - R_UNAGAMAMT_bulg + R_UAG2E_bulg  = 0
 M_uagmdalac_c_bulg:     - R_UAGPT5_bulg + R_PAPPT5_bulg  = 0
 M_uama_c_bulg:          R_UAMAS_bulg - R_UAMAGS_bulg  = 0
 M_uama_c_ther:          - R_UAMAGS_ther + R_UAMAS_ther  = 0
 M_uamag_c_bulg:         - R_UAAGLS1_bulg + R_UAMAGS_bulg - R_UAAGDS_bulg  = 0
 M_uamag_c_ther:         R_UAMAGS_ther - R_UAAGDS_ther - R_UAAGLS1_ther  = 0
 M_uamr_c_bulg:          R_UAPGR_bulg - R_UAMAS_bulg  = 0
 M_uamr_c_ther:          R_UAPGR_ther - R_UAMAS_ther  = 0
 M_udcp_c_bulg:          - R_UDCPK_bulg  = 0
 M_udcp_c_ther:          - R_UDCPK_ther  = 0
 M_udcpdp_c_bulg:        - R_UDCPDP_bulg + R_UDCPDPS_bulg + R_PGGT4_bulg
                         - 0.0002 R_biomass_LBUL10_bulg  = 0
 M_udcpdp_c_ther:        R_UDCPDPS_ther - R_UDCPDP_ther
                         - 0.000018 R_biomass_STR_ther + R_PGGT2_ther  = 0
 M_udcpp_c_bulg:         R_UDCPDP_bulg - R_PAPPT5_bulg + R_UDCPK_bulg
                         + R_GTAPGL_bulg - R_ACGAMT_bulg  = 0
 M_udcpp_c_ther:         R_UDCPDP_ther - R_PAPPT1_ther - R_ACGAMT_ther
                         + R_UDCPK_ther  = 0
 M_udp_c_bulg:           R_UAGPT5_bulg + 25 R_GTAglcS_bulg + 25 R_GRTAL_bulg
                         + 2 R_DAGGT_LPL_bulg + 3 R_CPSS_LPL_bulg
                         - R_NDPK2_bulg + R_UNAGAMAMT_bulg + R_URIDK1_bulg  = 0
 M_udp_c_ther:           R_UAGPT1_ther - R_RNDR4_ther + 2 R_DAGGT_LLA_ther
                         + 25 R_GRTAL_ther - R_NDPK2_ther + 9.8 R_GALTAL_ther
                         + R_URIDK1_ther + R_UGT1_LLA_ther + 25 R_GLTAL_ther
                         + 2 R_CPSS_LLA_ther + R_UGT2_LLA_ther  = 0
 M_udpg_c_bulg:          - 25 R_GTAglcS_bulg - 25 R_GRTAL_bulg - R_UDPG4E_bulg
                         - 2 R_DAGGT_LPL_bulg - 2 R_CPSS_LPL_bulg + R_GALU_bulg
                          = 0
 M_udpg_c_ther:          - R_UDPG4E_ther + R_GALU_ther - 2 R_DAGGT_LLA_ther
                         - 25 R_GRTAL_ther - R_UGT1_LLA_ther - 25 R_GLTAL_ther
                         - 2 R_CPSS_LLA_ther - R_UGT2_LLA_ther  = 0
 M_udpgal_c_bulg:        R_UDPG4E_bulg - 2 R_CPSS_LPL_bulg - R_UDPGALM_bulg
                          = 0
 M_udpgal_c_ther:        R_UDPG4E_ther + R_GALT_ther - R_UDPGALM_ther
                         - 9.8 R_GALTAL_ther - R_CPSS_LLA_ther  = 0
 M_udpgalfur_c_bulg:     R_UDPGALM_bulg  = 0
 M_udpgalfur_c_ther:     R_UDPGALM_ther  = 0
 M_ugmd_c_bulg:          - R_UGMDDS2_bulg + R_UAAGDS_bulg  = 0
 M_ugmd_c_ther:          R_UAAGDS_ther  = 0
 M_ugmdalac_c_bulg:      R_UGMDDS2_bulg - R_PAPPT5_bulg  = 0
 M_ump_c_bulg:           R_UPPRT_bulg + R_OMPDC_bulg + R_PAPPT5_bulg
                         + R_ACGAMT_bulg + R_CPSS_LPL_bulg + R_URIK2_bulg
                         + R_URIK1_bulg + R_URIK3_bulg - R_URIDK1_bulg  = 0
 M_ump_c_ther:           R_URIK1_ther + R_URIK3_ther + R_USHD_ther
                         + R_OMPDC_ther + R_PAPPT1_ther + R_ACGAMT_ther
                         - R_NTD2_ther + R_URIK2_ther - R_URIDK1_ther
                         + R_CPSS_LLA_ther + R_UPPRT_ther  = 0
 M_unGTAglc_c_bulg:      R_GTAglcS_bulg - R_GTAPGL_bulg  = 0
 M_unaga_c_bulg:         R_ACGAMT_bulg - R_UNAGAMAMT_bulg  = 0
 M_unaga_c_ther:         R_ACGAMT_ther  = 0
 M_unagama_c_bulg:       - R_CGPT_bulg + R_UNAGAMAMT_bulg  = 0
 M_unagamagp_c_bulg:     - R_GTAglcS_bulg + R_CGPT_bulg  = 0
 M_ura_c_bulg:           R_URAt2_bulg - R_UPPRT_bulg - R_YUMPS_bulg
                         + R_CSND_bulg + R_PYRNS1_bulg + R_DURIPP_bulg  = 0
 M_ura_c_ther:           R_CSND_ther + R_DURIPP_ther + R_URAt2_ther
                         - R_YUMPS_ther - R_UPPRT_ther  = 0
 M_ura_e:                - R_URAt2_bulg - R_EX_ura_e_ - R_URAt2_ther  = 0
 M_urcan_c_ther:         R_HISD1_ther - R_URCN_ther  = 0
 M_urea_c_ther:          R_AGMT_ther + R_UREAt_ther - R_UREA2_ther  = 0
 M_urea_e:               - R_EX_urea_e_ - R_UREAt_ther  = 0
 M_uri_c_bulg:           R_CYTD_bulg - R_PYRNS1_bulg - R_URIK2_bulg
                         - R_URIK1_bulg - R_URIK3_bulg  = 0
 M_uri_c_ther:           - R_URIK1_ther - R_URIK3_ther + R_NTD2_ther
                         + R_CYTD_ther - R_URIK2_ther  = 0
 M_utp_c_bulg:           - 0.22 R_RNAS_LPL_bulg - R_CTPS1_1_bulg - R_UAGDP_bulg
                         - R_CTPS2_bulg + R_NDPK2_bulg - R_RNTR4_bulg
                         - R_GALU_bulg  = 0
 M_utp_c_ther:           - R_GALT_ther - R_CTPS2_ther - R_GALU_ther
                         - R_RNTR4_ther + R_NDPK2_ther - 0.22 R_RNAS_LLA_ther
                         - R_UAGDP_ther - R_CTPS1_1_ther  = 0
 M_val_L_c_bulg:         - R_VALTRS_bulg - R_VALTA_bulg + M_val_L_c_bulg_src
                         + R_VALt6_bulg + 0.110536255502338 R_CASL_bulg
                         + 0.055 R_PEPTL_bulg  = 0
 M_val_L_c_ther:         R_VALabc_ther + R_VALt6_ther - R_VALTRS_ther
                         - R_VALTA_ther + 0.110536255502338 R_CASL_ther
                         + 0.055 R_PEPTL_ther  = 0
 M_val_L_e:              - R_VALabc_ther - R_VALt6_ther - R_VALt6_bulg
                         - R_EX_val_L_e_  = 0
 M_valtrna_c_bulg:       R_VALTRS_bulg - 0.06 R_PROTS_LPL_v60_bulg  = 0
 M_valtrna_c_ther:       R_VALTRS_ther - 0.052 R_PROTS_STR_ther  = 0
 M_xan_c_bulg:           R_PNS4_bulg + R_XANt2_bulg - R_XPPT_bulg
                         + R_PUNP7_bulg  = 0
 M_xan_c_ther:           R_PUNP7_ther - R_XPPT_ther  = 0
 M_xan_e:                - R_EX_xan_e_ - R_XANt2_bulg  = 0
 M_xmp_c_bulg:           - R_GMPS2_bulg + R_XPPT_bulg + R_IMPD_bulg  = 0
 M_xmp_c_ther:           - R_NTD10_ther + R_XPPT_ther - R_GMPS2_ther
                         + R_IMPD_ther  = 0
 M_xtsn_c_bulg:          - R_PNS4_bulg - R_PUNP7_bulg  = 0
 M_xtsn_c_ther:          R_NTD10_ther - R_PUNP7_ther  = 0
 M_xu5p_D_c_bulg:        - R_PKL_bulg + R_RPE_bulg  = 0
 M_xu5p_D_c_ther:        - R_TKT2_ther - R_TKT1_ther + R_RPE_ther  = 0
 M_zn2_c_bulg:           R_ZN2abc2_bulg  = 0
 M_zn2_c_ther:           R_ZN2abc2_ther  = 0
 M_zn2_e:                - R_ZN2abc2_ther - R_ZN2abc2_bulg  = 0
 ind_R_ASPt6_bulg:      bin_R_ASPt6_bulg = 0 -> R_ASPt6_bulg  = 0
 ind_R_CYSt6_bulg:      bin_R_CYSt6_bulg = 0 -> R_CYSt6_bulg  = 0
 ind_R_TYRt6_bulg:      bin_R_TYRt6_bulg = 0 -> R_TYRt6_bulg  = 0
 ind_R_LYSt6_bulg:      bin_R_LYSt6_bulg = 0 -> R_LYSt6_bulg  = 0
 ind_R_GLUt6_bulg:      bin_R_GLUt6_bulg = 0 -> R_GLUt6_bulg  = 0
 ind_R_GLNabc_bulg:     bin_R_GLNabc_bulg = 0 -> R_GLNabc_bulg  = 0
 ind_R_ARGt_bulg:       bin_R_ARGt_bulg = 0 -> R_ARGt_bulg  = 0
 ind_R_HISt6_bulg:      bin_R_HISt6_bulg = 0 -> R_HISt6_bulg  = 0
 ind_R_PHEt6_bulg:      bin_R_PHEt6_bulg = 0 -> R_PHEt6_bulg  = 0
 ind_R_GLYt6_bulg:      bin_R_GLYt6_bulg = 0 -> R_GLYt6_bulg  = 0
 ind_R_VALt6_bulg:      bin_R_VALt6_bulg = 0 -> R_VALt6_bulg  = 0
 ind_R_THRt6_bulg:      bin_R_THRt6_bulg = 0 -> R_THRt6_bulg  = 0
 ind_R_ALA_Lt6_bulg:    bin_R_ALA_Lt6_bulg = 0 -> R_ALA_Lt6_bulg  = 0
 ind_R_SERt6_bulg:      bin_R_SERt6_bulg = 0 -> R_SERt6_bulg  = 0
 ind_R_ILEt6_bulg:      bin_R_ILEt6_bulg = 0 -> R_ILEt6_bulg  = 0
 ind_R_ARGabc_bulg:     bin_R_ARGabc_bulg = 0 -> R_ARGabc_bulg  = 0
 ind_R_TRPt6_bulg:      bin_R_TRPt6_bulg = 0 -> R_TRPt6_bulg  = 0
 ind_R_ASNt6_bulg:      bin_R_ASNt6_bulg = 0 -> R_ASNt6_bulg  = 0
 ind_R_PROt6_bulg:      bin_R_PROt6_bulg = 0 -> R_PROt6_bulg  = 0
 ind_R_LEUt6_bulg:      bin_R_LEUt6_bulg = 0 -> R_LEUt6_bulg  = 0
 ind_R_gln_export_bulg: bin_R_gln_export_bulg = 0 -> R_gln_export_bulg  = 0
 ind_R_METt6_bulg:      bin_R_METt6_bulg = 0 -> R_METt6_bulg  = 0
 ind_R_VALabc_ther:     bin_R_VALabc_ther = 0 -> R_VALabc_ther  = 0
 ind_R_LEUt6_ther:      bin_R_LEUt6_ther = 0 -> R_LEUt6_ther  = 0
 ind_R_SERt6_ther:      bin_R_SERt6_ther = 0 -> R_SERt6_ther  = 0
 ind_R_VALt6_ther:      bin_R_VALt6_ther = 0 -> R_VALt6_ther  = 0
 ind_R_METabc_ther:     bin_R_METabc_ther = 0 -> R_METabc_ther  = 0
 ind_R_GLNabc_ther:     bin_R_GLNabc_ther = 0 -> R_GLNabc_ther  = 0
 ind_R_PROt6_ther:      bin_R_PROt6_ther = 0 -> R_PROt6_ther  = 0
 ind_R_ALA_Lt6_ther:    bin_R_ALA_Lt6_ther = 0 -> R_ALA_Lt6_ther  = 0
 ind_R_TRPt6_ther:      bin_R_TRPt6_ther = 0 -> R_TRPt6_ther  = 0
 ind_R_GLYt6_ther:      bin_R_GLYt6_ther = 0 -> R_GLYt6_ther  = 0
 ind_R_GLUt6_ther:      bin_R_GLUt6_ther = 0 -> R_GLUt6_ther  = 0
 ind_R_ASPabc_ther:     bin_R_ASPabc_ther = 0 -> R_ASPabc_ther  = 0
 ind_R_GLUabc_ther:     bin_R_GLUabc_ther = 0 -> R_GLUabc_ther  = 0
 ind_R_LYSt6_ther:      bin_R_LYSt6_ther = 0 -> R_LYSt6_ther  = 0
 ind_R_PHEt6_ther:      bin_R_PHEt6_ther = 0 -> R_PHEt6_ther  = 0
 ind_R_CYSt6_ther:      bin_R_CYSt6_ther = 0 -> R_CYSt6_ther  = 0
 ind_R_HISt6_ther:      bin_R_HISt6_ther = 0 -> R_HISt6_ther  = 0
 ind_R_LEUabc_ther:     bin_R_LEUabc_ther = 0 -> R_LEUabc_ther  = 0
 ind_R_TYRt6_ther:      bin_R_TYRt6_ther = 0 -> R_TYRt6_ther  = 0
 ind_R_ARGabc_ther:     bin_R_ARGabc_ther = 0 -> R_ARGabc_ther  = 0
 ind_R_ILEt6_ther:      bin_R_ILEt6_ther = 0 -> R_ILEt6_ther  = 0
 ind_R_THRt6_ther:      bin_R_THRt6_ther = 0 -> R_THRt6_ther  = 0
 ind_R_ARGt_ther:       bin_R_ARGt_ther = 0 -> R_ARGt_ther  = 0
 ind_R_ILEabc_ther:     bin_R_ILEabc_ther = 0 -> R_ILEabc_ther  = 0
 ind_R_METt6_ther:      bin_R_METt6_ther = 0 -> R_METt6_ther  = 0
Bounds
      R_FABM2_ther Free
      R_OCBT_ther Free
      R_ORPT_ther Free
      R_ICDHx_ther Free
      R_ADPTA_ther Free
      R_SDPTA_ther Free
      R_TPI_ther Free
      R_PUNP2_ther Free
      R_ENO_bulg Free
      R_ACONT_ther Free
      R_DIACTt_ther Free
      R_PPM_ther Free
      R_UDPG4E_ther Free
      R_ASPt6_bulg Free
      R_GK2_bulg Free
      R_METFR_bulg Free
      R_GALT_ther Free
      R_MTHFC_ther Free
      R_ASPTA5_ther Free
      R_L_LACt2_bulg Free
      R_H2Ot5_bulg Free
      R_PTPAT_ther Free
      R_ORNTAC_ther Free
      R_NH4t_ther Free
      R_GARFT_bulg Free
      R_ALALAC_bulg Free
      R_ACKr_ther Free
      R_TKT2_ther Free
 3.8103300794 <= R_EX_lys_L_e_ <= 7.07632729032
      R_GLUR_bulg Free
      R_THFGLUS_bulg Free
      R_NDPK1_bulg Free
      R_TPI_bulg Free
      R_TDPDRR_ther Free
      R_SHK3D_ther Free
      R_ADK2_ther Free
      R_FRDx_bulg Free
      R_HMGCOAR_bulg Free
      R_EX_xan_e_ >= -1000
      R_TYRt6_bulg Free
      R_VALt6_ther Free
      R_ACALD_bulg Free
      R_DAPDC_bulg Free
      R_CYTK1_ther Free
      R_EX_gua_e_ >= -1000
      R_EX_caspep_e_ = 0
      R_MTHFD_bulg Free
      R_LYSt6_bulg Free
      R_CYTK1_bulg Free
      R_ALCD19_D_bulg Free
-47.7440740575 <= R_EX_succ_e_ <= -25.7083475694
      R_HSAT_bulg Free
      R_GAPD_bulg Free
      R_G6PI_ther Free
      R_LCTSt6_bulg Free
      R_PGM_ther Free
      R_ACACT1r_bulg Free
      R_SERt6_ther Free
      R_PUNP3_bulg Free
      R_GLUDy_ther Free
      R_AHSERL3_ther Free
      R_ASPTA6_ther Free
      R_LCTSt6_ther = 0
      R_PMANM_ther Free
      R_PRASCS_bulg Free
      R_FRDx_ther Free
      R_ADK1_ther Free
      R_NDPK3_bulg Free
 1.83195065934 <= R_EX_ile_L_e_ <= 3.40219408162
      R_NNAMr_ther Free
      R_HSDy_ther Free
      R_PMANM_bulg Free
      R_ILEt6_ther Free
      R_PUNP7_ther Free
      R_EX_cas = 0
      R_ALCD19_L_bulg Free
      R_ATPS3r_ther Free
      R_HMGCOAS_ther Free
      R_AICART_bulg Free
      R_METFR_ther Free
      R_HSDy_bulg Free
      R_DHORD6_ther Free
      R_GLUt6_bulg Free
      R_2HXICDH_bulg Free
      R_LDH_L_ther Free
      R_ENO_ther Free
      R_MALt6_bulg Free
      R_GK1_bulg Free
      R_VALTA_ther Free
      R_PUNP5_bulg Free
      R_ARGt_bulg Free
      R_PUNP4_bulg Free
      R_EX_hxan_e_ >= -1000
      R_GK2_ther Free
      R_YUMPS_bulg Free
      R_PGMT_bulg Free
      R_MALt6_ther Free
      R_HISt6_bulg Free
      R_FABM2_bulg Free
      R_ASAD_bulg Free
      R_CYTK2_ther Free
      R_ASNt6_bulg Free
 7.84256470494 <= R_EX_leu_L_e_ <= 14.5647630235
      R_TDPDRR_bulg Free
      R_FOLt_bulg Free
      R_PUNP5_ther Free
      R_PRAI_ther Free
      R_NH4DIS_bulg Free
      R_G3PD1_ther Free
      R_ALAR_bulg Free
      R_ACt6_bulg Free
      R_DHFOR2_ther Free
 4.41093922626 <= R_EX_phe_L_e_ <= 8.19174427734
      R_DHORTS_ther Free
      R_TKT1_ther Free
      R_PUNP3_ther Free
      R_IPDDI_bulg Free
      R_NDPK6_bulg Free
      R_GALt6_bulg Free
      R_ATPS3r_bulg Free
      R_CLPNS_LPL_bulg Free
      R_CYSt6_ther Free
 7.14376934397 <= R_EX_trp_L_e_ <= 13.2670002102
      R_CYSt6_bulg Free
      R_TDPDRE_bulg Free
 282.870145714 <= R_biomass_LBUL10_bulg <= 525.330270612
      R_NDPK8_bulg Free
      R_UDPG4E_bulg Free
      R_PRPPS_bulg Free
      R_ADSL1r_ther Free
      R_HCO3E_bulg Free
      R_ILETA_ther Free
      R_EX_pnto_R_e_ >= -1000
      R_LCTSGALex_ther Free
 647.615305142 <= R_EX_for_e_ <= 1202.71413812
      R_ALAR_ther Free
      R_NDPK8_ther Free
      R_PTAr_bulg Free
      R_TRPt6_ther Free
      R_AGPR_ther Free
      R_VALTA_bulg Free
      R_TDPDRE_ther Free
      R_PROt6_bulg Free
      R_IPDDI_ther Free
      R_PTAr_ther Free
      R_GAPD_ther Free
      R_KARA2_ther Free
      R_AKGTA3_bulg Free
      R_PHEt6_bulg Free
      R_NH4DISex_bulg Free
      R_PSCVT_ther Free
      R_THRt6_ther Free
      R_DHORD1_ther Free
      R_GALU_ther Free
      R_METACH_ther Free
      R_CO2t_ther Free
      R_ACALDt_bulg Free
      R_GLNTAL_bulg Free
      R_NDPK4_ther Free
      R_DHAt_ther Free
      R_ADK2_bulg Free
      R_GLYt6_bulg Free
      R_FORt_bulg Free
      R_RPE_ther Free
      R_CITt6_ther Free
      R_ACKr_bulg Free
      R_DASYN_LLA_ther Free
      R_PFL_ther Free
      R_HISt6_ther Free
      R_SERAT_ther Free
      R_ORPT_bulg Free
      R_PHETA1_bulg Free
      R_EX_btn_e_ >= -1000
      R_EX_h2o_e_ Free
 5.24737908822 <= R_EX_met_L_e_ <= 9.7451325924
      R_NDPK2_ther Free
      R_ASPTA1_ther Free
      R_2HXICDH_ther Free
      R_PUNP4_ther Free
      R_THFGLUS_ther Free
      R_SERt6_bulg Free
      R_PUNP1_bulg Free
      R_PGAMT_ther Free
      R_G3PD1_bulg Free
      R_ALAD_Lr_ther Free
      R_VALt6_bulg Free
      R_DHFR_bulg Free
      R_DURIPP_ther Free
      R_EX_pi_e_ Free
      R_LEUTA_bulg Free
      R_PYDAMt_bulg Free
      R_PGI_bulg Free
      R_ALCD19_ther Free
      R_EX_pydam_e_ >= -1000
      R_NDPK7_ther Free
      R_NDPK3_ther Free
      R_ALA_Lt6_ther Free
      R_SUCCt6_bulg Free
      R_PRPPS_ther Free
      R_UREAt_ther Free
      R_O2t_ther Free
 576.164470461 <= R_biomass_STR_ther <= 1070.01973086
      R_ASPTA4_ther Free
      R_DHFR_ther Free
      R_ALCD2x_ther Free
      R_TYRTA_bulg Free
      R_G5SADs_ther Free
      R_NDPK9_bulg Free
      R_CPSS_LPL_bulg Free
      R_ADK1_bulg Free
      R_RPI_bulg Free
      R_SUCCt6_ther Free
      R_G6PI_bulg Free
      R_GLUt6_ther Free
      R_DAPE_ther Free
      R_PTPAT_bulg Free
      R_PGK_bulg Free
      R_HSAT_ther Free
      R_PUNP2_bulg Free
      R_EX_glc_e_ = 0
      R_OXFOt4_ther Free
      R_NDPK9_ther Free
      R_DHQD_ther Free
      R_FOLt_ther Free
      R_NDPK5_bulg Free
      R_DASYN_LPL_bulg Free
      R_ASPK_bulg Free
      R_DHORD6_bulg Free
      R_PGK_ther Free
      R_GTHRD_ther Free
      R_NNAMr_bulg Free
      R_FABM1_bulg Free
      R_GLUR_ther Free
      R_G6PDHy_bulg Free
      R_LYSt6_ther Free
      R_KARA1_ther Free
      R_SUCD3BC_bulg Free
      R_NDPK1_ther Free
      R_FBA_ther Free
      R_ALAALA_ther Free
      R_GARFT_ther Free
      R_PHEt6_ther Free
      R_AKGTA3_ther Free
      R_GCALDt_ther Free
      R_NDPK6_ther Free
      R_LDH_D_bulg Free
      R_ETOHt1_ther Free
      R_NDPK7_bulg Free
      R_FBA_bulg Free
      R_ECOAH2c_ther Free
      R_MTHFC_bulg Free
      R_THRt6_bulg Free
      R_CLPNS_LLA_ther Free
      R_NDPK2_bulg Free
      R_PKL_bulg Free
      R_ACt6_ther Free
      R_RPI_ther Free
      R_PUNP6_ther Free
      R_DHORTS_bulg Free
      R_PGMT_ther Free
      R_ARGSL_ther Free
      R_ALA_Lt6_bulg Free
      R_PGSA_LLA_ther Free
      R_ICDHy_ther = 0
      R_PGI_ther Free
      R_EX_h_e_ Free
-6357.79080065 <= R_EX_lcts_e_ <= -3423.42581573
      R_AHSERL4_ther Free
      R_NH3t_bulg Free
      R_YUMPS_ther Free
      R_DTMPK_bulg Free
      R_MAN6PI_ther Free
      R_CO2t_bulg Free
      R_FORt2_ther Free
      R_AICART_ther Free
      R_GALM_ther Free
      R_GALU_bulg Free
      R_L_LACt2_ther Free
      R_EX_nac_e_ >= -1000
      R_PUNP6_bulg Free
      R_EX_fru_e_ = 0
      R_TRPt6_bulg Free
 7910.0227157 <= R_EX_lac_L_e_ <= 14690.0421863
      R_EX_nh4_e_ Free
      R_THRA_ther Free
      R_HMGCOAR_ther Free
      R_APAT_ther Free
      R_TYRt6_ther Free
      R_EX_o2_e_ = 0
      R_PNTOt2_ther Free
      R_MTHFD_ther Free
      R_PPM2_ther Free
      R_LEUt6_ther Free
      R_ILETA_bulg Free
      R_PYDAMt_ther Free
      R_PUNP8_ther Free
      R_ASPK_ther Free
      R_P5CRr_ther Free
      R_PRASCS_ther Free
      R_PUNP7_bulg Free
      R_HMGCOAS_bulg Free
      R_LEUt6_bulg Free
      R_MAN6PI_bulg Free
      R_FABM1_ther Free
      R_PGSA_LPL_bulg Free
      R_ASAD_ther Free
      R_PROt6_ther Free
      R_LDH_L_bulg Free
      R_Kt2_ther Free
      R_IMPC_ther Free
      R_ALCD2x_bulg Free
      R_HCO3E_ther Free
      R_DHFOR2_bulg Free
      R_EX_ribflv_e_ >= -1000
 5834.79479723 <= R_EX_gal_e_ <= 10836.0474806
      R_IPPMIb_ther Free
      R_PRAI_bulg Free
      R_METt6_bulg Free
      R_GLYt6_ther Free
      R_DTMPK_ther Free
      R_CYTK2_bulg Free
      R_FABM_ther Free
      R_NDPK5_ther Free
      R_EX_thm_e_ >= -1000
      R_DURIPP_bulg Free
      R_ETOHt1_bulg Free
      R_H2Ot5_ther Free
      R_ALATA_L_ther Free
      R_GLYCt5_ther Free
      R_GHMT_ther Free
      R_ACOTA_ther Free
      R_GK1_ther Free
      R_LEUTA_ther Free
      R_PUNP1_ther Free
 6.07920095975 <= R_EX_val_L_e_ <= 11.2899446395
      R_PGM_bulg Free
 1.28288790489 <= R_EX_gln_L_e_ <= 2.38250610907
      R_IPPMIa_ther Free
      R_PNTOt2_bulg Free
      R_PGMT_B_ther Free
      R_ARGt_ther Free
      R_CPSS_LLA_ther Free
      R_GLNTAL_ther Free
      R_ASPTA2_ther Free
      R_EX_mpept_e_ = 0
      R_ASPR_bulg Free
      R_UAG2E_bulg Free
 -Inf <= R_gln_export_bulg <= 0
      R_EX_co2_e_ Free
      R_ADSL2r_ther Free
      R_ILEt6_bulg Free
      R_FUM_bulg Free
      R_ALATA_Lr_bulg Free
      R_RPE_bulg Free
 4.37879246716 <= R_EX_ser_L_e_ <= 8.1320431533
      R_IMPC_bulg Free
      R_NDPK4_bulg Free
      R_EX_mprot_e_ = 0
      R_PGAMT_bulg Free
      R_METt6_ther Free
 0 <= bin_R_ASPt6_bulg <= 1
 0 <= bin_R_CYSt6_bulg <= 1
 0 <= bin_R_TYRt6_bulg <= 1
 0 <= bin_R_LYSt6_bulg <= 1
 0 <= bin_R_GLUt6_bulg <= 1
 0 <= bin_R_GLNabc_bulg <= 1
 0 <= bin_R_ARGt_bulg <= 1
 0 <= bin_R_HISt6_bulg <= 1
 0 <= bin_R_PHEt6_bulg <= 1
 0 <= bin_R_GLYt6_bulg <= 1
 0 <= bin_R_VALt6_bulg <= 1
 0 <= bin_R_THRt6_bulg <= 1
 0 <= bin_R_ALA_Lt6_bulg <= 1
 0 <= bin_R_SERt6_bulg <= 1
 0 <= bin_R_ILEt6_bulg <= 1
 0 <= bin_R_ARGabc_bulg <= 1
 0 <= bin_R_TRPt6_bulg <= 1
 0 <= bin_R_ASNt6_bulg <= 1
 0 <= bin_R_PROt6_bulg <= 1
 0 <= bin_R_LEUt6_bulg <= 1
 0 <= bin_R_gln_export_bulg <= 1
 0 <= bin_R_METt6_bulg <= 1
 0 <= bin_R_VALabc_ther <= 1
 0 <= bin_R_LEUt6_ther <= 1
 0 <= bin_R_SERt6_ther <= 1
 0 <= bin_R_VALt6_ther <= 1
 0 <= bin_R_METabc_ther <= 1
 0 <= bin_R_GLNabc_ther <= 1
 0 <= bin_R_PROt6_ther <= 1
 0 <= bin_R_ALA_Lt6_ther <= 1
 0 <= bin_R_TRPt6_ther <= 1
 0 <= bin_R_GLYt6_ther <= 1
 0 <= bin_R_GLUt6_ther <= 1
 0 <= bin_R_ASPabc_ther <= 1
 0 <= bin_R_GLUabc_ther <= 1
 0 <= bin_R_LYSt6_ther <= 1
 0 <= bin_R_PHEt6_ther <= 1
 0 <= bin_R_CYSt6_ther <= 1
 0 <= bin_R_HISt6_ther <= 1
 0 <= bin_R_LEUabc_ther <= 1
 0 <= bin_R_TYRt6_ther <= 1
 0 <= bin_R_ARGabc_ther <= 1
 0 <= bin_R_ILEt6_ther <= 1
 0 <= bin_R_THRt6_ther <= 1
 0 <= bin_R_ARGt_ther <= 1
 0 <= bin_R_ILEabc_ther <= 1
 0 <= bin_R_METt6_ther <= 1
Binaries
 bin_R_ASPt6_bulg  bin_R_CYSt6_bulg  bin_R_TYRt6_bulg  bin_R_LYSt6_bulg 
 bin_R_GLUt6_bulg  bin_R_GLNabc_bulg  bin_R_ARGt_bulg  bin_R_HISt6_bulg 
 bin_R_PHEt6_bulg  bin_R_GLYt6_bulg  bin_R_VALt6_bulg  bin_R_THRt6_bulg 
 bin_R_ALA_Lt6_bulg  bin_R_SERt6_bulg  bin_R_ILEt6_bulg  bin_R_ARGabc_bulg 
 bin_R_TRPt6_bulg  bin_R_ASNt6_bulg  bin_R_PROt6_bulg  bin_R_LEUt6_bulg 
 bin_R_gln_export_bulg  bin_R_METt6_bulg  bin_R_VALabc_ther  bin_R_LEUt6_ther 
 bin_R_SERt6_ther  bin_R_VALt6_ther  bin_R_METabc_ther  bin_R_GLNabc_ther 
 bin_R_PROt6_ther  bin_R_ALA_Lt6_ther  bin_R_TRPt6_ther  bin_R_GLYt6_ther 
 bin_R_GLUt6_ther  bin_R_ASPabc_ther  bin_R_GLUabc_ther  bin_R_LYSt6_ther 
 bin_R_PHEt6_ther  bin_R_CYSt6_ther  bin_R_HISt6_ther  bin_R_LEUabc_ther 
 bin_R_TYRt6_ther  bin_R_ARGabc_ther  bin_R_ILEt6_ther  bin_R_THRt6_ther 
 bin_R_ARGt_ther  bin_R_ILEabc_ther  bin_R_METt6_ther 
End
