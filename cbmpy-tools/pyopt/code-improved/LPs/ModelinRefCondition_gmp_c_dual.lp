\ENCODING=ISO-8859-1
\Problem name: 1471131711.mip

Minimize
 rhs: - 0.0523306361177878 l1 + 37 l832 + 300 l873
Subject To
 dR_BiomassAuto:        31.84006918 pM_adp_c - 31.84006918 pM_atp_c
                        - 0.100827 pM_biomass_dna_c
                        - 0.21031 pM_biomass_glycogen_c
                        - 0.2327419 pM_biomass_ion_c
                        - 0.15719526 pM_biomass_lipid_c
                        - 0.025228 pM_biomass_lps_c
                        - 0.025228 pM_biomass_pept_p
                        - 0.0109472 pM_biomass_pigment_c
                        - 4.600079 pM_biomass_protein_c
                        - 0.530238 pM_biomass_rna_c
                        - 0.06843016 pM_biomass_solublep_c
                        - 27.23996718 pM_h2o_c + 31.84006918 pM_h_c
                        + 31.83248018 pM_pi_c - l1  = 1
 dR_FMETTRS:            - pM_10fthf_c + pM_fmettrna_c + pM_h_c - pM_mettrna_c
                        + pM_thf_c >= 0
 dR_BiomassSolubleP:    - 0.003258797 pM_10fthf_c - 0.003258797 pM_5mthf_c
                        - 0.0040771496 pM_accoa_c - 0.003258797 pM_adocbl_c
                        - 0.003258797 pM_amet_c + pM_biomass_solublep_c
                        - 0.3096733955 pM_cholphya_c - 0.003258797 pM_chor_c
                        - 0.0024477511 pM_coa_c - 0.003263181 pM_fad_c
                        - 0.003258797 pM_gthrd_c - 0.003258797 pM_hemeA_c
                        - 0.003258797 pM_hemeO_c - 0.0004569184 pM_malcoa_c
                        - 0.003263181 pM_mlthf_c - 0.0261098323 pM_nad_c
                        - 0.0006527385 pM_nadh_c - 0.0016318828 pM_nadp_c
                        - 0.0048955022 pM_nadph_c - 0.0032636194 pM_pheme_c
                        - 0.0240303398 pM_phllqne_c - 0.4861891306 pM_ptrc_c
                        - 0.0032636194 pM_ribflv_c - 0.0985574197 pM_spmd_c
                        - 0.0014360189 pM_succoa_c - 0.0032636194 pM_thf_c
                        - 0.0007822574 pM_udcpdp_c >= 0
 dR_AICART:             - pM_10fthf_c - pM_aicar_c + pM_fprica_c + pM_thf_c
                         = 0
 dR_FTHFD:              - pM_10fthf_c + pM_for_c - pM_h2o_c + pM_h_c + pM_thf_c
                        >= 0
 dR_MTHFC:              pM_10fthf_c - pM_h2o_c + pM_h_c - pM_methf_c  = 0
 dR_FTHFLi:             pM_10fthf_c + pM_adp_c - pM_atp_c - pM_for_c + pM_pi_c
                        - pM_thf_c >= 0
 dR_GARFT:              - pM_10fthf_c + pM_fgam_c - pM_gar_c + pM_h_c
                        + pM_thf_c  = 0
 dR_MGDG160:            - pM_12dgr160_c + pM_h_c + pM_mgdg160_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_SQD2_160:           - pM_12dgr160_c - pM_h2o_c + pM_h_c + pM_sqdg160_c
                        + pM_udp_c - pM_udpsq_c >= 0
 dR_PAPA160:            pM_12dgr160_c - pM_h2o_c - pM_pa160_c + pM_pi_c >= 0
 dR_PAPA161:            pM_12dgr161_c - pM_h2o_c - pM_pa161_c + pM_pi_c >= 0
 dR_MGDG161:            - pM_12dgr161_c + pM_h_c + pM_mgdg161_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_SQD2_161:           - pM_12dgr161_c - pM_h2o_c + pM_h_c + pM_sqdg161_c
                        + pM_udp_c - pM_udpsq_c >= 0
 dR_MGDG180:            - pM_12dgr180_c + pM_h_c + pM_mgdg180_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_SQD2_180:           - pM_12dgr180_c - pM_h2o_c + pM_h_c + pM_sqdg180_c
                        + pM_udp_c - pM_udpsq_c >= 0
 dR_PAPA180:            pM_12dgr180_c - pM_h2o_c - pM_pa180_c + pM_pi_c >= 0
 dR_PAPA181_9:          pM_12dgr181_9_c - pM_h2o_c - pM_pa181_9_c + pM_pi_c
                        >= 0
 dR_MGDG181_9:          - pM_12dgr181_9_c + pM_h_c + pM_mgdg181_9_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_SQD2_181_9:         - pM_12dgr181_9_c - pM_h2o_c + pM_h_c + pM_sqdg181_9_c
                        + pM_udp_c - pM_udpsq_c >= 0
 dR_MGDG181:            - pM_12dgr181_c + pM_h_c + pM_mgdg181_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_SQD2_181:           - pM_12dgr181_c - pM_h2o_c + pM_h_c + pM_sqdg181_c
                        + pM_udp_c - pM_udpsq_c >= 0
 dR_PAPA181:            pM_12dgr181_c - pM_h2o_c - pM_pa181_c + pM_pi_c >= 0
 dR_MGDG182_9_12:       - pM_12dgr182_9_12_c + pM_h_c + pM_mgdg182_9_12_c
                        + pM_udp_c - pM_udpgal_c >= 0
 dR_SQD2_182_9_12:      - pM_12dgr182_9_12_c - pM_h2o_c + pM_h_c
                        + pM_sqdg182_9_12_c + pM_udp_c - pM_udpsq_c >= 0
 dR_PAPA182_9_12:       pM_12dgr182_9_12_c - pM_h2o_c - pM_pa182_9_12_c
                        + pM_pi_c >= 0
 dR_PAPA183_6_9_12:     pM_12dgr183_6_9_12_c - pM_h2o_c - pM_pa183_6_9_12_c
                        + pM_pi_c >= 0
 dR_MGDG183_6_9_12:     - pM_12dgr183_6_9_12_c + pM_h_c + pM_mgdg183_6_9_12_c
                        + pM_udp_c - pM_udpgal_c >= 0
 dR_SQD2_183_6_9_12:    - pM_12dgr183_6_9_12_c - pM_h2o_c + pM_h_c
                        + pM_sqdg183_6_9_12_c + pM_udp_c - pM_udpsq_c >= 0
 dR_MGDG183_9_12_15:    - pM_12dgr183_9_12_15_c + pM_h_c + pM_mgdg183_9_12_15_c
                        + pM_udp_c - pM_udpgal_c >= 0
 dR_SQD2_183_9_12_15:   - pM_12dgr183_9_12_15_c - pM_h2o_c + pM_h_c
                        + pM_sqdg183_9_12_15_c + pM_udp_c - pM_udpsq_c >= 0
 dR_PAPA183_9_12_15:    pM_12dgr183_9_12_15_c - pM_h2o_c - pM_pa183_9_12_15_c
                        + pM_pi_c >= 0
 dR_PAPA184_6_9_12_15:  pM_12dgr184_6_9_12_15_c - pM_h2o_c
                        - pM_pa184_6_9_12_15_c + pM_pi_c >= 0
 dR_SQD2_184_6_9_12_15: - pM_12dgr184_6_9_12_15_c - pM_h2o_c + pM_h_c
                        + pM_sqdg184_6_9_12_15_c + pM_udp_c - pM_udpsq_c >= 0
 dR_MGDG184_6_9_12_15:  - pM_12dgr184_6_9_12_15_c + pM_h_c
                        + pM_mgdg184_6_9_12_15_c + pM_udp_c - pM_udpgal_c >= 0
 dR_LCARS:              pM_12ppd_S_c - pM_h_c - pM_lald_L_c + pM_nad_c
                        - pM_nadh_c  = 0
 dR_GAPDi_nadp:         - pM_13dpg_c + pM_g3p_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c + pM_pi_c >= 0
 dR_PGK:                pM_13dpg_c - pM_3pg_c + pM_adp_c - pM_atp_c  = 0
 dR_DPGM:               - pM_13dpg_c + pM_23dpg_c + pM_h_c  = 0
 dR_GAPDi:              pM_13dpg_c - pM_g3p_c + pM_h_c - pM_nad_c + pM_nadh_c
                        - pM_pi_c >= 0
 dR_LYSDC:              pM_15dap_c + pM_co2_c - pM_h_c - pM_lys_L_c >= 0
 dR_AGPAT161:           - pM_1hdec9eg3p_c + pM_ACP_c - pM_hdeACP_c + pM_pa161_c
                        >= 0
 dR_G3PAT161:           pM_1hdec9eg3p_c + pM_ACP_c - pM_glyc3p_c - pM_hdeACP_c
                        >= 0
 dR_AGPAT160:           - pM_1hdecg3p_c + pM_ACP_c + pM_pa160_c - pM_palmACP_c
                        >= 0
 dR_G3PAT160:           pM_1hdecg3p_c + pM_ACP_c - pM_glyc3p_c - pM_palmACP_c
                        >= 0
 dR_AGPAT181:           - pM_1odec11eg3p_c + pM_ACP_c - pM_octeACP_c
                        + pM_pa181_c >= 0
 dR_G3PAT181:           pM_1odec11eg3p_c + pM_ACP_c - pM_glyc3p_c
                        - pM_octeACP_c >= 0
 dR_G3PAT183_6_9_12:    pM_1odec691215eg3p_c + pM_ACP_c - pM_glyc3p_c
                        - pM_octe_6_9_12_15_ACP_c >= 0
 dR_AGPAT184_6_9_12_15: - pM_1odec691215eg3p_c + pM_ACP_c
                        - pM_octe_6_9_12_15_ACP_c + pM_pa184_6_9_12_15_c >= 0
 dR_AGPAT183_6_9_12:    - pM_1odec6912eg3p_c + pM_ACP_c - pM_octe_6_9_12_ACP_c
                        + pM_pa183_6_9_12_c >= 0
 dR_G3PAT184_6_9_12_15: pM_1odec6912eg3p_c + pM_ACP_c - pM_glyc3p_c
                        - pM_octe_6_9_12_ACP_c >= 0
 dR_G3PAT183_9_12_15:   pM_1odec91215eg3p_c + pM_ACP_c - pM_glyc3p_c
                        - pM_octe_9_12_15_ACP_c >= 0
 dR_AGPAT183_9_12_15:   - pM_1odec91215eg3p_c + pM_ACP_c
                        - pM_octe_9_12_15_ACP_c + pM_pa183_9_12_15_c >= 0
 dR_AGPAT182_9_12:      - pM_1odec912eg3p_c + pM_ACP_c - pM_octe_9_12_ACP_c
                        + pM_pa182_9_12_c >= 0
 dR_G3PAT182_9_12:      pM_1odec912eg3p_c + pM_ACP_c - pM_glyc3p_c
                        - pM_octe_9_12_ACP_c >= 0
 dR_G3PAT181_9:         pM_1odec9eg3p_c + pM_ACP_c - pM_glyc3p_c
                        - pM_octe_9_ACP_c >= 0
 dR_AGPAT181_9:         - pM_1odec9eg3p_c + pM_ACP_c - pM_octe_9_ACP_c
                        + pM_pa181_9_c >= 0
 dR_AGPAT180:           - pM_1odecg3p_c + pM_ACP_c - pM_ocdcaACP_c + pM_pa180_c
                        >= 0
 dR_G3PAT180:           pM_1odecg3p_c + pM_ACP_c - pM_glyc3p_c - pM_ocdcaACP_c
                        >= 0
 dR_HPROb:              - pM_1p3h5c_c + pM_4hpro_LT_c - 2 pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_HPROa:              - pM_1p3h5c_c + pM_4hpro_LT_c - 2 pM_h_c + pM_nad_c
                        - pM_nadh_c  = 0
 dR_PHCD:               - pM_1p3h5c_c + pM_e4hglu_c - 2 pM_h2o_c + pM_h_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_P5CR:               - pM_1pyr5c_c - 2 pM_h_c + pM_nadp_c - pM_nadph_c
                        + pM_pro_L_c >= 0
 dR_P5CD:               - pM_1pyr5c_c + pM_glu_L_c - 2 pM_h2o_c + pM_h_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_PROD2:              pM_1pyr5c_c - pM_fad_c + pM_fadh2_c + pM_h_c
                        - pM_pro_L_c >= 0
 dR_G5SADs:             pM_1pyr5c_c - pM_glu5sa_c + pM_h2o_c + pM_h_c >= 0
 dR_P5CRx:              - pM_1pyr5c_c - 2 pM_h_c + pM_nad_c - pM_nadh_c
                        + pM_pro_L_c >= 0
 dR_23DBtpp:            pM_23DB_c - pM_23DB_p  = 0
 dR_23BDD:              pM_23DB_c - pM_actn_R_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_EX_23DB_e:          - pM_23DB_e >= 0
 dR_23DBtex:            - pM_23DB_e + pM_23DB_p  = 0
 dR_DHDPRy:             - pM_23dhdp_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        + pM_thdp_c >= 0
 dR_DHDPS:              pM_23dhdp_c - pM_aspsa_c + 2 pM_h2o_c + pM_h_c
                        - pM_pyr_c >= 0
 dR_KARA1i:             pM_23dhmb_c - pM_alac_S_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_DHAD1:              - pM_23dhmb_c + pM_3mob_c + pM_h2o_c >= 0
 dR_DHAD2:              - pM_23dhmp_c + pM_3mop_c + pM_h2o_c >= 0
 dR_KARA2i:             pM_23dhmp_c - pM_2ahbut_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_TOCOPHS1:           - pM_23dmphol_c + pM_gtocophe_c  = 0
 dR_MPBQ:               pM_23dmphol_c - pM_2m6phol_c + pM_ahcys_c - pM_amet_c
                        + pM_h_c  = 0
 dR_ADSL2r:             - pM_25aics_c + pM_aicar_c + pM_fum_c  = 0
 dR_PRASCS:             pM_25aics_c - pM_5aizc_c + pM_adp_c - pM_asp_L_c
                        - pM_atp_c + pM_h_c + pM_pi_c  = 0
 dR_DHPPDA:             - pM_25dhpp_c + pM_5apru_c - pM_h2o_c - pM_h_c
                        + pM_nh4_c >= 0
 dR_GTPCII:             pM_25dhpp_c + pM_for_c - pM_gtp_c - 3 pM_h2o_c
                        + 2 pM_h_c + pM_ppi_c >= 0
 dR_SDPDS:              pM_26dap_LL_c - pM_h2o_c - pM_sl26da_c + pM_succ_c >= 0
 dR_MNXR33992:          pM_26dap_LL_c + pM_akg_c - pM_glu_L_c - pM_h2o_c
                        - pM_h_c - pM_thdp_c >= 0
 dR_DAPE:               - pM_26dap_LL_c + pM_26dap_M_c  = 0
 dR_DAPDC:              - pM_26dap_M_c + pM_co2_c - pM_h_c + pM_lys_L_c >= 0
 dR_UAAGDS:             - pM_26dap_M_c + pM_adp_c - pM_atp_c + pM_h_c + pM_pi_c
                        - pM_uamag_c + pM_ugmd_c >= 0
 dR_ACHBS:              pM_2ahbut_c - pM_2obut_c + pM_co2_c - pM_h_c - pM_pyr_c
                        >= 0
 dR_HPPK:               pM_2ahhmd_c - pM_2ahhmp_c + pM_amp_c - pM_atp_c
                        + pM_h_c >= 0
 dR_FOLD3:              - pM_2ahhmd_c - pM_4abz_c + pM_dhpt_c + pM_ppi_c >= 0
 dR_DHNPA:              pM_2ahhmp_c - pM_dhnpt_c + pM_gcald_c >= 0
 dR_DHPS:               - pM_2ahhmp_c - pM_4abz_c + pM_dhpt_c + pM_h2o_c >= 0
 dR_IGPS:               - pM_2cpr5p_c + pM_3ig3p_c + pM_co2_c + pM_h2o_c
                        - pM_h_c >= 0
 dR_PRAIi:              pM_2cpr5p_c - pM_pran_c >= 0
 dR_DHQS:               - pM_2dda7p_c + pM_3dhq_c + pM_pi_c >= 0
 dR_DDPA:               pM_2dda7p_c - pM_e4p_c - pM_h2o_c - pM_pep_c + pM_pi_c
                        >= 0
 dR_DPR:                - pM_2dhp_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        + pM_pant_R_c >= 0
 dR_MOHMT:              pM_2dhp_c - pM_3mob_c - pM_h2o_c - pM_mlthf_c
                        + pM_thf_c >= 0
 dR_DRPA:               - pM_2dr5p_c + pM_acald_c + pM_g3p_c >= 0
 dR_DRBK:               pM_2dr5p_c + pM_adp_c - pM_atp_c - pM_drib_c + pM_h_c
                        >= 0
 dR_GLXCL:              pM_2h3oppan_c + pM_co2_c - 2 pM_glx_c - pM_h_c >= 0
 dR_TRSARr:             - pM_2h3oppan_c + pM_glyc_R_c - pM_h_c + pM_nad_c
                        - pM_nadh_c  = 0
 dR_PPOR:               pM_2hyoxplac_c + pM_co2_c - pM_o2_c - pM_phpyr_c >= 0
 dR_IPPMIa:             pM_2ippm_c - pM_3c2hmp_c + pM_h2o_c  = 0
 dR_IPPMIb:             - pM_2ippm_c + pM_3c3hmp_c - pM_h2o_c  = 0
 dR_TOCOPHS2:           - pM_2m6phol_c + pM_dtocophe_c >= 0
 dR_HGPHT:              pM_2m6phol_c + pM_co2_c - pM_h_c - pM_hgentis_c
                        - pM_phdp_c + pM_ppi_c >= 0
 dR_TMPPP:              - pM_2mahmp_c - pM_4mpetz_c - pM_h_c + pM_ppi_c
                        + pM_thmmp_c >= 0
 dR_MEPCT:              - pM_2me4p_c + pM_4c2me_c - pM_ctp_c - pM_h_c
                        + pM_ppi_c >= 0
 dR_DXPRIi:             pM_2me4p_c - pM_dxyl5p_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_MECDPDH_syn:        - pM_2mecdp_c + pM_h2mb4p_c + pM_h2o_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_MECDPS:             pM_2mecdp_c - pM_2p4c2me_c + pM_cmp_c >= 0
 dR_THRD_L:             pM_2obut_c + pM_nh4_c - pM_thr_L_c >= 0
 dR_OPHBDC:             pM_2oph_c - pM_3ophb_c + pM_co2_c - pM_h_c >= 0
 dR_CDPMEK:             pM_2p4c2me_c - pM_4c2me_c + pM_adp_c - pM_atp_c
                        + pM_h_c >= 0
 dR_ENO:                - pM_2pg_c + pM_h2o_c + pM_pep_c  = 0
 dR_PGM:                - pM_2pg_c + pM_3pg_c  = 0
 dR_RBCh:               pM_2pglyc_c + pM_3pg_c + 2 pM_h_c - pM_o2_c
                        - pM_rb15bp_c >= 0
 dR_PGLYCP:             - pM_2pglyc_c + pM_glyclt_c - pM_h2o_c + pM_pi_c >= 0
 dR_SUCBZS:             - pM_2shchc_c + pM_h2o_c + pM_sucbz_c >= 0
 dR_SHCHCS2:            pM_2shchc_c - pM_ichor_c + pM_pyr_c - pM_ssaltpp_c
                        + pM_thmpp_c >= 0
 dR_34DHOXPEGOX:        - pM_34dhmald_c + pM_34dhoxpeg_c - pM_h_c + pM_nad_c
                        - pM_nadh_c  = 0
 dR_TYRTA:              pM_34hpp_c - pM_akg_c + pM_glu_L_c - pM_tyr_L_c  = 0
 dR_34HPPOR:            - pM_34hpp_c + pM_co2_c + pM_hgentis_c - pM_o2_c >= 0
 dR_PPND:               pM_34hpp_c + pM_co2_c - pM_nad_c + pM_nadh_c
                        - pM_pphn_c >= 0
 dR_GUACYC:             pM_35cgmp_c - pM_gtp_c + pM_ppi_c >= 0
 dR_IPMD:               - pM_3c2hmp_c + pM_3c4mop_c + pM_h_c - pM_nad_c
                        + pM_nadh_c >= 0
 dR_IPPS:               pM_3c3hmp_c - pM_3mob_c - pM_accoa_c + pM_coa_c
                        - pM_h2o_c + pM_h_c >= 0
 dR_OMCDC:              - pM_3c4mop_c + pM_4mop_c + pM_co2_c - pM_h_c >= 0
 dR_DHQTi:              - pM_3dhq_c + pM_3dhsk_c + pM_h2o_c >= 0
 dR_SHK3D:              - pM_3dhsk_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        + pM_skm_c >= 0
 dR_3HAD40:             - pM_3haACP_c + pM_but2eACP_c + pM_h2o_c >= 0
 dR_3OAR40:             pM_3haACP_c - pM_actACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_PHBS_syn:           - pM_3hbcoa_R_c + pM_PHB_c + pM_coa_c - pM_phbg_c >= 0
 dR_AACOAR_syn:         pM_3hbcoa_R_c - pM_aacoa_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_3OAR121:            pM_3hcddec5eACP_c - pM_3ocddec5eACP_c - pM_h_c
                        + pM_nadp_c - pM_nadph_c >= 0
 dR_3HAD121:            - pM_3hcddec5eACP_c + pM_h2o_c + pM_t3c5ddeceACP_c >= 0
 dR_3OAR141:            pM_3hcmrs7eACP_c - pM_3ocmrs7eACP_c - pM_h_c
                        + pM_nadp_c - pM_nadph_c >= 0
 dR_3HAD141:            - pM_3hcmrs7eACP_c + pM_h2o_c + pM_t3c7mrseACP_c >= 0
 dR_3OAR161:            pM_3hcpalm9eACP_c - pM_3ocpalm9eACP_c - pM_h_c
                        + pM_nadp_c - pM_nadph_c >= 0
 dR_3HAD161:            - pM_3hcpalm9eACP_c + pM_h2o_c + pM_t3c9palmeACP_c >= 0
 dR_3HAD181:            - pM_3hcvac11eACP_c + pM_h2o_c + pM_t3c11vaceACP_c >= 0
 dR_3OAR181:            pM_3hcvac11eACP_c - pM_3ocvac11eACP_c - pM_h_c
                        + pM_nadp_c - pM_nadph_c >= 0
 dR_3HAD120:            - pM_3hddecACP_c + pM_h2o_c + pM_tddec2eACP_c >= 0
 dR_3OAR120:            pM_3hddecACP_c - pM_3oddecACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_3HAD100:            - pM_3hdecACP_c + pM_h2o_c + pM_tdec2eACP_c >= 0
 dR_3OAR100:            pM_3hdecACP_c - pM_3odecACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_3HAD60:             - pM_3hhexACP_c + pM_h2o_c + pM_thex2eACP_c >= 0
 dR_3OAR60:             pM_3hhexACP_c - pM_3ohexACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_HIBD_syn:           - pM_3hmp_c + pM_h_c + pM_mmtsa_c - pM_nad_c
                        + pM_nadh_c >= 0
 dR_3OAR140:            pM_3hmrsACP_c - pM_3omrsACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_UAGAAT:             - pM_3hmrsACP_c + pM_ACP_c + pM_u3aga_c - pM_uacgam_c
                         = 0
 dR_U23GAAT:            - pM_3hmrsACP_c + pM_ACP_c + pM_h_c + pM_u23ga_c
                        - pM_u3hga_c >= 0
 dR_3HAD140:            - pM_3hmrsACP_c + pM_h2o_c + pM_tmrs2eACP_c >= 0
 dR_3OAR80:             pM_3hoctACP_c - pM_3ooctACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_3HAD80:             - pM_3hoctACP_c + pM_h2o_c + pM_toct2eACP_c >= 0
 dR_3OAR180:            pM_3hoctaACP_c - pM_3ooctdACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_3HAD180:            - pM_3hoctaACP_c + pM_h2o_c + pM_toctd2eACP_c >= 0
 dR_3HAD160:            - pM_3hpalmACP_c + pM_h2o_c + pM_tpalm2eACP_c >= 0
 dR_3OAR160:            pM_3hpalmACP_c - pM_3opalmACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c  = 0
 dR_TRPS3:              - pM_3ig3p_c + pM_g3p_c + pM_indole_c >= 0
 dR_TRPS1:              - pM_3ig3p_c + pM_g3p_c + pM_h2o_c - pM_ser_L_c
                        + pM_trp_L_c >= 0
 dR_VALTA:              pM_3mob_c - pM_akg_c + pM_glu_L_c - pM_val_L_c  = 0
 dR_ILETA:              pM_3mop_c - pM_akg_c + pM_glu_L_c - pM_ile_L_c  = 0
 dR_3OAS121:            pM_3ocddec5eACP_c + pM_ACP_c - pM_cdec3eACP_c
                        + pM_co2_c - pM_h_c - pM_malACP_c >= 0
 dR_3OAS141:            pM_3ocmrs7eACP_c + pM_ACP_c - pM_cddec5eACP_c
                        + pM_co2_c - pM_h_c - pM_malACP_c >= 0
 dR_3OAS161:            pM_3ocpalm9eACP_c + pM_ACP_c + pM_co2_c - pM_h_c
                        - pM_malACP_c - pM_tdeACP_c >= 0
 dR_3OAS181:            pM_3ocvac11eACP_c + pM_ACP_c + pM_co2_c - pM_h_c
                        - pM_hdeACP_c - pM_malACP_c >= 0
 dR_3OAS120:            pM_3oddecACP_c + pM_ACP_c + pM_co2_c - pM_dcaACP_c
                        - pM_h_c - pM_malACP_c >= 0
 dR_3OAS100:            pM_3odecACP_c + pM_ACP_c + pM_co2_c - pM_h_c
                        - pM_malACP_c - pM_ocACP_c >= 0
 dR_3OAS60:             pM_3ohexACP_c + pM_ACP_c - pM_butACP_c + pM_co2_c
                        - pM_h_c - pM_malACP_c >= 0
 dR_3OAS140:            pM_3omrsACP_c + pM_ACP_c + pM_co2_c - pM_ddcaACP_c
                        - pM_h_c - pM_malACP_c >= 0
 dR_3OAS80:             pM_3ooctACP_c + pM_ACP_c + pM_co2_c - pM_h_c
                        - pM_hexACP_c - pM_malACP_c >= 0
 dR_3OAS180:            pM_3ooctdACP_c + pM_ACP_c + pM_co2_c - pM_h_c
                        - pM_malACP_c - pM_palmACP_c >= 0
 dR_3OAS160:            pM_3opalmACP_c + pM_ACP_c + pM_co2_c - pM_h_c
                        - pM_malACP_c - pM_myrsACP_c >= 0
 dR_GLYCK:              pM_3pg_c + pM_adp_c - pM_atp_c - pM_glyc_R_c + pM_h_c
                        >= 0
 dR_RBPC:               2 pM_3pg_c - pM_co2_c - pM_h2o_c + 2 pM_h_c
                        - pM_rb15bp_c >= 0
 dR_PGCD:               - pM_3pg_c + pM_3php_c + pM_h_c - pM_nad_c + pM_nadh_c
                        >= 0
 dR_CHORS:              - pM_3psme_c + pM_chor_c + pM_pi_c >= 0
 dR_PSCVTi:             pM_3psme_c - pM_pep_c + pM_pi_c - pM_skm5p_c >= 0
 dR_NABTNO:             pM_4aabutn_c - pM_h2o_c + 2 pM_h_c - pM_n4abutn_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_GLUDC:              pM_4abut_c + pM_co2_c - pM_glu_L_c - pM_h_c >= 0
 dR_ABTA:               - pM_4abut_c - pM_akg_c + pM_glu_L_c + pM_sucsal_c >= 0
 dR_ABUTD:              pM_4abut_c - pM_4abutn_c - pM_h2o_c + 2 pM_h_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_GUDBUTNAH:          pM_4abut_c - pM_4gudbutn_c - pM_h2o_c + pM_urea_c >= 0
 dR_ADCL:               pM_4abz_c - pM_4adcho_c + pM_h_c + pM_pyr_c >= 0
 dR_ADCS:               pM_4adcho_c - pM_chor_c - pM_gln_L_c + pM_glu_L_c >= 0
 dR_AHMMPS:             pM_4ahmmp_c - pM_air_c + pM_gcald_c - 2 pM_h_c
                        + pM_pi_c >= 0
 dR_AMID:               - pM_4gudbd_c + pM_4gudbutn_c - pM_h2o_c + pM_nh4_c
                        >= 0
 dR_GDBTALDH:           pM_4gudbutn_c - pM_gdbtal_c - pM_h2o_c + pM_h_c
                        - pM_nad_c + pM_nadh_c  = 0
 dR_DDPGA:              - pM_4h2oglt_c + pM_glx_c + pM_pyr_c  = 0
 dR_EHGLAT:             pM_4h2oglt_c - pM_akg_c - pM_e4hglu_c + pM_glu_L_c >= 0
 dR_THZSN:              pM_4hba_c + pM_4mhetz_c + pM_co2_c - pM_cys_L_c
                        - pM_dxyl_c - pM_fdxo_2_2_c + pM_fdxr_2_2_c + pM_h2o_c
                        + pM_h_c + pM_nh4_c + pM_pyr_c - pM_tyr_L_c >= 0
 dR_THZPSN:             pM_4hba_c + pM_4mpetz_c + pM_ala_L_c + pM_amp_c
                        - pM_atp_c + pM_co2_c - pM_cys_L_c - pM_dxyl5p_c
                        + pM_h2o_c + pM_h_c + pM_ppi_c - pM_tyr_L_c >= 0
 dR_PROAKGOX1r:         pM_4hpro_LT_c - pM_akg_c + pM_co2_c - pM_o2_c
                        - pM_pro_L_c + pM_succ_c  = 0
 dR_LEUTAi:             - pM_4mop_c + pM_akg_c - pM_glu_L_c + pM_leu_L_c >= 0
 dR_LEUTA:              pM_4mop_c - pM_akg_c + pM_glu_L_c - pM_leu_L_c  = 0
 dR_ASPKi:              pM_4pasp_c + pM_adp_c - pM_asp_L_c - pM_atp_c >= 0
 dR_ASADi:              - pM_4pasp_c + pM_aspsa_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c + pM_pi_c >= 0
 dR_PPNCL3:             - pM_4ppan_c + pM_4ppcys_c + pM_amp_c - pM_atp_c
                        - pM_cys_L_c + pM_h_c + pM_ppi_c >= 0
 dR_PPNCL:              - pM_4ppan_c + pM_4ppcys_c + pM_cdp_c - pM_ctp_c
                        - pM_cys_L_c + pM_h_c + pM_pi_c >= 0
 dR_PNTK:               pM_4ppan_c + pM_adp_c - pM_atp_c + pM_h_c - pM_pnto_R_c
                        >= 0
 dR_PPCDC:              - pM_4ppcys_c + pM_co2_c - pM_h_c + pM_pan4p_c >= 0
 dR_RBFSb:              pM_4r5au_c - 2 pM_dmlz_c + pM_ribflv_c >= 0
 dR_RBFSa:              - pM_4r5au_c - pM_db4p_c + pM_dmlz_c + 2 pM_h2o_c
                        + pM_pi_c >= 0
 dR_PMDPHT:             pM_4r5au_c - pM_5aprbu_c - pM_h2o_c + pM_pi_c >= 0
 dR_AIRCr:              pM_5aizc_c - pM_air_c - pM_co2_c + pM_h_c  = 0
 dR_PPBNGS:             - 2 pM_5aop_c + 2 pM_h2o_c + pM_h_c + pM_ppbng_c >= 0
 dR_G1SATi:             pM_5aop_c - pM_glu1sa_c >= 0
 dR_APRAUR:             pM_5aprbu_c - pM_5apru_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_FOMETR:             - pM_5fthf_c + pM_h2o_c - pM_h_c + pM_methf_c  = 0
 dR_MNXR7081:           - pM_5g2oxpt_c + pM_co2_c + pM_gdbtal_c >= 0
 dR_ARGDH:              pM_5g2oxpt_c - pM_arg_L_c - pM_fad_c + pM_fadh2_c
                        - pM_h2o_c + pM_nh4_c >= 0
 dR_MTRI:               - pM_5mdr1p_c + pM_5mdru1p_c >= 0
 dR_MTAP:               pM_5mdr1p_c - pM_5mta_c + pM_ade_c - pM_pi_c >= 0
 dR_DM_5mdru1p:         - pM_5mdru1p_c >= 0
 dR_SPMS:               pM_5mta_c - pM_ametam_c + pM_h_c - pM_ptrc_c
                        + pM_spmd_c >= 0
 dR_MTHFR2:             pM_5mthf_c - 2 pM_h_c - pM_mlthf_c + pM_nad_c
                        - pM_nadh_c >= 0
 dR_METS:               - pM_5mthf_c + pM_h_c - pM_hcys_L_c + pM_met_L_c
                        + pM_thf_c >= 0
 dR_GND:                - pM_6pgc_c + pM_co2_c - pM_nadp_c + pM_nadph_c
                        + pM_ru5p_D_c >= 0
 dR_PGL:                pM_6pgc_c - pM_6pgl_c - pM_h2o_c + pM_h_c >= 0
 dR_G6PBDH:             pM_6pgl_c - pM_g6p_B_c + pM_h_c - pM_nadp_c
                        + pM_nadph_c >= 0
 dR_G6PDH2:             pM_6pgl_c - pM_g6p_c + pM_h_c - pM_nadp_c + pM_nadph_c
                        >= 0
 dR_PTHPS:              pM_6pthp_c - pM_ahdt_c + pM_pppi_c >= 0
 dR_AOXSr:              pM_8aonn_c - pM_ala_L_c + pM_co2_c + pM_coa_c - pM_h_c
                        - pM_pmcoa_c  = 0
 dR_ACOATA_syn:         - pM_ACP_c + pM_acACP_c - pM_accoa_c + pM_coa_c  = 0
 dR_KAS14:              pM_ACP_c - pM_acACP_c + pM_actACP_c + pM_co2_c - pM_h_c
                        - pM_malACP_c >= 0
 dR_MCOATA:             - pM_ACP_c + pM_coa_c + pM_malACP_c - pM_malcoa_c  = 0
 dR_ADPT2:              pM_C04051_c - pM_aicar_c - pM_ppi_c + pM_prpp_c  = 0
 dR_DM_PHB:             - pM_PHB_c >= 0
 dR_BKTIO_syn:          pM_aacoa_c - 2 pM_accoa_c + pM_coa_c >= 0
 dR_R05219:             pM_ac_c + pM_ahcys_c - pM_amet_c - pM_h2o_c + pM_h_c
                        - pM_pre5_c + pM_pre6a_c >= 0
 dR_UHGADA:             pM_ac_c - pM_h2o_c - pM_u3aga_c + pM_u3hga_c >= 0
 dR_ALDD2xr:            pM_ac_c - pM_acald_c - pM_h2o_c + 2 pM_h_c - pM_nad_c
                        + pM_nadh_c  = 0
 dR_ACt2rpp:            pM_ac_c - pM_ac_p + pM_h_c - pM_h_p  = 0
 dR_CYSS:               pM_ac_c - pM_acser_c + pM_cys_L_c - pM_h2s_c + pM_h_c
                        >= 0
 dR_METACH:             pM_ac_c - pM_achms_c - pM_h2s_c + pM_h_c + pM_hcys_L_c
                        >= 0
 dR_ACKr:               - pM_ac_c + pM_actp_c + pM_adp_c - pM_atp_c  = 0
 dR_ACS:                - pM_ac_c + pM_accoa_c + pM_amp_c - pM_atp_c - pM_coa_c
                        + pM_ppi_c >= 0
 dR_EX_ac_e:            - pM_ac_e >= 0
 dR_ACtex:              - pM_ac_e + pM_ac_p  = 0
 dR_ALCD2y:             pM_acald_c - pM_etoh_c + pM_h_c - pM_nadp_c
                        + pM_nadph_c  = 0
 dR_PDH:                pM_accoa_c + pM_co2_c - pM_coa_c - pM_nad_c + pM_nadh_c
                        - pM_pyr_c >= 0
 dR_ACCOACr:            - pM_accoa_c + pM_adp_c - pM_atp_c + pM_h_c - pM_hco3_c
                        + pM_malcoa_c + pM_pi_c  = 0
 dR_PTAr:               - pM_accoa_c + pM_actp_c + pM_coa_c - pM_pi_c  = 0
 dR_PDHbr:              pM_accoa_c - pM_adhlam_c - pM_coa_c + pM_dhlam_c  = 0
 dR_ACGS:               - pM_accoa_c + pM_acglu_c + pM_coa_c - pM_glu_L_c
                        + pM_h_c >= 0
 dR_HSERTA:             - pM_accoa_c + pM_achms_c + pM_coa_c - pM_hom_L_c >= 0
 dR_KAS15:              - pM_accoa_c + pM_actACP_c + pM_co2_c + pM_coa_c
                        - pM_h_c - pM_malACP_c >= 0
 dR_CS:                 - pM_accoa_c + pM_cit_c + pM_coa_c - pM_h2o_c + pM_h_c
                        - pM_oaa_c >= 0
 dR_POR_syn:            pM_accoa_c + pM_co2_c - pM_coa_c - 2 pM_fdxo_2_2_c
                        + 2 pM_fdxr_2_2_c + pM_h_c - pM_pyr_c >= 0
 dR_SERATi:             - pM_accoa_c + pM_acser_c + pM_coa_c - pM_ser_L_c >= 0
 dR_G1PACT:             - pM_accoa_c + pM_acgam1p_c + pM_coa_c - pM_gam1p_c
                        + pM_h_c >= 0
 dR_AGPR:               pM_acg5p_c - pM_acg5sa_c + pM_h_c - pM_nadp_c
                        + pM_nadph_c - pM_pi_c  = 0
 dR_ACGK:               pM_acg5p_c - pM_acglu_c + pM_adp_c - pM_atp_c >= 0
 dR_ACOTA:              pM_acg5sa_c - pM_acorn_c - pM_akg_c + pM_glu_L_c  = 0
 dR_UAGDP:              - pM_acgam1p_c - pM_h_c + pM_ppi_c + pM_uacgam_c
                        - pM_utp_c >= 0
 dR_ORNTACr:            pM_acglu_c - pM_acorn_c - pM_glu_L_c + pM_orn_c  = 0
 dR_UAG2EMA:            pM_acmana_c - pM_h2o_c + pM_h_c - pM_uacgam_c
                        + pM_udp_c  = 0
 dR_ACONT2:             - pM_acon_C_c - pM_h2o_c + pM_icit_c  = 0
 dR_ACONT1:             pM_acon_C_c - pM_cit_c + pM_h2o_c  = 0
 dR_ACLDC:              pM_actn_R_c - pM_alac_S_c + pM_co2_c - pM_h_c  = 0
 dR_CYRDAAT:            pM_adcobdam_c - pM_atp_c - pM_co1dam_c - pM_h_c
                        + pM_pppi_c >= 0
 dR_ADCYRS:             - pM_adcobdam_c + pM_adcobhex_c + 4 pM_adp_c
                        - 4 pM_atp_c - 4 pM_gln_L_c + 4 pM_glu_L_c - 4 pM_h2o_c
                        + 4 pM_h_c + 4 pM_pi_c >= 0
 dR_ADCPS2:             - pM_adcobhex_c + pM_adocbip_c + pM_adp_c - pM_applp_c
                        - pM_atp_c + pM_h_c + pM_pi_c >= 0
 dR_ADCPS1:             - pM_adcobhex_c + pM_adocbi_c + pM_adp_c - pM_appl_c
                        - pM_atp_c + pM_h_c + pM_pi_c >= 0
 dR_PUNP1:              pM_ade_c - pM_adn_c - pM_pi_c + pM_r1p_c  = 0
 dR_ADPT:               - pM_ade_c + pM_amp_c + pM_ppi_c - pM_prpp_c >= 0
 dR_AMPN:               pM_ade_c - pM_amp_c - pM_h2o_c + pM_r5p_c >= 0
 dR_PDHa:               pM_adhlam_c + pM_co2_c - pM_h_c - pM_lpam_c - pM_pyr_c
                        >= 0
 dR_NTD7:               pM_adn_c - pM_amp_c - pM_h2o_c + pM_pi_c >= 0
 dR_AHCi:               pM_adn_c - pM_ahcys_c - pM_h2o_c + pM_hcys_L_c >= 0
 dR_ACBIPGT:            - pM_adocbip_c + pM_agdpcbi_c - pM_gtp_c - pM_h_c
                        + pM_ppi_c >= 0
 dR_ADOCBLS:            pM_adocbl_c - pM_agdpcbi_c + pM_gmp_c + pM_h_c
                        - pM_rdmbzi_c >= 0
 dR_CTPS2:              pM_adp_c - pM_atp_c + pM_ctp_c - pM_gln_L_c
                        + pM_glu_L_c - pM_h2o_c + 2 pM_h_c + pM_pi_c - pM_utp_c
                        >= 0
 dR_CTPS1:              pM_adp_c - pM_atp_c + pM_ctp_c + 2 pM_h_c - pM_nh4_c
                        + pM_pi_c - pM_utp_c >= 0
 dR_TMK:                pM_adp_c - pM_atp_c + pM_h_c - pM_thm_c + pM_thmmp_c
                        >= 0
 dR_FE2abcpp:           pM_adp_c - pM_atp_c + pM_fe2_c - pM_fe2_p - pM_h2o_c
                        + pM_h_c + pM_pi_c >= 0
 dR_COBALT2abcpp:       pM_adp_c - pM_atp_c - pM_cobalt2_c + pM_cobalt2_p
                        - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_PPK2r:              pM_adp_c - pM_atp_c - pM_ppi_c + pM_pppi_c  = 0
 dR_ATPS4rpp:           - 3 pM_adp_c + 3 pM_atp_c + 3 pM_h2o_c + 11 pM_h_c
                        - 14 pM_h_p - 3 pM_pi_c >= 0
 dR_SULabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_so4_c - pM_so4_p >= 0
 dR_GGPS:               pM_adp_c - pM_adpglc_c + pM_glcglycp_c - pM_glyc3p_c
                        + pM_h_c >= 0
 dR_SUCRabcpp_syn:      pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_sucr_c - pM_sucr_p >= 0
 dR_ZN2abcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        - pM_zn2_c + pM_zn2_p >= 0
 dR_PRFGS:              pM_adp_c - pM_atp_c - pM_fgam_c + pM_fpram_c
                        - pM_gln_L_c + pM_glu_L_c - pM_h2o_c + pM_h_c + pM_pi_c
                        >= 0
 dR_GLUK_syn:           pM_adp_c - pM_atp_c + pM_g6p_B_c - pM_glc_B_c + pM_h_c
                        >= 0
 dR_NDPK10:             pM_adp_c - pM_atp_c - pM_didp_c + pM_ditp_c  = 0
 dR_BCT1_syn:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_hco3_c
                        - pM_hco3_p + pM_pi_c >= 0
 dR_CBMK:               pM_adp_c - pM_atp_c + pM_cbp_c - pM_co2_c + 2 pM_h_c
                        - pM_nh4_c >= 0
 dR_HISabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_his_L_c
                        - pM_his_L_p + pM_pi_c >= 0
 dR_PIuabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + 2 pM_pi_c
                        - pM_pi_p >= 0
 dR_DPCOAK:             pM_adp_c - pM_atp_c + pM_coa_c - pM_dpcoa_c + pM_h_c
                        >= 0
 dR_UMPK:               pM_adp_c - pM_atp_c + pM_udp_c - pM_ump_c  = 0
 dR_THFGLUS:            pM_adp_c - pM_atp_c - pM_glu_L_c + pM_h_c + pM_pi_c
                        - pM_thf_c + pM_thfglu_c  = 0
 dR_HSK:                pM_adp_c - pM_atp_c + pM_h_c - pM_hom_L_c + pM_phom_c
                        >= 0
 dR_NDPK9:              pM_adp_c - pM_atp_c - pM_idp_c + pM_itp_c  = 0
 dR_NDPK8:              pM_adp_c - pM_atp_c - pM_dadp_c + pM_datp_c  = 0
 dR_NDPK7:              pM_adp_c - pM_atp_c - pM_dcdp_c + pM_dctp_c  = 0
 dR_NDPK6:              pM_adp_c - pM_atp_c - pM_dudp_c + pM_dutp_c  = 0
 dR_NDPK5:              pM_adp_c - pM_atp_c - pM_dgdp_c + pM_dgtp_c  = 0
 dR_NDPK4:              pM_adp_c - pM_atp_c - pM_dtdp_c + pM_dttp_c  = 0
 dR_NDPK3:              pM_adp_c - pM_atp_c - pM_cdp_c + pM_ctp_c  = 0
 dR_NDPK2:              pM_adp_c - pM_atp_c - pM_udp_c + pM_utp_c  = 0
 dR_NDPK1:              pM_adp_c - pM_atp_c - pM_gdp_c + pM_gtp_c  = 0
 dR_PFK:                pM_adp_c - pM_atp_c - pM_f6p_c + pM_fdp_c + pM_h_c >= 0
 dR_COBALT2abcpp_syn:   pM_adp_c - pM_atp_c + pM_cobalt2_c - pM_cobalt2_p
                        - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_PROabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_pro_L_c - pM_pro_L_p >= 0
 dR_TMPKr:              pM_adp_c - pM_atp_c - pM_thmmp_c + pM_thmpp_c  = 0
 dR_GLUCYS:             pM_adp_c - pM_atp_c - pM_cys_L_c - pM_glu_L_c
                        + pM_glucys_c + pM_h_c + pM_pi_c >= 0
 dR_CU2abcpp_syn:       pM_adp_c - pM_atp_c + pM_cu2_c - pM_cu2_p - pM_h2o_c
                        + pM_h_c + pM_pi_c >= 0
 dR_SHKK:               pM_adp_c - pM_atp_c + pM_h_c + pM_skm5p_c - pM_skm_c
                        >= 0
 dR_BiomassProtein:     4.32 pM_adp_c - 0.0849441933 pM_ala_L_c
                        - 0.0504230471 pM_arg_L_c - 0.040560173 pM_asn_L_c
                        - 0.0503230488 pM_asp_L_c - 4.32 pM_atp_c
                        + pM_biomass_protein_c - 0.0099637419 pM_cys_L_c
                        - 0.0556555659 pM_gln_L_c - 0.0605880899 pM_glu_L_c
                        - 0.0738726444 pM_gly_c - 4.32 pM_h2o_c + 4.32 pM_h_c
                        - 0.0186190281 pM_his_L_c - 0.0628010954 pM_ile_L_c
                        - 0.1141306486 pM_leu_L_c - 0.0418688462 pM_lys_L_c
                        - 0.0196255325 pM_met_L_c - 0.040158006 pM_phe_L_c
                        + 4.32 pM_pi_c - 0.0513273794 pM_pro_L_c
                        - 0.0580707418 pM_ser_L_c - 0.0549512302 pM_thr_L_c
                        - 0.0154990816 pM_trp_L_c - 0.0291864553 pM_tyr_L_c
                        - 0.0674314506 pM_val_L_c >= 0
 dR_R05224:             2 pM_adp_c - 2 pM_atp_c + pM_fdxo_2_2_c - pM_fdxr_2_2_c
                        - 2 pM_gln_L_c + 2 pM_glu_L_c - pM_h2o_c + pM_h_c
                        + pM_hgbam_c - pM_hgbyr_c + pM_ppi_c >= 0
 dR_HEX1:               pM_adp_c - pM_atp_c + pM_g6p_c - pM_glc_D_c + pM_h_c
                        >= 0
 dR_HEX7:               pM_adp_c - pM_atp_c + pM_f6p_c - pM_fru_c + pM_h_c >= 0
 dR_UAMAS:              pM_adp_c - pM_ala_L_c - pM_atp_c + pM_h_c + pM_pi_c
                        + pM_uama_c - pM_uamr_c >= 0
 dR_ALAabcpp:           pM_adp_c + pM_ala_L_c - pM_ala_L_p - pM_atp_c
                        - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_TMDK1:              pM_adp_c - pM_atp_c + pM_dtmp_c + pM_h_c - pM_thymd_c
                        >= 0
 dR_RNDR1:              - pM_adp_c + pM_dadp_c + pM_h2o_c + pM_trdox_c
                        - pM_trdrd_c >= 0
 dR_CU2abcu_syn:        pM_adp_c - pM_atp_c - pM_cu2_c + pM_cu2_u - pM_h2o_c
                        + pM_h_c + pM_pi_c >= 0
 dR_PYK:                - pM_adp_c + pM_atp_c - pM_h_c - pM_pep_c + pM_pyr_c
                        >= 0
 dR_GLYabcpp:           pM_adp_c - pM_atp_c + pM_gly_c - pM_gly_p - pM_h2o_c
                        + pM_h_c + pM_pi_c >= 0
 dR_GTHS:               pM_adp_c - pM_atp_c - pM_glucys_c - pM_gly_c
                        + pM_gthrd_c + pM_h_c + pM_pi_c >= 0
 dR_PRUK:               pM_adp_c - pM_atp_c + pM_h_c + pM_rb15bp_c
                        - pM_ru5p_D_c >= 0
 dR_ADSK:               pM_adp_c - pM_aps_c - pM_atp_c + pM_h_c + pM_paps_c
                        >= 0
 dR_DTMPK:              pM_adp_c - pM_atp_c + pM_dtdp_c - pM_dtmp_c  = 0
 dR_NTP1:               pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_NI2abcpp_syn:       pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_ni2_c
                        - pM_ni2_p + pM_pi_c >= 0
 dR_SERabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_ser_L_c - pM_ser_L_p >= 0
 dR_GLCS1:              pM_adp_c - pM_adpglc_c + pM_glycogen_c + pM_h_c >= 0
 dR_CYTK1:              pM_adp_c - pM_atp_c + pM_cdp_c - pM_cmp_c  = 0
 dR_ALAALAr:            pM_adp_c - 2 pM_ala_D_c + pM_alaala_c - pM_atp_c
                        + pM_h_c + pM_pi_c  = 0
 dR_GLCGLYCabcpp_syn:   pM_adp_c - pM_atp_c + pM_glcglyc_c - pM_glcglyc_p
                        - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_GALK:               pM_adp_c - pM_atp_c + pM_gal1p_c - pM_gal_c + pM_h_c
                        >= 0
 dR_GK1:                pM_adp_c - pM_atp_c + pM_gdp_c - pM_gmp_c  = 0
 dR_PTRCabcpp:          pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_ptrc_c - pM_ptrc_p >= 0
 dR_BiomassMixo:        16.88006918 pM_adp_c - 16.88006918 pM_atp_c
                        - 0.100827 pM_biomass_dna_c
                        - 0.2327419 pM_biomass_ion_c
                        - 0.15719526 pM_biomass_lipid_c
                        - 0.025228 pM_biomass_lps_c
                        - 0.025228 pM_biomass_pept_p
                        - 0.0109472 pM_biomass_pigment_c
                        - 4.600079 pM_biomass_protein_c
                        - 0.530238 pM_biomass_rna_c
                        - 0.06843016 pM_biomass_solublep_c
                        - 12.27996718 pM_h2o_c + 16.88006918 pM_h_c
                        + 16.87248018 pM_pi_c >= 0
 dR_THRPS:              pM_adp_c - pM_atp_c + pM_h_c - pM_thr_L_c + pM_thrp_c
                        >= 0
 dR_UREAabcpp:          pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_urea_c - pM_urea_p >= 0
 dR_UAMAGS:             pM_adp_c - pM_atp_c - pM_glu_D_c + pM_h_c + pM_pi_c
                        - pM_uama_c + pM_uamag_c >= 0
 dR_DHFS:               pM_adp_c - pM_atp_c + pM_dhf_c - pM_dhpt_c - pM_glu_L_c
                        + pM_h_c + pM_pi_c >= 0
 dR_LEUabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_leu_L_c
                        - pM_leu_L_p + pM_pi_c >= 0
 dR_GLNTRAT:            pM_adp_c - pM_atp_c - pM_gln_L_c + pM_glntrna_c
                        + pM_glu_L_c - pM_glutrna_gln_c - pM_h2o_c + pM_h_c
                        + pM_pi_c >= 0
 dR_URIDK2r:            pM_adp_c - pM_atp_c + pM_dudp_c - pM_dump_c  = 0
 dR_ZNabcpp:            pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_zn2_c - pM_zn2_p >= 0
 dR_ARGabcpp:           pM_adp_c + pM_arg_L_c - pM_arg_L_p - pM_atp_c
                        - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_MG2tpp:             pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_mg2_c
                        - pM_mg2_p + pM_pi_c >= 0
 dR_UGLDDS2:            pM_adp_c - pM_alaala_c - pM_atp_c + pM_fdxo_2_2_c
                        - pM_fdxr_2_2_c + pM_pi_c - pM_uGgl_c + pM_uGgla_c >= 0
 dR_COCHL:              pM_adp_c - pM_atp_c + pM_co2dam_c - pM_cobalt2_c
                        - pM_h2o_c + 2 pM_h_c - pM_hgbam_c + pM_pi_c >= 0
 dR_MOBDabcpp:          pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_mobd_c
                        - pM_mobd_p + pM_pi_c >= 0
 dR_BiomassDNA:         3.36 pM_adp_c - 3.36 pM_atp_c + pM_biomass_dna_c
                        - 0.2609519276 pM_datp_c - 0.2382496752 pM_dctp_c
                        - 0.2389439337 pM_dgtp_c - 0.2618544636 pM_dttp_c
                        - 3.36 pM_h2o_c + 3.36 pM_h_c + 3.36 pM_pi_c + pM_ppi_c
                        >= 0
 dR_UGMDDS:             pM_adp_c - pM_alaala_c - pM_atp_c + pM_h_c + pM_pi_c
                        - pM_ugmd_c + pM_ugmda_c >= 0
 dR_NADK:               pM_adp_c - pM_atp_c + pM_h_c - pM_nad_c + pM_nadp_c
                        >= 0
 dR_RBK:                pM_adp_c - pM_atp_c + pM_h_c + pM_r5p_c - pM_rib_D_c
                        >= 0
 dR_BiomassRNA:         2.41 pM_adp_c - 2.6668620129 pM_atp_c
                        + pM_biomass_rna_c - 0.2340458436 pM_ctp_c
                        - 0.2500575213 pM_gtp_c - 2.41 pM_h2o_c + 2.41 pM_h_c
                        + 2.41 pM_pi_c + pM_ppi_c - 0.2590346222 pM_utp_c >= 0
 dR_ATPSu:              - 3 pM_adp_c + 3 pM_atp_c + 3 pM_h2o_c + 11 pM_h_c
                        - 14 pM_h_u - 3 pM_pi_c >= 0
 dR_CBPS:               2 pM_adp_c - 2 pM_atp_c + pM_cbp_c - pM_gln_L_c
                        + pM_glu_L_c - pM_h2o_c + 2 pM_h_c - pM_hco3_c
                        + pM_pi_c >= 0
 dR_ADK1:               2 pM_adp_c - pM_amp_c - pM_atp_c  = 0
 dR_CYPHYS:             2 pM_adp_c - pM_arg_L_c - pM_asp_L_c - 2 pM_atp_c
                        + pM_cyanphy_c + 2 pM_h_c + 2 pM_pi_c - pM_precyanphy_c
                        >= 0
 dR_BiomassHetero:      17.70006718 pM_adp_c - 17.70006918 pM_atp_c
                        - 0.100827 pM_biomass_dna_c
                        - 0.2327419 pM_biomass_ion_c
                        - 0.15719526 pM_biomass_lipid_c
                        - 0.025228 pM_biomass_lps_c
                        - 0.025228 pM_biomass_pept_p
                        - 0.0109472 pM_biomass_pigment_c
                        - 4.600079 pM_biomass_protein_c
                        - 0.530238 pM_biomass_rna_c
                        - 0.06843016 pM_biomass_solublep_c
                        - 13.09996718 pM_h2o_c + 17.70006718 pM_h_c
                        + 17.69248018 pM_pi_c >= 0
 dR_LYSabcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_lys_L_c
                        - pM_lys_L_p + pM_pi_c >= 0
 dR_PFK_2:              pM_adp_c - pM_atp_c + pM_h_c - pM_tag6p_D_c
                        + pM_tagdp_D_c >= 0
 dR_PFK_3:              pM_adp_c - pM_atp_c + pM_h_c + pM_s17bp_c - pM_s7p_c
                        >= 0
 dR_CA2abcpp:           pM_adp_c - pM_atp_c + pM_ca2_c - pM_ca2_p - pM_h2o_c
                        + pM_h_c + pM_pi_c >= 0
 dR_GLYK:               pM_adp_c - pM_atp_c + pM_glyc3p_c - pM_glyc_c + pM_h_c
                        >= 0
 dR_NO3abcpp:           pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_no3_c
                        - pM_no3_p + pM_pi_c >= 0
 dR_PRAGS:              pM_adp_c - pM_atp_c + pM_gar_c - pM_gly_c + pM_h_c
                        + pM_pi_c - pM_pram_c >= 0
 dR_SPMDabcpp:          pM_adp_c - pM_atp_c - pM_h2o_c + pM_h_c + pM_pi_c
                        + pM_spmd_c - pM_spmd_p >= 0
 dR_ADK_syn:            pM_adp_c - pM_atp_c + pM_dadp_c - pM_damp_c  = 0
 dR_GLNabcpp:           pM_adp_c - pM_atp_c + pM_gln_L_c - pM_gln_L_p
                        - pM_h2o_c + pM_h_c + pM_pi_c >= 0
 dR_GLNS:               pM_adp_c - pM_atp_c + pM_gln_L_c - pM_glu_L_c + pM_h_c
                        - pM_nh4_c + pM_pi_c >= 0
 dR_BiomassPept:        0.83 pM_adp_c - 0.83 pM_atp_c + pM_biomass_pept_p
                        - 0.83 pM_h2o_c + 0.83 pM_h_c - pM_peptido_syn_p
                        + 0.83 pM_pi_c >= 0
 dR_RBFK:               pM_adp_c - pM_atp_c + pM_fmn_c + pM_h_c - pM_ribflv_c
                        >= 0
 dR_DBTSr:              pM_adp_c - pM_atp_c - pM_co2_c - pM_dann_c + pM_dtbt_c
                        + 3 pM_h_c + pM_pi_c  = 0
 dR_SUCOAS:             pM_adp_c - pM_atp_c - pM_coa_c + pM_pi_c - pM_succ_c
                        + pM_succoa_c  = 0
 dR_MPML:               pM_adp_c - pM_atp_c - pM_h2o_c + 3 pM_h_c - pM_mg2_c
                        + pM_mppp9_c + pM_pi_c - pM_ppp9_c  = 0
 dR_FE3abc:             pM_adp_c - pM_atp_c + pM_fe3_c - pM_fe3_p - pM_h2o_c
                        + pM_h_c + pM_pi_c >= 0
 dR_GLU5K:              pM_adp_c - pM_atp_c + pM_glu5p_c - pM_glu_L_c >= 0
 dR_PRAIS:              pM_adp_c + pM_air_c - pM_atp_c - pM_fpram_c + 2 pM_h_c
                        + pM_pi_c >= 0
 dR_GLGC:               pM_adpglc_c - pM_atp_c - pM_g1p_c - pM_h_c + pM_ppi_c
                        >= 0
 dR_NADN:               pM_adprib_c - pM_h2o_c + pM_h_c - pM_nad_c + pM_ncam_c
                        >= 0
 dR_ADPRDP:             - pM_adprib_c + pM_amp_c - pM_h2o_c + 2 pM_h_c
                        + pM_r5p_c >= 0
 dR_ARGDC:              pM_agm_c - pM_arg_L_c + pM_co2_c - pM_h_c >= 0
 dR_AGMT:               - pM_agm_c - pM_h2o_c + pM_ptrc_c + pM_urea_c >= 0
 dR_PC6YM:              2 pM_ahcys_c - 2 pM_amet_c + pM_co2_c + 2 pM_h_c
                        - pM_pre6b_c + pM_pre8_c >= 0
 dR_DNAMT_syn:          pM_ahcys_c - pM_amet_c + pM_dna5mtc_c - pM_dna_c
                        + pM_h_c >= 0
 dR_PC11M:              pM_ahcys_c - pM_amet_c + pM_h_c - pM_pre4_c + pM_pre5_c
                        >= 0
 dR_PC20M:              pM_ahcys_c - pM_amet_c - pM_dscl_c + pM_h_c
                        + pM_pre3a_c >= 0
 dR_UPP3MT:             2 pM_ahcys_c - 2 pM_amet_c + pM_dscl_c + pM_h_c
                        - pM_uppg3_c >= 0
 dR_PC17M:              pM_ahcys_c - pM_amet_c + 3 pM_h_c - pM_pre3b_c
                        + pM_pre4_c >= 0
 dR_TOCOPHOM1:          pM_ahcys_c - pM_amet_c + pM_atocophe_c - pM_gtocophe_c
                        + pM_h_c  = 0
 dR_TOCOPHOM2:          pM_ahcys_c - pM_amet_c + pM_btocophe_c - pM_dtocophe_c
                        + pM_h_c  = 0
 dR_DMTPHT:             pM_ahcys_c - pM_amet_c - pM_dmtphllqne_c + pM_h_c
                        + pM_phllqne_c >= 0
 dR_MPOMT:              pM_ahcys_c - pM_amet_c - pM_mppp9_c + pM_mppp9om_c >= 0
 dR_DNTPPA:             - pM_ahdt_c + pM_dhpmp_c - pM_h2o_c + pM_h_c + pM_ppi_c
                        >= 0
 dR_GTPCI:              pM_ahdt_c + pM_for_c - pM_gtp_c - pM_h2o_c + pM_h_c
                        >= 0
 dR_IG3PS:              pM_aicar_c + pM_eig3p_c - pM_gln_L_c + pM_glu_L_c
                        + pM_h_c - pM_prlp_c >= 0
 dR_GLUDy:              pM_akg_c - pM_glu_L_c - pM_h2o_c + pM_h_c - pM_nadp_c
                        + pM_nadph_c + pM_nh4_c  = 0
 dR_ICDHy:              pM_akg_c + pM_co2_c - pM_icit_c - pM_nadp_c
                        + pM_nadph_c >= 0
 dR_AKGtpp:             pM_akg_c - pM_akg_p + pM_h_c - pM_h_p  = 0
 dR_ASPTA:              - pM_akg_c - pM_asp_L_c + pM_glu_L_c + pM_oaa_c  = 0
 dR_PHETA1:             - pM_akg_c + pM_glu_L_c - pM_phe_L_c + pM_phpyr_c  = 0
 dR_ORNTA:              - pM_akg_c + pM_glu5sa_c + pM_glu_L_c - pM_orn_c >= 0
 dR_SDPTA:              - pM_akg_c + pM_glu_L_c - pM_sl26da_c + pM_sl2a6o_c
                         = 0
 dR_OXGDC2:             - pM_akg_c + pM_co2_c - pM_h_c + pM_ssaltpp_c
                        - pM_thmpp_c >= 0
 dR_GLMS_syn:           - pM_akg_c + 2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c
                        - pM_gln_L_c + 2 pM_glu_L_c - 2 pM_h_c >= 0
 dR_OXGDC:              - pM_akg_c + pM_co2_c - pM_h_c + pM_sucsal_c >= 0
 dR_GLUSx:              - pM_akg_c - pM_gln_L_c + 2 pM_glu_L_c - pM_h_c
                        + pM_nad_c - pM_nadh_c >= 0
 dR_HSTPT:              pM_akg_c - pM_glu_L_c + pM_hisp_c - pM_imacp_c >= 0
 dR_CYSTA:              - pM_akg_c - pM_cys_L_c + pM_glu_L_c + pM_mercppyr_c
                         = 0
 dR_AKGtex:             - pM_akg_e + pM_akg_p  = 0
 dR_EX_akg_e:           - pM_akg_e >= 0
 dR_PANTS:              - pM_ala_B_c + pM_amp_c - pM_atp_c + pM_h_c
                        - pM_pant_R_c + pM_pnto_R_c + pM_ppi_c >= 0
 dR_ASP1DC:             pM_ala_B_c - pM_asp_L_c + pM_co2_c - pM_h_c >= 0
 dR_BAMPPALDOX:         pM_ala_B_c - pM_bamppald_c - pM_h2o_c + 2 pM_h_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_ALAR:               pM_ala_D_c - pM_ala_L_c  = 0
 dR_SPT_syn2:           - pM_ala_L_c - pM_glx_c + pM_gly_c + pM_pyr_c  = 0
 dR_ALATRS:             - pM_ala_L_c + pM_alatrna_c + pM_amp_c - pM_atp_c
                        + pM_ppi_c - pM_trnaala_c >= 0
 dR_ALAD_syn:           pM_ala_L_c + pM_h2o_c - pM_h_c + pM_nad_c - pM_nadh_c
                        - pM_nh4_c - pM_pyr_c >= 0
 dR_EX_ala_L_e:         - pM_ala_L_e >= 0
 dR_ALAtex:             - pM_ala_L_e + pM_ala_L_p  = 0
 dR_ACLS:               pM_alac_S_c + pM_co2_c - pM_h_c - 2 pM_pyr_c >= 0
 dR_BTS4:               - pM_amet_c + pM_btn_c + pM_dad_5_c - pM_dtbt_c
                        + pM_h_c + pM_met_L_c - pM_s_c >= 0
 dR_METAT:              pM_amet_c - pM_atp_c - pM_h2o_c - pM_met_L_c + pM_pi_c
                        + pM_ppi_c >= 0
 dR_ADMDC:              - pM_amet_c + pM_ametam_c + pM_co2_c - pM_h_c >= 0
 dR_CPPPGO2:            - 2 pM_amet_c + 2 pM_co2_c - pM_cpppg3_c + 2 pM_dad_5_c
                        + 2 pM_met_L_c + pM_pppg9_c >= 0
 dR_SERTRS:             pM_amp_c - pM_atp_c + pM_ppi_c - pM_ser_L_c
                        + pM_sertrna_c - pM_trnaser_c >= 0
 dR_GLUTRS:             pM_amp_c - pM_atp_c - pM_glu_L_c + pM_glutrna_c
                        + pM_ppi_c - pM_trnaglu_c >= 0
 dR_THRTRS:             pM_amp_c - pM_atp_c + pM_ppi_c - pM_thr_L_c
                        + pM_thrtrna_c - pM_trnathr_c >= 0
 dR_PROTRS:             pM_amp_c - pM_atp_c + pM_ppi_c - pM_pro_L_c
                        + pM_protrna_c - pM_trnapro_c >= 0
 dR_TRPTRS:             pM_amp_c - pM_atp_c + pM_ppi_c - pM_trnatrp_c
                        - pM_trp_L_c + pM_trptrna_c >= 0
 dR_ADSL1r:             pM_amp_c - pM_dcamp_c + pM_fum_c  = 0
 dR_GMPS:               pM_amp_c - pM_atp_c + pM_gmp_c + 2 pM_h_c - pM_nh4_c
                        + pM_ppi_c - pM_xmp_c >= 0
 dR_ASPTRS:             pM_amp_c - pM_asp_L_c + pM_asptrna_c - pM_atp_c
                        + pM_ppi_c - pM_trnaasp_c >= 0
 dR_BPNT:               pM_amp_c - pM_h2o_c - pM_pap_c + pM_pi_c >= 0
 dR_PRPPS:              pM_amp_c - pM_atp_c + pM_h_c + pM_prpp_c - pM_r5p_c
                         = 0
 dR_SUCBZL:             pM_amp_c - pM_atp_c - pM_coa_c + pM_ppi_c + pM_sbzcoa_c
                        - pM_sucbz_c >= 0
 dR_NADS2:              pM_amp_c - pM_atp_c - pM_dnad_c - pM_gln_L_c
                        + pM_glu_L_c - pM_h2o_c + pM_h_c + pM_nad_c + pM_ppi_c
                        >= 0
 dR_PPS:                pM_amp_c - pM_atp_c - pM_h2o_c + 2 pM_h_c + pM_pep_c
                        + pM_pi_c - pM_pyr_c >= 0
 dR_ASNTRS:             pM_amp_c - pM_asn_L_c + pM_asntrna_c - pM_atp_c
                        + pM_ppi_c - pM_trnaasn_c >= 0
 dR_ARGTRS:             pM_amp_c - pM_arg_L_c + pM_argtrna_c - pM_atp_c
                        + pM_ppi_c - pM_trnaarg_c >= 0
 dR_LYSTRS:             pM_amp_c - pM_atp_c - pM_lys_L_c + pM_lystrna_c
                        + pM_ppi_c - pM_trnalys_c >= 0
 dR_ILETRS:             pM_amp_c - pM_atp_c - pM_ile_L_c + pM_iletrna_c
                        + pM_ppi_c - pM_trnaile_c >= 0
 dR_LEUTRS:             pM_amp_c - pM_atp_c - pM_leu_L_c + pM_leutrna_c
                        + pM_ppi_c - pM_trnaleu_c >= 0
 dR_GLYTRS:             pM_amp_c - pM_atp_c - pM_gly_c + pM_glytrna_c
                        + pM_ppi_c - pM_trnagly_c >= 0
 dR_ARGSS:              pM_amp_c + pM_argsuc_c - pM_asp_L_c - pM_atp_c
                        - pM_citr_L_c + pM_h_c + pM_ppi_c >= 0
 dR_HISTRS:             pM_amp_c - pM_atp_c - pM_his_L_c + pM_histrna_c
                        + pM_ppi_c - pM_trnahis_c >= 0
 dR_ASNS1:              pM_amp_c + pM_asn_L_c - pM_asp_L_c - pM_atp_c
                        - pM_gln_L_c + pM_glu_L_c - pM_h2o_c + pM_h_c
                        + pM_ppi_c >= 0
 dR_CYSTRS:             pM_amp_c - pM_atp_c - pM_cys_L_c + pM_cystrna_c
                        + pM_ppi_c - pM_trnacys_c >= 0
 dR_METTRS:             pM_amp_c - pM_atp_c - pM_met_L_c + pM_mettrna_c
                        + pM_ppi_c - pM_trnamet_c >= 0
 dR_TYRTRS:             pM_amp_c - pM_atp_c + pM_ppi_c - pM_trnatyr_c
                        - pM_tyr_L_c + pM_tyrtrna_c >= 0
 dR_GMPS2:              pM_amp_c - pM_atp_c - pM_gln_L_c + pM_glu_L_c
                        + pM_gmp_c - pM_h2o_c + 2 pM_h_c + pM_ppi_c - pM_xmp_c
                        >= 0
 dR_VALTRS:             pM_amp_c - pM_atp_c + pM_ppi_c - pM_trnaval_c
                        - pM_val_L_c + pM_valtrna_c >= 0
 dR_PHETRS:             pM_amp_c - pM_atp_c - pM_phe_L_c + pM_phetrna_c
                        + pM_ppi_c - pM_trnaphe_c  = 0
 dR_ANS:                pM_anth_c - pM_chor_c - pM_gln_L_c + pM_glu_L_c
                        + pM_h_c + pM_pyr_c >= 0
 dR_ANPRT:              - pM_anth_c + pM_ppi_c + pM_pran_c - pM_prpp_c >= 0
 dR_ANS2:               pM_anth_c - pM_chor_c + pM_h2o_c + pM_h_c - pM_nh4_c
                        + pM_pyr_c >= 0
 dR_THRPD:              pM_applp_c + pM_co2_c - pM_h_c - pM_thrp_c >= 0
 dR_SADT:               pM_aps_c - pM_atp_c - pM_h_c + pM_ppi_c - pM_so4_c >= 0
 dR_CYPHYC:             pM_arg_L_c + pM_asp_L_c - pM_cyanphy_c - 2 pM_h2o_c
                        + pM_precyanphy_c >= 0
 dR_ARGN:               - pM_arg_L_c - pM_h2o_c + pM_orn_c + pM_urea_c >= 0
 dR_ARGD:               - pM_arg_L_c + pM_citr_L_c - pM_h2o_c + pM_nh4_c >= 0
 dR_ARGSL:              pM_arg_L_c - pM_argsuc_c + pM_fum_c  = 0
 dR_ARGtex:             - pM_arg_L_e + pM_arg_L_p  = 0
 dR_EX_arg_L_e:         - pM_arg_L_e >= 0
 dR_ASNN:               - pM_asn_L_c + pM_asp_L_c - pM_h2o_c + pM_nh4_c >= 0
 dR_ASPCT:              - pM_asp_L_c + pM_cbasp_c - pM_cbp_c + pM_h_c + pM_pi_c
                        >= 0
 dR_ASPO1:              - pM_asp_L_c + pM_h2o2_c - pM_h2o_c + pM_nh4_c
                        - pM_o2_c + pM_oaa_c >= 0
 dR_ASPO6:              - pM_asp_L_c + pM_h2o2_c + pM_h_c + pM_iasp_c - pM_o2_c
                        >= 0
 dR_ADSS:               - pM_asp_L_c + pM_dcamp_c + pM_gdp_c - pM_gtp_c
                        + 2 pM_h_c - pM_imp_c + pM_pi_c >= 0
 dR_HSDx:               pM_aspsa_c + pM_h_c - pM_hom_L_c - pM_nad_c + pM_nadh_c
                         = 0
 dR_BiomassPigment:     - 0.0240463315 pM_atocophe_c
                        - 0.3126096171 pM_bcarote_c + pM_biomass_pigment_c
                        - 0.00132253 pM_btocophe_c - 0.0204399298 pM_dtocophe_c
                        - 0.2300314236 pM_echin_c - 0.14745323 pM_gcarote_c
                        - 0.00132253 pM_gtocophe_c - 0.265419468 pM_zeax_c >= 0
 dR_FMNAT:              - pM_atp_c + pM_fad_c - pM_fmn_c - pM_h_c + pM_ppi_c
                        >= 0
 dR_NNATr:              - pM_atp_c + pM_dnad_c - pM_h_c - pM_nicrnt_c
                        + pM_ppi_c  = 0
 dR_NMNAT:              - pM_atp_c - pM_h_c + pM_nad_c - pM_nmn_c + pM_ppi_c
                        >= 0
 dR_ATPPRT:             - pM_atp_c + pM_ppi_c + pM_prbatp_c - pM_prpp_c >= 0
 dR_PTPAT:              - pM_atp_c + pM_dpcoa_c - pM_h_c - pM_pan4p_c
                        + pM_ppi_c  = 0
 dR_ADNCYC:             - pM_atp_c + pM_camp_c + pM_ppi_c >= 0
 dR_GCATENEC:           pM_bcarote_c - pM_gcarote_c >= 0
 dR_BCAROKE:            - pM_bcarote_c + pM_echin_c + pM_h2o_c - pM_o2_c >= 0
 dR_BCAROHX:            - pM_bcarote_c + pM_bcryptox_c + pM_h2o_c - pM_h_c
                        + pM_nadp_c - pM_nadph_c - pM_o2_c >= 0
 dR_BCAROHX2:           - pM_bcryptox_c + pM_h2o_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c - pM_o2_c + pM_zeax_c >= 0
 dR_GLCP2:              - pM_bglycogen_c + pM_g1p_c - pM_pi_c >= 0
 dR_GLBRAN2:            pM_bglycogen_c - pM_glycogen_c >= 0
 dR_BILIRED:            pM_bilirub_c - pM_biliverd_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_HOXG:               pM_biliverd_c + pM_co_c + pM_fe2_c + 3 pM_h2o_c
                        - 5 pM_h_c + 3 pM_nadp_c - 3 pM_nadph_c - 3 pM_o2_c
                        - pM_pheme_c >= 0
 dR_BiomassGlycogen:    pM_biomass_glycogen_c - pM_glycogen_c >= 0
 dR_BiomassIon:         pM_biomass_ion_c - 0.0193862815 pM_ca2_c
                        - 0.0130157054 pM_cobalt2_c - 0.0129241877 pM_cu2_c
                        - 0.0292856594 pM_fe2_c - 0.0292856594 pM_fe3_c
                        - 0.7269855578 pM_k_c - 0.0325394783 pM_mg2_c
                        - 0.0130144164 pM_mn2_c - 0.0130157054 pM_mobd_c
                        - 0.016271243 pM_na1_c - 0.0488094323 pM_nh4_c
                        - 0.016271243 pM_pi_c - 0.016271243 pM_so4_c
                        - 0.0129241877 pM_zn2_c >= 0
 dR_BiomassLipid:       pM_biomass_lipid_c - 0.13211 pM_dgdg160_c
                        - 0.009757 pM_dgdg161_c - 0.00226 pM_dgdg180_c
                        - 0.002953 pM_dgdg181_9_c - 0.002953 pM_dgdg181_c
                        - 0.033727 pM_dgdg182_9_12_c
                        - 0.073361 pM_dgdg183_6_9_12_c
                        - 0.000963 pM_dgdg183_9_12_15_c
                        - 0.00398 pM_dgdg184_6_9_12_15_c
                        - 0.192938 pM_mgdg160_c - 0.01123 pM_mgdg161_c
                        - 0.000497 pM_mgdg180_c - 0.004602 pM_mgdg181_9_c
                        - 0.004602 pM_mgdg181_c - 0.061434 pM_mgdg182_9_12_c
                        - 0.091587 pM_mgdg183_6_9_12_c
                        - 0.00046 pM_mgdg183_9_12_15_c
                        - 0.002761 pM_mgdg184_6_9_12_15_c - 0.05087 pM_pg160_c
                        - 0.000867 pM_pg161_c - 0.001744 pM_pg180_c
                        - 0.003246 pM_pg181_9_c - 0.003246 pM_pg181_c
                        - 0.032967 pM_pg182_9_12_c - 0.000006 pM_pg183_6_9_12_c
                        - 0.008236 pM_pg183_9_12_15_c
                        - 0.000006 pM_pg184_6_9_12_15_c - 0.151958 pM_sqdg160_c
                        - 0.017555 pM_sqdg161_c - 0.003336 pM_sqdg180_c
                        - 0.01952 pM_sqdg181_9_c - 0.01952 pM_sqdg181_c
                        - 0.051887 pM_sqdg182_9_12_c
                        - 0.001606 pM_sqdg183_6_9_12_c
                        - 0.002224 pM_sqdg183_9_12_15_c
                        - 0.000006 pM_sqdg184_6_9_12_15_c >= 0
 dR_BiomassLPS:         pM_biomass_lps_c - pM_lipidAds_c >= 0
 dR_EAR40y:             - pM_but2eACP_c + pM_butACP_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_CA2t3pp:            - pM_ca2_c + pM_ca2_p + pM_h_c - pM_h_p >= 0
 dR_CA2tex:             - pM_ca2_e + pM_ca2_p  = 0
 dR_EX_ca2_e:           - pM_ca2_e  = 0
 dR_DHORTS:             pM_cbasp_c - pM_dhor_S_c - pM_h2o_c + pM_h_c  = 0
 dR_OCBT:               - pM_cbp_c + pM_citr_L_c + pM_h_c - pM_orn_c + pM_pi_c
                         = 0
 dR_EAR121y:            pM_cddec5eACP_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        - pM_t3c5ddeceACP_c >= 0
 dR_T2DECAI:            pM_cdec3eACP_c - pM_tdec2eACP_c  = 0
 dR_PYK4:               - pM_cdp_c + pM_ctp_c - pM_h_c - pM_pep_c + pM_pyr_c
                        >= 0
 dR_RNDR3:              - pM_cdp_c + pM_dcdp_c + pM_h2o_c + pM_trdox_c
                        - pM_trdrd_c >= 0
 dR_PGSA161:            - pM_cdpdhdec9eg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp161_c >= 0
 dR_DASYN161:           pM_cdpdhdec9eg_c - pM_ctp_c - pM_h_c - pM_pa161_c
                        + pM_ppi_c >= 0
 dR_PGSA160:            - pM_cdpdhdecg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp160_c >= 0
 dR_DASYN160:           pM_cdpdhdecg_c - pM_ctp_c - pM_h_c - pM_pa160_c
                        + pM_ppi_c >= 0
 dR_PGSA181:            - pM_cdpdodec11eg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp181_c >= 0
 dR_DASYN181:           pM_cdpdodec11eg_c - pM_ctp_c - pM_h_c - pM_pa181_c
                        + pM_ppi_c >= 0
 dR_DASYN184_6_9_12_15: pM_cdpdodec691215eg_c - pM_ctp_c - pM_h_c
                        - pM_pa184_6_9_12_15_c + pM_ppi_c >= 0
 dR_PGSA184_6_9_12_15:  - pM_cdpdodec691215eg_c + pM_cmp_c - pM_glyc3p_c
                        + pM_h_c + pM_pgp184_6_9_12_15_c >= 0
 dR_PGSA183_6_9_12:     - pM_cdpdodec6912eg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp183_6_9_12_c >= 0
 dR_DASYN183_6_9_12:    pM_cdpdodec6912eg_c - pM_ctp_c - pM_h_c
                        - pM_pa183_6_9_12_c + pM_ppi_c >= 0
 dR_PGSA183_9_12_15:    - pM_cdpdodec91215eg_c + pM_cmp_c - pM_glyc3p_c
                        + pM_h_c + pM_pgp183_9_12_15_c >= 0
 dR_DASYN183_9_12_15:   pM_cdpdodec91215eg_c - pM_ctp_c - pM_h_c
                        - pM_pa183_9_12_15_c + pM_ppi_c >= 0
 dR_PGSA182_9_12:       - pM_cdpdodec912eg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp182_9_12_c >= 0
 dR_DASYN182_9_12:      pM_cdpdodec912eg_c - pM_ctp_c - pM_h_c
                        - pM_pa182_9_12_c + pM_ppi_c >= 0
 dR_PGSA181_9:          - pM_cdpdodec9eg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp181_9_c >= 0
 dR_DASYN181_9:         pM_cdpdodec9eg_c - pM_ctp_c - pM_h_c - pM_pa181_9_c
                        + pM_ppi_c >= 0
 dR_PGSA180:            - pM_cdpdodecg_c + pM_cmp_c - pM_glyc3p_c + pM_h_c
                        + pM_pgp180_c >= 0
 dR_DASYN180:           pM_cdpdodecg_c - pM_ctp_c - pM_h_c - pM_pa180_c
                        + pM_ppi_c >= 0
 dR_MNXR33817:          pM_cdpglc_c - pM_ctp_c - pM_g1p_c - pM_h_c + pM_ppi_c
                        >= 0
 dR_AMPTASECG:          - pM_cgly_c + pM_cys_L_c + pM_gly_c - pM_h2o_c >= 0
 dR_GTHRDH_syn:         pM_cgly_c + pM_glu_L_c - pM_gthrd_c - pM_h2o_c >= 0
 dR_CHPHYS:             - pM_chlld_c + pM_cholphya_c - pM_h_c - pM_phdp_c
                        + pM_ppi_c >= 0
 dR_POR:                pM_chlld_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        - pM_pchlld_c >= 0
 dR_CHORM:              - pM_chor_c + pM_pphn_c >= 0
 dR_ICHORS:             - pM_chor_c + pM_ichor_c  = 0
 dR_CITtpp:             pM_cit_c - pM_cit_p + pM_h_c - pM_h_p  = 0
 dR_EX_cit_e:           - pM_cit_e >= 0
 dR_CITtex:             - pM_cit_e + pM_cit_p  = 0
 dR_NTPP4:              pM_cmp_c - pM_ctp_c - pM_h2o_c + pM_h_c + pM_ppi_c >= 0
 dR_NTD4:               - pM_cmp_c + pM_cytd_c - pM_h2o_c + pM_pi_c >= 0
 dR_CYRDAR:             2 pM_co1dam_c - 2 pM_co2dam_c + pM_h_c + pM_nad_c
                        - pM_nadh_c >= 0
 dR_GLYCL:              pM_co2_c - pM_gly_c + pM_mlthf_c - pM_nad_c + pM_nadh_c
                        + pM_nh4_c - pM_thf_c >= 0
 dR_PPNDH:              pM_co2_c + pM_h2o_c - pM_h_c + pM_phpyr_c - pM_pphn_c
                        >= 0
 dR_CPPPGO:             2 pM_co2_c - pM_cpppg3_c + 2 pM_h2o_c - 2 pM_h_c
                        - pM_o2_c + pM_pppg9_c >= 0
 dR_DHNAPHT:            pM_co2_c - pM_dhna_c + pM_dmtphllqne_c
                        - 2 pM_fdxo_2_2_c + 2 pM_fdxr_2_2_c + pM_h_c
                        - pM_phdp_c + pM_ppi_c >= 0
 dR_HCO3E:              - pM_co2_c - pM_h2o_c + pM_h_c + pM_hco3_c  = 0
 dR_PPC:                - pM_co2_c - pM_h2o_c + pM_h_c + pM_oaa_c - pM_pep_c
                        + pM_pi_c >= 0
 dR_FDH6pp:             pM_co2_c - 2 pM_ficytc6_p + 2 pM_focytc6_p - pM_for_c
                        + pM_h_c >= 0
 dR_CYNTAH:             2 pM_co2_c - pM_cynt_c - 3 pM_h_c - pM_hco3_c
                        + pM_nh4_c >= 0
 dR_DXPS:               pM_co2_c + pM_dxyl5p_c - pM_g3p_c - pM_h_c - pM_pyr_c
                        >= 0
 dR_CO2tpp:             - pM_co2_c + pM_co2_p >= 0
 dR_PDX5PS:             pM_co2_c - pM_dxyl5p_c + 2 pM_h2o_c + pM_h_c - pM_nad_c
                        + pM_nadh_c + pM_pdx5p_c - pM_phthr_c + pM_pi_c >= 0
 dR_ME1:                pM_co2_c - pM_mal_L_c - pM_nadp_c + pM_nadph_c
                        + pM_pyr_c >= 0
 dR_UPPDC2:             4 pM_co2_c + pM_cpppg1_c - 4 pM_h_c - pM_uppg1_c >= 0
 dR_UPPDC1:             4 pM_co2_c + pM_cpppg3_c - 4 pM_h_c - pM_uppg3_c >= 0
 dR_OMPDC:              pM_co2_c - pM_h_c - pM_orot5p_c + pM_ump_c >= 0
 dR_GLYCL_2:            - pM_co2_c + pM_gly_c - pM_mlthf_c + pM_nad_c
                        - pM_nadh_c - pM_nh4_c + pM_thf_c >= 0
 dR_UREA:               pM_co2_c - pM_h2o_c - 2 pM_h_c + 2 pM_nh4_c - pM_urea_c
                        >= 0
 dR_OXADC:              pM_co2_c + pM_for_c - pM_h_c - pM_oxa_c >= 0
 dR_NNDPR:              pM_co2_c - 2 pM_h_c + pM_nicrnt_c + pM_ppi_c
                        - pM_prpp_c - pM_quln_c >= 0
 dR_EX_co2_e:           - pM_co2_e >= 0
 dR_CO2tex:             - pM_co2_e + pM_co2_p  = 0
 dR_NDH1_3u:            - pM_co2_p - pM_h2o_c - 3 pM_h_c + 3 pM_h_u + pM_hco3_c
                        + pM_nadp_c - pM_nadph_c - pM_pq_u + pM_pqh2_u >= 0
 dR_NDH1_4pp:           - pM_co2_p - pM_h2o_c - 3 pM_h_c + 3 pM_h_p + pM_hco3_c
                        + pM_nadp_c - pM_nadph_c - pM_pq_p + pM_pqh2_p >= 0
 dR_Cotpp:              pM_co_c - pM_co_p  = 0
 dR_COt:                - pM_co_e + pM_co_p  = 0
 dR_EX_co_e:            - pM_co_e >= 0
 dR_NPHS:               pM_coa_c + pM_dhna_c - pM_sbzcoa_c >= 0
 dR_THDPS:              pM_coa_c - pM_h2o_c + pM_sl2a6o_c - pM_succoa_c
                        - pM_thdp_c >= 0
 dR_EX_cobalt2_e:       - pM_cobalt2_e  = 0
 dR_COBALT2tex:         - pM_cobalt2_e + pM_cobalt2_p  = 0
 dR_CYTDH:              pM_csn_c - pM_cytd_c - pM_h2o_c + pM_rib_D_c >= 0
 dR_CSND:               - pM_csn_c - pM_h2o_c - pM_h_c + pM_nh4_c + pM_ura_c
                        >= 0
 dR_DCTPD2:             - pM_ctp_c - pM_h2o_c - pM_h_c + pM_nh4_c + pM_utp_c
                        >= 0
 dR_EX_cu2_e:           - pM_cu2_e  = 0
 dR_CU2tex:             - pM_cu2_e + pM_cu2_p  = 0
 dR_DM_cyanphy:         - pM_cyanphy_c >= 0
 dR_CYNTt2pp:           pM_cynt_c - pM_cynt_p + pM_h_c - pM_h_p >= 0
 dR_EX_cynt_e:          - pM_cynt_e >= 0
 dR_CYNTtex:            - pM_cynt_e + pM_cynt_p  = 0
 dR_CYSTL:              - pM_cyst_L_c - pM_h2o_c + pM_hcys_L_c + pM_nh4_c
                        + pM_pyr_c >= 0
 dR_DB4PS:              pM_db4p_c + pM_for_c + pM_h_c - pM_ru5p_D_c >= 0
 dR_EAR100y:            pM_dcaACP_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        - pM_tdec2eACP_c >= 0
 dR_DCTPD:              - pM_dctp_c + pM_dutp_c - pM_h2o_c - pM_h_c + pM_nh4_c
                        >= 0
 dR_EAR120y:            pM_ddcaACP_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        - pM_tddec2eACP_c >= 0
 dR_DPPS:               pM_decdp_c - pM_frdp_c - 7 pM_ipdp_c + 7 pM_ppi_c >= 0
 dR_UDPDPS:             - pM_decdp_c - pM_ipdp_c + pM_ppi_c + pM_udcpdp_c >= 0
 dR_DGDG160:            pM_dgdg160_c - pM_mgdg160_c + pM_udp_c - pM_udpgal_c
                        >= 0
 dR_DGDG161:            pM_dgdg161_c - pM_mgdg161_c + pM_udp_c - pM_udpgal_c
                        >= 0
 dR_DGDG180:            pM_dgdg180_c - pM_mgdg180_c + pM_udp_c - pM_udpgal_c
                        >= 0
 dR_DGDG181_9:          pM_dgdg181_9_c - pM_mgdg181_9_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_DGDG181:            pM_dgdg181_c - pM_mgdg181_c + pM_udp_c - pM_udpgal_c
                        >= 0
 dR_DGDG182_9_12:       pM_dgdg182_9_12_c - pM_mgdg182_9_12_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_DGDG183_6_9_12:     pM_dgdg183_6_9_12_c - pM_mgdg183_6_9_12_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_DGDG183_9_12_15:    pM_dgdg183_9_12_15_c - pM_mgdg183_9_12_15_c + pM_udp_c
                        - pM_udpgal_c >= 0
 dR_DGDG184_6_9_12_15:  pM_dgdg184_6_9_12_15_c - pM_mgdg184_6_9_12_15_c
                        + pM_udp_c - pM_udpgal_c >= 0
 dR_RNDR2:              pM_dgdp_c - pM_gdp_c + pM_h2o_c + pM_trdox_c
                        - pM_trdrd_c >= 0
 dR_NTPTP1:             pM_dgsn_c - pM_dgtp_c - pM_h2o_c + pM_pppi_c >= 0
 dR_GLYCDx:             pM_dha_c - pM_glyc_c + pM_h_c - pM_nad_c + pM_nadh_c
                        >= 0
 dR_FBA3:               pM_dhap_c + pM_e4p_c - pM_s17bp_c  = 0
 dR_FBA2:               pM_dhap_c - pM_f1p_c + pM_glyald_c  = 0
 dR_G3PD:               pM_dhap_c - pM_fad_c + pM_fadh2_c - pM_glyc3p_c >= 0
 dR_TPI:                - pM_dhap_c + pM_g3p_c  = 0
 dR_G3PD2:              pM_dhap_c - pM_glyc3p_c + pM_h_c - pM_nad_c + pM_nadh_c
                         = 0
 dR_MGSA:               - pM_dhap_c + pM_mthgxl_c + pM_pi_c >= 0
 dR_QULNS:              - pM_dhap_c + 2 pM_h2o_c - pM_iasp_c + pM_pi_c
                        + pM_quln_c >= 0
 dR_FBA:                pM_dhap_c - pM_fdp_c + pM_g3p_c  = 0
 dR_TMDS:               pM_dhf_c + pM_dtmp_c - pM_dump_c - pM_mlthf_c >= 0
 dR_DHFR:               - pM_dhf_c - pM_h_c + pM_nadp_c - pM_nadph_c + pM_thf_c
                         = 0
 dR_PDHcr:              - pM_dhlam_c + pM_h_c + pM_lpam_c - pM_nad_c
                        + pM_nadh_c  = 0
 dR_DNMPPA:             pM_dhnpt_c - pM_dhpmp_c - pM_h2o_c + pM_pi_c >= 0
 dR_DHORD:              - pM_dhor_S_c + pM_h2o2_c - pM_o2_c + pM_orot_c  = 0
 dR_DMATT:              - pM_dmpp_c + pM_grdp_c - pM_ipdp_c + pM_ppi_c >= 0
 dR_IPDDI:              pM_dmpp_c - pM_ipdp_c  = 0
 dR_DMPPS_syn:          pM_dmpp_c - pM_h2mb4p_c + pM_h2o_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_DM_dna5mtc:         - pM_dna5mtc_c >= 0
 dR_Sink_dna:           - pM_dna_c  = 0
 dR_SQLC2:              pM_dpterol_c - pM_h2o_c - pM_sql_c >= 0
 dR_SQLC:               pM_dptne_c - pM_sql_c >= 0
 dR_TDPDRE:             - pM_dtdp4d6dg_c + pM_dtdp4d6dm_c >= 0
 dR_TDPGDH:             pM_dtdp4d6dg_c - pM_dtdpglu_c + pM_h2o_c >= 0
 dR_TDPDRR:             - pM_dtdp4d6dm_c + pM_dtdprmn_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_MNXR6118:           pM_dtdpgal_c - pM_dtdpglu_c  = 0
 dR_G1PTT:              pM_dtdpglu_c - pM_dttp_c - pM_g1p_c - pM_h_c + pM_ppi_c
                        >= 0
 dR_RNDR4:              pM_dudp_c + pM_h2o_c + pM_trdox_c - pM_trdrd_c
                        - pM_udp_c >= 0
 dR_MPOMOR:             pM_dvpchlld_c + 2 pM_h2o_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c - pM_o2_c - pM_omppp9_c >= 0
 dR_MPOMOR2:            pM_dvpchlld_c + pM_h_c - pM_nadp_c + pM_nadph_c
                        - pM_omppp9_c >= 0
 dR_DVOCHR:             - pM_dvpchlld_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        + pM_pchlld_c >= 0
 dR_TALA:               pM_e4p_c + pM_f6p_c - pM_g3p_c - pM_s7p_c  = 0
 dR_TKT2:               - pM_e4p_c + pM_f6p_c + pM_g3p_c - pM_xu5p_D_c  = 0
 dR_IGPDH:              - pM_eig3p_c + pM_h2o_c + pM_imacp_c >= 0
 dR_ETOHtpp:            pM_etoh_c - pM_etoh_p  = 0
 dR_ETOHtex:            - pM_etoh_e + pM_etoh_p  = 0
 dR_EX_etoh_e:          - pM_etoh_e >= 0
 dR_GF6PTA:             - pM_f6p_c + pM_gam6p_c - pM_gln_L_c + pM_glu_L_c >= 0
 dR_PGI:                pM_f6p_c - pM_g6p_c  = 0
 dR_SPS:                - pM_f6p_c + pM_h_c + pM_suc6p_c + pM_udp_c - pM_udpg_c
                        >= 0
 dR_PGI2:               pM_f6p_c - pM_g6p_B_c  = 0
 dR_MAN6PI:             pM_f6p_c - pM_man6p_c  = 0
 dR_FBP:                pM_f6p_c - pM_fdp_c - pM_h2o_c + pM_pi_c >= 0
 dR_SUCD2pp:            pM_fad_c - pM_fadh2_c - pM_pq_p + pM_pqh2_p >= 0
 dR_SUCD2u:             pM_fad_c - pM_fadh2_c - pM_pq_u + pM_pqh2_u >= 0
 dR_SUCD4:              pM_fad_c - pM_fadh2_c - pM_q8_c + pM_q8h2_c  = 0
 dR_SUCD1:              - pM_fad_c + pM_fadh2_c + pM_fum_c - pM_succ_c  = 0
 dR_PSI_2:              - 2 pM_fdxo_2_2_c + 2 pM_fdxr_2_2_c + 2 pM_ficytc6_u
                        - 2 pM_focytc6_u - 2 pM_photon_c >= 0
 dR_DES_6B:             2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c + 2 pM_h2o_c
                        - 2 pM_h_c - pM_o2_c + pM_octe_6_9_12_ACP_c
                        - pM_octe_9_12_ACP_c >= 0
 dR_SULR_syn:           6 pM_fdxo_2_2_c - 6 pM_fdxr_2_2_c + 3 pM_h2o_c
                        + pM_h2s_c - 8 pM_h_c - pM_so3_c >= 0
 dR_FNOR:               2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c - pM_h_c - pM_nadp_c
                        + pM_nadph_c  = 0
 dR_FTR:                2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c - 2 pM_h_c
                        - pM_trdox_c + pM_trdrd_c >= 0
 dR_NOR_syn:            6 pM_fdxo_2_2_c - 6 pM_fdxr_2_2_c + 2 pM_h2o_c
                        - 8 pM_h_c + pM_nh4_c - pM_no2_c >= 0
 dR_NAR_syn:            2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c + pM_h2o_c - 2 pM_h_c
                        + pM_no2_c - pM_no3_c >= 0
 dR_FQR:                2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c - 2 pM_h_c - pM_pq_u
                        + pM_pqh2_u >= 0
 dR_DES_9:              2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c + 2 pM_h2o_c
                        - 2 pM_h_c - pM_o2_c - pM_ocdcaACP_c + pM_octe_9_ACP_c
                        >= 0
 dR_PSI:                - 2 pM_fdxo_2_2_c + 2 pM_fdxr_2_2_c + 2 pM_pcox_u
                        - 2 pM_pcrd_u - 2 pM_photon_c >= 0
 dR_DES_6:              2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c + 2 pM_h2o_c
                        - 2 pM_h_c - pM_o2_c + pM_octe_6_9_12_15_ACP_c
                        - pM_octe_9_12_15_ACP_c >= 0
 dR_DES_12:             2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c + 2 pM_h2o_c
                        - 2 pM_h_c - pM_o2_c + pM_octe_9_12_ACP_c
                        - pM_octe_9_ACP_c >= 0
 dR_DES_15:             2 pM_fdxo_2_2_c - 2 pM_fdxr_2_2_c + 2 pM_h2o_c
                        - 2 pM_h_c - pM_o2_c + pM_octe_9_12_15_ACP_c
                        - pM_octe_9_12_ACP_c >= 0
 dR_FCLT:               - pM_fe2_c + 2 pM_h_c + pM_pheme_c - pM_ppp9_c >= 0
 dR_FE2tex:             - pM_fe2_e + pM_fe2_p  = 0
 dR_EX_fe2_e:           - pM_fe2_e  = 0
 dR_EX_fe3_e:           - pM_fe3_e  = 0
 dR_FE3tex:             - pM_fe3_e + pM_fe3_p  = 0
 dR_MNXR15637:          pM_fgam_c - pM_gar_c - pM_h2o_c + 2 pM_h_c - pM_methf_c
                        + pM_thf_c >= 0
 dR_CBFC2pp:            - 2 pM_ficytc6_p + 2 pM_focytc6_p - 2 pM_h_c + 4 pM_h_p
                        + pM_pq_p - pM_pqh2_p >= 0
 dR_CYO1bpp_syn:        2 pM_ficytc6_p - 2 pM_focytc6_p + pM_h2o_p - 4 pM_h_c
                        + 2 pM_h_p - 0.5 pM_o2_p >= 0
 dR_CYO1b_syn:          2 pM_ficytc6_u - 2 pM_focytc6_u + pM_h2o_u - 4 pM_h_c
                        + 2 pM_h_u - 0.5 pM_o2_u >= 0
 dR_CBFC2:              - 2 pM_ficytc6_u + 2 pM_focytc6_u - 2 pM_h_c + 4 pM_h_u
                        + pM_pq_u - pM_pqh2_u >= 0
 dR_ACP1_FMN:           - pM_fmn_c - pM_h2o_c + pM_pi_c + pM_ribflv_c >= 0
 dR_IMPC:               pM_fprica_c - pM_h2o_c - pM_imp_c  = 0
 dR_SQLS:               - 2 pM_frdp_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        + 2 pM_ppi_c + pM_sql_c >= 0
 dR_FRTT:               - pM_frdp_c + pM_ggdp_c - pM_ipdp_c + pM_ppi_c >= 0
 dR_GRTT:               pM_frdp_c - pM_grdp_c - pM_ipdp_c + pM_ppi_c >= 0
 dR_HEMEOS:             - pM_frdp_c - pM_h2o_c + pM_hemeO_c - pM_pheme_c
                        + pM_ppi_c >= 0
 dR_FRUt3:              pM_fru_c - pM_fru_p + pM_h_c - pM_h_p >= 0
 dR_FRUtex:             - pM_fru_e + pM_fru_p  = 0
 dR_EX_fru_e:           - pM_fru_e >= 0
 dR_FUM:                - pM_fum_c - pM_h2o_c + pM_mal_L_c  = 0
 dR_FUMtpp:             pM_fum_c - pM_fum_p + pM_h_c - pM_h_p  = 0
 dR_FUMtex:             - pM_fum_e + pM_fum_p  = 0
 dR_EX_fum_e:           - pM_fum_e >= 0
 dR_UGLT:               pM_g1p_c - pM_gal1p_c - pM_udpg_c + pM_udpgal_c  = 0
 dR_GLCP:               pM_g1p_c - pM_glycogen_c - pM_pi_c >= 0
 dR_GALU:               - pM_g1p_c - pM_h_c + pM_ppi_c + pM_udpg_c - pM_utp_c
                         = 0
 dR_PGMT:               - pM_g1p_c + pM_g6p_c  = 0
 dR_TKT1:               pM_g3p_c - pM_r5p_c + pM_s7p_c - pM_xu5p_D_c  = 0
 dR_MI1PS:              - pM_g6p_c + pM_mi1p_D_c >= 0
 dR_PGAMT:              - pM_gam1p_c + pM_gam6p_c  = 0
 dR_GCALDD:             - pM_gcald_c + pM_glyclt_c - pM_h2o_c + 2 pM_h_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_LYCOPC:             pM_gcarote_c - pM_lycop_c >= 0
 dR_PYK3:               - pM_gdp_c + pM_gtp_c - pM_h_c - pM_pep_c + pM_pyr_c
                        >= 0
 dR_MAN1PT2:            - pM_gdp_c + pM_gdpmann_c - pM_h_c - pM_man1p_c
                        + pM_pi_c >= 0
 dR_GFUCS:              - pM_gdpddman_c + pM_gdpfuc_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c >= 0
 dR_GMAND:              pM_gdpddman_c - pM_gdpmann_c + pM_h2o_c >= 0
 dR_MAN1PT:             pM_gdpmann_c - pM_gtp_c - pM_h_c - pM_man1p_c
                        + pM_ppi_c >= 0
 dR_PHYTES:             - 2 pM_ggdp_c + pM_ppi_c + pM_prephytedp_c >= 0
 dR_GGDPR:              - pM_ggdp_c - 3 pM_h_c + 3 pM_nadp_c - 3 pM_nadph_c
                        + pM_phdp_c >= 0
 dR_MALT:               2 pM_glc_D_c - pM_h2o_c - pM_malt_c >= 0
 dR_GLCt2pp:            pM_glc_D_c - pM_glc_D_p + pM_h_c - pM_h_p  = 0
 dR_GLCtex:             - pM_glc_D_e + pM_glc_D_p  = 0
 dR_EX_glc_e:           - pM_glc_D_e >= 0
 dR_GGPP:               pM_glcglyc_c - pM_glcglycp_c - pM_h2o_c + pM_pi_c >= 0
 dR_GLCGLYCtex:         - pM_glcglyc_e + pM_glcglyc_p  = 0
 dR_EX_glcglyc_e:       - pM_glcglyc_e >= 0
 dR_GLUN:               - pM_gln_L_c + pM_glu_L_c - pM_h2o_c + pM_nh4_c >= 0
 dR_GLUPRT:             - pM_gln_L_c + pM_glu_L_c - pM_h2o_c + pM_ppi_c
                        + pM_pram_c - pM_prpp_c >= 0
 dR_GLNtex:             - pM_gln_L_e + pM_gln_L_p  = 0
 dR_EX_gln_L_e:         - pM_gln_L_e >= 0
 dR_GLUTRR:             pM_glu1sa_c - pM_glutrna_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c + pM_trnaglu_c >= 0
 dR_G5SD:               - pM_glu5p_c + pM_glu5sa_c - pM_h_c + pM_nadp_c
                        - pM_nadph_c + pM_pi_c >= 0
 dR_PUTA3:              - pM_glu5sa_c + pM_glu_L_c - pM_h2o_c + 2 pM_h_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_GLUR:               - pM_glu_D_c + pM_glu_L_c  = 0
 dR_GLUt4pp:            pM_glu_L_c - pM_glu_L_p + pM_na1_c - pM_na1_p >= 0
 dR_GLUttrappp:         pM_glu_L_c - pM_glu_L_p + pM_h_c - pM_h_p >= 0
 dR_GLUtex:             - pM_glu_L_e + pM_glu_L_p  = 0
 dR_EX_glu_L_e:         - pM_glu_L_e >= 0
 dR_GLYCTO1:            pM_glx_c - pM_glyclt_c + pM_h2o2_c - pM_o2_c >= 0
 dR_SPT_syn:            - pM_glx_c + pM_gly_c + pM_hpyr_c - pM_ser_L_c  = 0
 dR_GLYCTO_syn:         pM_glx_c - pM_glyclt_c + pM_h_c - pM_nad_c + pM_nadh_c
                        >= 0
 dR_GLXO3r:             - pM_glx_c + pM_h2o2_c - pM_h2o_c + pM_h_c - pM_o2_c
                        + pM_oxa_c  = 0
 dR_GHMT2r:             pM_gly_c + pM_h2o_c + pM_mlthf_c - pM_ser_L_c
                        - pM_thf_c  = 0
 dR_GLYtex:             - pM_gly_e + pM_gly_p  = 0
 dR_EX_gly_e:           - pM_gly_e >= 0
 dR_GLYALDDr:           - pM_glyald_c + pM_glyc_R_c - pM_h2o_c + 2 pM_h_c
                        - pM_nad_c + pM_nadh_c  = 0
 dR_ALCD19:             - pM_glyald_c + pM_glyc_c - pM_h_c + pM_nad_c
                        - pM_nadh_c  = 0
 dR_HPYRR2_syn:         pM_glyc_R_c - pM_h_c - pM_hpyr_c + pM_nad_c - pM_nadh_c
                        >= 0
 dR_Sink_glycogen:      - pM_glycogen_c >= 0
 dR_EX_gmp_c:           - pM_gmp_c >= 0
 dR_NTPP2:              pM_gmp_c - pM_gtp_c - pM_h2o_c + pM_h_c + pM_ppi_c >= 0
 dR_GTHP:               pM_gthox_c - 2 pM_gthrd_c - pM_h2o2_c + 2 pM_h2o_c  = 0
 dR_LALDO:              - pM_gthrd_c + pM_h_c - pM_lald_D_c + pM_lgt_S_c
                        - pM_nad_c + pM_nadh_c  = 0
 dR_LGTHL:              - pM_gthrd_c + pM_lgt_S_c - pM_mthgxl_c  = 0
 dR_GLYOX:              pM_gthrd_c - pM_h2o_c + pM_h_c + pM_lac_D_c
                        - pM_lgt_S_c >= 0
 dR_H2tpp:              pM_h2_c - pM_h2_p  = 0
 dR_H2ASE_syn:          pM_h2_c - pM_h_c + pM_nadp_c - pM_nadph_c  = 0
 dR_H2tex:              - pM_h2_e + pM_h2_p  = 0
 dR_EX_h2_e:            - pM_h2_e >= 0
 dR_IPDPS_syn:          - pM_h2mb4p_c + pM_h2o_c - pM_h_c + pM_ipdp_c
                        + pM_nadp_c - pM_nadph_c >= 0
 dR_PYDXNO:             pM_h2o2_c - pM_o2_c + pM_pydx_c - pM_pydxn_c  = 0
 dR_CAT:                - 2 pM_h2o2_c + 2 pM_h2o_c + pM_o2_c >= 0
 dR_PYDXO:              2 pM_h2o2_c - 2 pM_h2o_c - pM_nh4_c - 0.5 pM_o2_c
                        + pM_pydam_c - pM_pydx_c  = 0
 dR_PYAM5PO:            pM_h2o2_c - pM_h2o_c + pM_nh4_c - pM_o2_c - pM_pyam5p_c
                        + pM_pydx5p_c >= 0
 dR_PDX5PO:             pM_h2o2_c - pM_o2_c - pM_pdx5p_c + pM_pydx5p_c  = 0
 dR_HEMEAS:             - pM_h2o_c + 4 pM_h_c + pM_hemeA_c - pM_hemeO_c >= 0
 dR_SSALyr:             - pM_h2o_c + 2 pM_h_c - pM_nadp_c + pM_nadph_c
                        + pM_succ_c - pM_sucsal_c  = 0
 dR_PGPP182_9_12:       - pM_h2o_c + pM_pg182_9_12_c - pM_pgp182_9_12_c
                        + pM_pi_c >= 0
 dR_PGPP184_6_9_12_15:  - pM_h2o_c + pM_pg184_6_9_12_15_c
                        - pM_pgp184_6_9_12_15_c + pM_pi_c >= 0
 dR_IMACTD:             - pM_h2o_c + 2 pM_h_c + pM_im4ac_c - pM_im4act_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_MEHLER:             pM_h2o_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        - 0.5 pM_o2_c >= 0
 dR_PGPP183_9_12_15:    - pM_h2o_c + pM_pg183_9_12_15_c - pM_pgp183_9_12_15_c
                        + pM_pi_c >= 0
 dR_PRAMPC:             - pM_h2o_c - pM_prbamp_c + pM_prfp_c >= 0
 dR_USHD:               - pM_h2o_c + 2 pM_h_c + pM_lipidX_c - pM_u23ga_c
                        + pM_ump_c >= 0
 dR_PGPP181_9:          - pM_h2o_c + pM_pg181_9_c - pM_pgp181_9_c + pM_pi_c
                        >= 0
 dR_TRPS2:              pM_h2o_c - pM_indole_c - pM_ser_L_c + pM_trp_L_c  = 0
 dR_ZCAROTDH2:          2 pM_h2o_c - pM_h_c + pM_lycop_c + pM_nadp_c
                        - pM_nadph_c - pM_neuspn_c - pM_o2_c >= 0
 dR_SPP:                - pM_h2o_c + pM_pi_c - pM_suc6p_c + pM_sucr_c >= 0
 dR_SBP:                - pM_h2o_c + pM_pi_c - pM_s17bp_c + pM_s7p_c >= 0
 dR_IMPD:               - pM_h2o_c + pM_h_c - pM_imp_c - pM_nad_c + pM_nadh_c
                        + pM_xmp_c >= 0
 dR_NTPP8:              - pM_h2o_c + pM_h_c + pM_ppi_c + pM_ump_c - pM_utp_c
                        >= 0
 dR_NTPP9:              - pM_h2o_c + pM_h_c + pM_imp_c - pM_itp_c + pM_ppi_c
                        >= 0
 dR_SQD1:               pM_h2o_c - pM_h_c - pM_so3_c - pM_udpg_c + pM_udpsq_c
                        >= 0
 dR_ZCAROTDH1:          2 pM_h2o_c - pM_h_c + pM_nadp_c - pM_nadph_c
                        + pM_neuspn_c - pM_o2_c - pM_zcarote_c >= 0
 dR_PPA:                - pM_h2o_c + pM_h_c + 2 pM_pi_c - pM_ppi_c >= 0
 dR_HMBS:               - pM_h2o_c + pM_hmbil_c + 4 pM_nh4_c - 4 pM_ppbng_c
                        >= 0
 dR_THRS:               - pM_h2o_c - pM_phom_c + pM_pi_c + pM_thr_L_c >= 0
 dR_AMID3:              - pM_h2o_c - pM_iad_c + pM_ind3ac_c + pM_nh4_c >= 0
 dR_AMID2:              - pM_h2o_c + pM_nh4_c + pM_pac_c - pM_pad_c >= 0
 dR_ALDD20x:            - pM_h2o_c + 2 pM_h_c - pM_id3acald_c + pM_ind3ac_c
                        - pM_nad_c + pM_nadh_c >= 0
 dR_MPOMC2:             - pM_h2o_c + pM_h_c + pM_hmppp9_c - pM_mppp9om_c
                        - pM_nadp_c + pM_nadph_c >= 0
 dR_MPOMC1:             pM_h2o_c - pM_h_c + pM_hmppp9_c - pM_mppp9om_c
                        + pM_nadp_c - pM_nadph_c - pM_o2_c >= 0
 dR_UDCPDP:             - pM_h2o_c + pM_h_c + pM_pi_c - pM_udcpdp_c
                        + pM_udcpp_c >= 0
 dR_UPP3S:              pM_h2o_c - pM_hmbil_c + pM_uppg3_c >= 0
 dR_PPPGO:              3 pM_h2o_c - 1.5 pM_o2_c + pM_ppp9_c - pM_pppg9_c >= 0
 dR_MI1PP:              - pM_h2o_c + pM_inost_c - pM_mi1p_D_c + pM_pi_c >= 0
 dR_MPOMMM:             2 pM_h2o_c - pM_h_c - pM_hmppp9_c + pM_nadp_c
                        - pM_nadph_c - pM_o2_c + pM_omppp9_c >= 0
 dR_PGPP160:            - pM_h2o_c + pM_pg160_c - pM_pgp160_c + pM_pi_c >= 0
 dR_PGPP161:            - pM_h2o_c + pM_pg161_c - pM_pgp161_c + pM_pi_c >= 0
 dR_PGPP180:            - pM_h2o_c + pM_pg180_c - pM_pgp180_c + pM_pi_c >= 0
 dR_PGPP181:            - pM_h2o_c + pM_pg181_c - pM_pgp181_c + pM_pi_c >= 0
 dR_HISTP:              - pM_h2o_c - pM_hisp_c + pM_histd_c + pM_pi_c >= 0
 dR_HISTD:              - pM_h2o_c + 3 pM_h_c + pM_his_L_c - pM_histd_c
                        - 2 pM_nad_c + 2 pM_nadh_c >= 0
 dR_H2Otpp:             pM_h2o_c - pM_h2o_p  = 0
 dR_PGPP183_6_9_12:     - pM_h2o_c + pM_pg183_6_9_12_c - pM_pgp183_6_9_12_c
                        + pM_pi_c >= 0
 dR_H2Otu_syn:          pM_h2o_c - pM_h2o_u  = 0
 dR_PRATPP:             - pM_h2o_c + pM_h_c + pM_ppi_c + pM_prbamp_c
                        - pM_prbatp_c >= 0
 dR_H2Otex:             - pM_h2o_e + pM_h2o_p  = 0
 dR_EX_h2o_e:           - pM_h2o_e  = 0
 dR_CYTBDpp:            pM_h2o_p - 0.5 pM_o2_p + pM_pq_p - pM_pqh2_p >= 0
 dR_CYO1b2pp_syn:       pM_h2o_p - 4 pM_h_c + 2 pM_h_p - 0.5 pM_o2_p
                        + 2 pM_pcox_p - 2 pM_pcrd_p >= 0
 dR_CYTBDu:             pM_h2o_u - 0.5 pM_o2_u + pM_pq_u - pM_pqh2_u >= 0
 dR_CYO1b2_syn:         pM_h2o_u - 4 pM_h_c + 2 pM_h_u - 0.5 pM_o2_u
                        + 2 pM_pcox_u - 2 pM_pcrd_u >= 0
 dR_PSII:               - pM_h2o_u - 2 pM_h_c + 2 pM_h_u + 0.5 pM_o2_u
                        - 2 pM_photon_c - pM_pq_u + pM_pqh2_u >= 0
 dR_PPTGS_Syn:          pM_h_c + pM_peptido_syn_c - pM_uaagmda_c + pM_udcpdp_c
                        >= 0
 dR_EAR141y:            - pM_h_c + pM_nadp_c - pM_nadph_c - pM_t3c7mrseACP_c
                        + pM_tdeACP_c >= 0
 dR_EAR160y:            - pM_h_c + pM_nadp_c - pM_nadph_c + pM_palmACP_c
                        - pM_tpalm2eACP_c >= 0
 dR_PAPSR:              2 pM_h_c + pM_pap_c - pM_paps_c + pM_so3_c + pM_trdox_c
                        - pM_trdrd_c >= 0
 dR_PYRtpp:             pM_h_c - pM_h_p + pM_pyr_c - pM_pyr_p  = 0
 dR_NAt3pp:             pM_h_c - pM_h_p - pM_na1_c + pM_na1_p  = 0
 dR_UAGPT3:             pM_h_c + pM_uaagmda_c - pM_uacgam_c - pM_uagmda_c
                        + pM_udp_c >= 0
 dR_MNt2pp:             pM_h_c - pM_h_p + pM_mn2_c - pM_mn2_p >= 0
 dR_SPMDt2pp:           pM_h_c - pM_h_p + pM_spmd_c - pM_spmd_p >= 0
 dR_NDH1_1u:            - 4 pM_h_c + 3 pM_h_u + pM_nadp_c - pM_nadph_c
                        - pM_pq_u + pM_pqh2_u >= 0
 dR_NDH1_1p:            - 4 pM_h_c + 3 pM_h_p + pM_nadp_c - pM_nadph_c
                        - pM_pq_p + pM_pqh2_p >= 0
 dR_PTRCt2pp:           pM_h_c - pM_h_p + pM_ptrc_c - pM_ptrc_p >= 0
 dR_PHYTFDH2:           pM_h_c - pM_nadp_c + pM_nadph_c - pM_phytof_c
                        + pM_zcarote_c >= 0
 dR_PHYTFDH1:           pM_h_c - pM_nad_c + pM_nadh_c - pM_phytof_c
                        + pM_zcarote_c >= 0
 dR_EAR60y:             - pM_h_c + pM_hexACP_c + pM_nadp_c - pM_nadph_c
                        - pM_thex2eACP_c >= 0
 dR_NH4tpp:             pM_h_c - pM_h_p + pM_nh4_c - pM_nh4_p >= 0
 dR_NMNS:               - pM_h_c - pM_ncam_c + pM_nmn_c + pM_ppi_c - pM_prpp_c
                        >= 0
 dR_PYK5:               - pM_h_c - pM_idp_c + pM_itp_c - pM_pep_c + pM_pyr_c
                        >= 0
 dR_PYK2:               - pM_h_c - pM_pep_c + pM_pyr_c - pM_udp_c + pM_utp_c
                        >= 0
 dR_Ktpp:               pM_h_c - pM_h_p + pM_k_c - pM_k_p >= 0
 dR_PC6AR:              - pM_h_c + pM_nadp_c - pM_nadph_c - pM_pre6a_c
                        + pM_pre6b_c >= 0
 dR_UAPGR:              - pM_h_c + pM_nadp_c - pM_nadph_c - pM_uaccg_c
                        + pM_uamr_c >= 0
 dR_TRDR:               - pM_h_c + pM_nadp_c - pM_nadph_c - pM_trdox_c
                        + pM_trdrd_c >= 0
 dR_R05217:             - pM_h_c - 0.5 pM_o2_c - pM_pre3a_c + pM_pre3b_c >= 0
 dR_LDH_D:              pM_h_c - pM_lac_D_c - pM_nad_c + pM_nadh_c + pM_pyr_c
                         = 0
 dR_EAR161y:            - pM_h_c + pM_hdeACP_c + pM_nadp_c - pM_nadph_c
                        - pM_t3c9palmeACP_c >= 0
 dR_SUCCtpp:            pM_h_c - pM_h_p + pM_succ_c - pM_succ_p  = 0
 dR_NDH1_2u:            - 4 pM_h_c + 3 pM_h_u + pM_nad_c - pM_nadh_c - pM_pq_u
                        + pM_pqh2_u >= 0
 dR_NDH2_pp:            - pM_h_c + pM_nad_c - pM_nadh_c - pM_pq_p + pM_pqh2_p
                        >= 0
 dR_EAR180y:            - pM_h_c + pM_nadp_c - pM_nadph_c + pM_ocdcaACP_c
                        - pM_toctd2eACP_c >= 0
 dR_MALtpp:             pM_h_c - pM_h_p + pM_mal_L_c - pM_mal_L_p  = 0
 dR_LPADSS:             pM_h_c + pM_lipidAds_c - pM_lipidX_c - pM_u23ga_c
                        + pM_udp_c >= 0
 dR_EAR140y:            - pM_h_c + pM_myrsACP_c + pM_nadp_c - pM_nadph_c
                        - pM_tmrs2eACP_c >= 0
 dR_EAR181y:            - pM_h_c + pM_nadp_c - pM_nadph_c + pM_octeACP_c
                        - pM_t3c11vaceACP_c >= 0
 dR_MDH:                pM_h_c - pM_mal_L_c - pM_nad_c + pM_nadh_c + pM_oaa_c
                         = 0
 dR_CBFCpp:             - 2 pM_h_c + 4 pM_h_p - 2 pM_pcox_p + 2 pM_pcrd_p
                        + pM_pq_p - pM_pqh2_p >= 0
 dR_PHYTEDH1:           pM_h_c - pM_nad_c + pM_nadh_c - pM_phytoe_c
                        + pM_phytof_c >= 0
 dR_NDH2_syn:           - pM_h_c + pM_nad_c - pM_nadh_c - pM_pq_u + pM_pqh2_u
                        >= 0
 dR_NADH5:              - pM_h_c + pM_nad_c - pM_nadh_c - pM_q8_c + pM_q8h2_c
                        >= 0
 dR_CBFCu:              - 2 pM_h_c + 4 pM_h_u - 2 pM_pcox_u + 2 pM_pcrd_u
                        + pM_pq_u - pM_pqh2_u >= 0
 dR_MPOMMM2:            pM_h_c - pM_hmppp9_c - pM_nadp_c + pM_nadph_c
                        + pM_omppp9_c >= 0
 dR_PC8XM:              - pM_h_c + pM_hgbyr_c - pM_pre8_c >= 0
 dR_NDH1_2p:            - 4 pM_h_c + 3 pM_h_p + pM_nad_c - pM_nadh_c - pM_pq_p
                        + pM_pqh2_p >= 0
 dR_EAR80y:             - pM_h_c + pM_nadp_c - pM_nadph_c + pM_ocACP_c
                        - pM_toct2eACP_c >= 0
 dR_PHYTEDH2:           pM_h_c - pM_nadp_c + pM_nadph_c - pM_phytoe_c
                        + pM_phytof_c >= 0
 dR_Htex:               - pM_h_e + pM_h_p  = 0
 dR_EX_h_e:             - pM_h_e  = 0
 dR_H2CO3_NAt_syn:      pM_hco3_c - pM_hco3_p + pM_na1_c - pM_na1_p >= 0
 dR_EX_hco3_e:          - pM_hco3_e - l832  = 0
 dR_HCO3tex:            - pM_hco3_e + pM_hco3_p  = 0
 dR_EX_his_L_e:         - pM_his_L_e >= 0
 dR_HIStex:             - pM_his_L_e + pM_his_L_p  = 0
 dR_EX_k_e:             - pM_k_e  = 0
 dR_Kex:                - pM_k_e + pM_k_p  = 0
 dR_LEUtex:             - pM_leu_L_e + pM_leu_L_p  = 0
 dR_EX_leu_L_e:         - pM_leu_L_e >= 0
 dR_EX_lys_L_e:         - pM_lys_L_e >= 0
 dR_LYStex:             - pM_lys_L_e + pM_lys_L_p  = 0
 dR_MALtex:             - pM_mal_L_e + pM_mal_L_p  = 0
 dR_EX_mal_L_e:         - pM_mal_L_e >= 0
 dR_PMANM:              - pM_man1p_c + pM_man6p_c  = 0
 dR_MTHFD:              pM_methf_c - pM_mlthf_c - pM_nadp_c + pM_nadph_c >= 0
 dR_MTHFD2i:            pM_methf_c - pM_mlthf_c - pM_nad_c + pM_nadh_c >= 0
 dR_EX_mg2_e:           - pM_mg2_e  = 0
 dR_MG2tex:             - pM_mg2_e + pM_mg2_p  = 0
 dR_EX_mn2_e:           - pM_mn2_e  = 0
 dR_MNtex:              - pM_mn2_e + pM_mn2_p  = 0
 dR_MOBDtex:            - pM_mobd_e + pM_mobd_p  = 0
 dR_EX_mobd_e:          - pM_mobd_e  = 0
 dR_EX_na1_e:           - pM_na1_e  = 0
 dR_NAtex:              - pM_na1_e + pM_na1_p  = 0
 dR_NADTRHD:            - pM_nad_c + pM_nadh_c + pM_nadp_c - pM_nadph_c >= 0
 dR_SER_AL:             pM_nh4_c + pM_pyr_c - pM_ser_L_c >= 0
 dR_EX_nh4_e:           - pM_nh4_e >= 0
 dR_NH4tex:             - pM_nh4_e + pM_nh4_p  = 0
 dR_Ni2tex:             - pM_ni2_e + pM_ni2_p  = 0
 dR_EX_ni2_e:           - pM_ni2_e  = 0
 dR_EX_no3_e:           - pM_no3_e  = 0
 dR_NO3tex:             - pM_no3_e + pM_no3_p  = 0
 dR_O2tu:               pM_o2_c - pM_o2_u  = 0
 dR_O2tc:               pM_o2_c - pM_o2_p  = 0
 dR_EX_o2_e:            - pM_o2_e  = 0
 dR_O2tex:              - pM_o2_e + pM_o2_p  = 0
 dR_ORPT:               - pM_orot5p_c + pM_orot_c - pM_ppi_c + pM_prpp_c  = 0
 dR_UAGCVT:             - pM_pep_c + pM_pi_c + pM_uaccg_c - pM_uacgam_c >= 0
 dR_PPTGF_Syn:          - pM_peptido_syn_c + pM_peptido_syn_p >= 0
 dR_Sink_phbg:          - pM_phbg_c  = 0
 dR_PHOtpp:             pM_photon_c - pM_photon_p  = 0
 dR_PHOtex:             - pM_photon_e + pM_photon_p  = 0
 dR_EX_photon_e:        - pM_photon_e - l873 + u873  = 0
 dR_PHYTES2:            pM_phytoe_c + pM_ppi_c - pM_prephytedp_c >= 0
 dR_EX_pi_e:            - pM_pi_e  = 0
 dR_PItex:              - pM_pi_e + pM_pi_p  = 0
 dR_UPPRT:              pM_ppi_c - pM_prpp_c + pM_ump_c - pM_ura_c >= 0
 dR_Sink_precyanphy:    - pM_precyanphy_c  = 0
 dR_PRMICIi:            - pM_prfp_c + pM_prlp_c >= 0
 dR_PROtex:             - pM_pro_L_e + pM_pro_L_p  = 0
 dR_EX_pro_L_e:         - pM_pro_L_e >= 0
 dR_PTRCtex:            - pM_ptrc_e + pM_ptrc_p  = 0
 dR_EX_ptrc_e:          - pM_ptrc_e >= 0
 dR_PYRtex:             - pM_pyr_e + pM_pyr_p  = 0
 dR_EX_pyr_e:           - pM_pyr_e >= 0
 dR_PPM:                - pM_r1p_c + pM_r5p_c  = 0
 dR_RPI:                - pM_r5p_c + pM_ru5p_D_c  = 0
 dR_Sink_rdmbzi:        - pM_rdmbzi_c  = 0
 dR_RPE:                - pM_ru5p_D_c + pM_xu5p_D_c  = 0
 dR_EX_ser_L_e:         - pM_ser_L_e >= 0
 dR_SERtex:             - pM_ser_L_e + pM_ser_L_p  = 0
 dR_EX_so4_e:           - pM_so4_e  = 0
 dR_SO4tex:             - pM_so4_e + pM_so4_p  = 0
 dR_EX_spmd_e:          - pM_spmd_e >= 0
 dR_SPMDtex:            - pM_spmd_e + pM_spmd_p  = 0
 dR_EX_succ_e:          - pM_succ_e >= 0
 dR_SUCCtex:            - pM_succ_e + pM_succ_p  = 0
 dR_SUCRtex:            - pM_sucr_e + pM_sucr_p  = 0
 dR_EX_sucr_e:          - pM_sucr_e >= 0
 dR_PAPPT3:             pM_uagmda_c - pM_udcpp_c - pM_ugmda_c + pM_ump_c >= 0
 dR_UDPG4E:             - pM_udpg_c + pM_udpgal_c  = 0
 dR_UREAtex:            - pM_urea_e + pM_urea_p  = 0
 dR_EX_urea_e:          - pM_urea_e >= 0
 dR_EX_zn2_e:           - pM_zn2_e  = 0
 dR_Zn2tex:             - pM_zn2_e + pM_zn2_p  = 0
Bounds
      pM_10fthf_c Free
      pM_12dgr160_c Free
      pM_12dgr161_c Free
      pM_12dgr180_c Free
      pM_12dgr181_9_c Free
      pM_12dgr181_c Free
      pM_12dgr182_9_12_c Free
      pM_12dgr183_6_9_12_c Free
      pM_12dgr183_9_12_15_c Free
      pM_12dgr184_6_9_12_15_c Free
      pM_12ppd_S_c Free
      pM_13dpg_c Free
      pM_15dap_c Free
      pM_1hdec9eg3p_c Free
      pM_1hdecg3p_c Free
      pM_1odec11eg3p_c Free
      pM_1odec691215eg3p_c Free
      pM_1odec6912eg3p_c Free
      pM_1odec91215eg3p_c Free
      pM_1odec912eg3p_c Free
      pM_1odec9eg3p_c Free
      pM_1odecg3p_c Free
      pM_1p3h5c_c Free
      pM_1pyr5c_c Free
      pM_23DB_c Free
      pM_23DB_e Free
      pM_23DB_p Free
      pM_23dhdp_c Free
      pM_23dhmb_c Free
      pM_23dhmp_c Free
      pM_23dmphol_c Free
      pM_23dpg_c Free
      pM_25aics_c Free
      pM_25dhpp_c Free
      pM_26dap_LL_c Free
      pM_26dap_M_c Free
      pM_2ahbut_c Free
      pM_2ahhmd_c Free
      pM_2ahhmp_c Free
      pM_2cpr5p_c Free
      pM_2dda7p_c Free
      pM_2dhp_c Free
      pM_2dr5p_c Free
      pM_2h3oppan_c Free
      pM_2hyoxplac_c Free
      pM_2ippm_c Free
      pM_2m6phol_c Free
      pM_2mahmp_c Free
      pM_2me4p_c Free
      pM_2mecdp_c Free
      pM_2obut_c Free
      pM_2oph_c Free
      pM_2p4c2me_c Free
      pM_2pg_c Free
      pM_2pglyc_c Free
      pM_2shchc_c Free
      pM_34dhmald_c Free
      pM_34dhoxpeg_c Free
      pM_34hpp_c Free
      pM_35cgmp_c Free
      pM_3c2hmp_c Free
      pM_3c3hmp_c Free
      pM_3c4mop_c Free
      pM_3dhq_c Free
      pM_3dhsk_c Free
      pM_3haACP_c Free
      pM_3hbcoa_R_c Free
      pM_3hcddec5eACP_c Free
      pM_3hcmrs7eACP_c Free
      pM_3hcpalm9eACP_c Free
      pM_3hcvac11eACP_c Free
      pM_3hddecACP_c Free
      pM_3hdecACP_c Free
      pM_3hhexACP_c Free
      pM_3hmp_c Free
      pM_3hmrsACP_c Free
      pM_3hoctACP_c Free
      pM_3hoctaACP_c Free
      pM_3hpalmACP_c Free
      pM_3ig3p_c Free
      pM_3mob_c Free
      pM_3mop_c Free
      pM_3ocddec5eACP_c Free
      pM_3ocmrs7eACP_c Free
      pM_3ocpalm9eACP_c Free
      pM_3ocvac11eACP_c Free
      pM_3oddecACP_c Free
      pM_3odecACP_c Free
      pM_3ohexACP_c Free
      pM_3omrsACP_c Free
      pM_3ooctACP_c Free
      pM_3ooctdACP_c Free
      pM_3opalmACP_c Free
      pM_3ophb_c Free
      pM_3pg_c Free
      pM_3php_c Free
      pM_3psme_c Free
      pM_4aabutn_c Free
      pM_4abut_c Free
      pM_4abutn_c Free
      pM_4abz_c Free
      pM_4adcho_c Free
      pM_4ahmmp_c Free
      pM_4c2me_c Free
      pM_4gudbd_c Free
      pM_4gudbutn_c Free
      pM_4h2oglt_c Free
      pM_4hba_c Free
      pM_4hpro_LT_c Free
      pM_4mhetz_c Free
      pM_4mop_c Free
      pM_4mpetz_c Free
      pM_4pasp_c Free
      pM_4ppan_c Free
      pM_4ppcys_c Free
      pM_4r5au_c Free
      pM_5aizc_c Free
      pM_5aop_c Free
      pM_5aprbu_c Free
      pM_5apru_c Free
      pM_5fthf_c Free
      pM_5g2oxpt_c Free
      pM_5mdr1p_c Free
      pM_5mdru1p_c Free
      pM_5mta_c Free
      pM_5mthf_c Free
      pM_6pgc_c Free
      pM_6pgl_c Free
      pM_6pthp_c Free
      pM_8aonn_c Free
      pM_ACP_c Free
      pM_C04051_c Free
      pM_PHB_c Free
      pM_aacoa_c Free
      pM_acACP_c Free
      pM_ac_c Free
      pM_ac_e Free
      pM_ac_p Free
      pM_acald_c Free
      pM_accoa_c Free
      pM_acg5p_c Free
      pM_acg5sa_c Free
      pM_acgam1p_c Free
      pM_acglu_c Free
      pM_achms_c Free
      pM_acmana_c Free
      pM_acon_C_c Free
      pM_acorn_c Free
      pM_acser_c Free
      pM_actACP_c Free
      pM_actn_R_c Free
      pM_actp_c Free
      pM_adcobdam_c Free
      pM_adcobhex_c Free
      pM_ade_c Free
      pM_adhlam_c Free
      pM_adn_c Free
      pM_adocbi_c Free
      pM_adocbip_c Free
      pM_adocbl_c Free
      pM_adp_c Free
      pM_adpglc_c Free
      pM_adprib_c Free
      pM_agdpcbi_c Free
      pM_agm_c Free
      pM_ahcys_c Free
      pM_ahdt_c Free
      pM_aicar_c Free
      pM_air_c Free
      pM_akg_c Free
      pM_akg_e Free
      pM_akg_p Free
      pM_ala_B_c Free
      pM_ala_D_c Free
      pM_ala_L_c Free
      pM_ala_L_e Free
      pM_ala_L_p Free
      pM_alaala_c Free
      pM_alac_S_c Free
      pM_alatrna_c Free
      pM_amet_c Free
      pM_ametam_c Free
      pM_amp_c Free
      pM_anth_c Free
      pM_appl_c Free
      pM_applp_c Free
      pM_aps_c Free
      pM_arg_L_c Free
      pM_arg_L_e Free
      pM_arg_L_p Free
      pM_argsuc_c Free
      pM_argtrna_c Free
      pM_asn_L_c Free
      pM_asntrna_c Free
      pM_asp_L_c Free
      pM_aspsa_c Free
      pM_asptrna_c Free
      pM_atocophe_c Free
      pM_atp_c Free
      pM_bamppald_c Free
      pM_bcarote_c Free
      pM_bcryptox_c Free
      pM_bglycogen_c Free
      pM_bilirub_c Free
      pM_biliverd_c Free
      pM_biomass_dna_c Free
      pM_biomass_glycogen_c Free
      pM_biomass_ion_c Free
      pM_biomass_lipid_c Free
      pM_biomass_lps_c Free
      pM_biomass_pept_p Free
      pM_biomass_pigment_c Free
      pM_biomass_protein_c Free
      pM_biomass_rna_c Free
      pM_biomass_solublep_c Free
      pM_btn_c Free
      pM_btocophe_c Free
      pM_but2eACP_c Free
      pM_butACP_c Free
      pM_ca2_c Free
      pM_ca2_e Free
      pM_ca2_p Free
      pM_camp_c Free
      pM_cbasp_c Free
      pM_cbp_c Free
      pM_cddec5eACP_c Free
      pM_cdec3eACP_c Free
      pM_cdp_c Free
      pM_cdpdhdec9eg_c Free
      pM_cdpdhdecg_c Free
      pM_cdpdodec11eg_c Free
      pM_cdpdodec691215eg_c Free
      pM_cdpdodec6912eg_c Free
      pM_cdpdodec91215eg_c Free
      pM_cdpdodec912eg_c Free
      pM_cdpdodec9eg_c Free
      pM_cdpdodecg_c Free
      pM_cdpglc_c Free
      pM_cgly_c Free
      pM_chlld_c Free
      pM_cholphya_c Free
      pM_chor_c Free
      pM_cit_c Free
      pM_cit_e Free
      pM_cit_p Free
      pM_citr_L_c Free
      pM_cmp_c Free
      pM_co1dam_c Free
      pM_co2_c Free
      pM_co2_e Free
      pM_co2_p Free
      pM_co2dam_c Free
      pM_co_c Free
      pM_co_e Free
      pM_co_p Free
      pM_coa_c Free
      pM_cobalt2_c Free
      pM_cobalt2_e Free
      pM_cobalt2_p Free
      pM_cpppg1_c Free
      pM_cpppg3_c Free
      pM_csn_c Free
      pM_ctp_c Free
      pM_cu2_c Free
      pM_cu2_e Free
      pM_cu2_p Free
      pM_cu2_u Free
      pM_cyanphy_c Free
      pM_cynt_c Free
      pM_cynt_e Free
      pM_cynt_p Free
      pM_cys_L_c Free
      pM_cyst_L_c Free
      pM_cystrna_c Free
      pM_cytd_c Free
      pM_dad_5_c Free
      pM_dadp_c Free
      pM_damp_c Free
      pM_dann_c Free
      pM_datp_c Free
      pM_db4p_c Free
      pM_dcaACP_c Free
      pM_dcamp_c Free
      pM_dcdp_c Free
      pM_dctp_c Free
      pM_ddcaACP_c Free
      pM_decdp_c Free
      pM_dgdg160_c Free
      pM_dgdg161_c Free
      pM_dgdg180_c Free
      pM_dgdg181_9_c Free
      pM_dgdg181_c Free
      pM_dgdg182_9_12_c Free
      pM_dgdg183_6_9_12_c Free
      pM_dgdg183_9_12_15_c Free
      pM_dgdg184_6_9_12_15_c Free
      pM_dgdp_c Free
      pM_dgsn_c Free
      pM_dgtp_c Free
      pM_dha_c Free
      pM_dhap_c Free
      pM_dhf_c Free
      pM_dhlam_c Free
      pM_dhna_c Free
      pM_dhnpt_c Free
      pM_dhor_S_c Free
      pM_dhpmp_c Free
      pM_dhpt_c Free
      pM_didp_c Free
      pM_ditp_c Free
      pM_dmlz_c Free
      pM_dmpp_c Free
      pM_dmtphllqne_c Free
      pM_dna5mtc_c Free
      pM_dna_c Free
      pM_dnad_c Free
      pM_dpcoa_c Free
      pM_dpterol_c Free
      pM_dptne_c Free
      pM_drib_c Free
      pM_dscl_c Free
      pM_dtbt_c Free
      pM_dtdp4d6dg_c Free
      pM_dtdp4d6dm_c Free
      pM_dtdp_c Free
      pM_dtdpgal_c Free
      pM_dtdpglu_c Free
      pM_dtdprmn_c Free
      pM_dtmp_c Free
      pM_dtocophe_c Free
      pM_dttp_c Free
      pM_dudp_c Free
      pM_dump_c Free
      pM_dutp_c Free
      pM_dvpchlld_c Free
      pM_dxyl5p_c Free
      pM_dxyl_c Free
      pM_e4hglu_c Free
      pM_e4p_c Free
      pM_echin_c Free
      pM_eig3p_c Free
      pM_etoh_c Free
      pM_etoh_e Free
      pM_etoh_p Free
      pM_f1p_c Free
      pM_f6p_c Free
      pM_fad_c Free
      pM_fadh2_c Free
      pM_fdp_c Free
      pM_fdxo_2_2_c Free
      pM_fdxr_2_2_c Free
      pM_fe2_c Free
      pM_fe2_e Free
      pM_fe2_p Free
      pM_fe3_c Free
      pM_fe3_e Free
      pM_fe3_p Free
      pM_fgam_c Free
      pM_ficytc6_p Free
      pM_ficytc6_u Free
      pM_fmettrna_c Free
      pM_fmn_c Free
      pM_focytc6_p Free
      pM_focytc6_u Free
      pM_for_c Free
      pM_fpram_c Free
      pM_fprica_c Free
      pM_frdp_c Free
      pM_fru_c Free
      pM_fru_e Free
      pM_fru_p Free
      pM_fum_c Free
      pM_fum_e Free
      pM_fum_p Free
      pM_g1p_c Free
      pM_g3p_c Free
      pM_g6p_B_c Free
      pM_g6p_c Free
      pM_gal1p_c Free
      pM_gal_c Free
      pM_gam1p_c Free
      pM_gam6p_c Free
      pM_gar_c Free
      pM_gcald_c Free
      pM_gcarote_c Free
      pM_gdbtal_c Free
      pM_gdp_c Free
      pM_gdpddman_c Free
      pM_gdpfuc_c Free
      pM_gdpmann_c Free
      pM_ggdp_c Free
      pM_glc_B_c Free
      pM_glc_D_c Free
      pM_glc_D_e Free
      pM_glc_D_p Free
      pM_glcglyc_c Free
      pM_glcglyc_e Free
      pM_glcglyc_p Free
      pM_glcglycp_c Free
      pM_gln_L_c Free
      pM_gln_L_e Free
      pM_gln_L_p Free
      pM_glntrna_c Free
      pM_glu1sa_c Free
      pM_glu5p_c Free
      pM_glu5sa_c Free
      pM_glu_D_c Free
      pM_glu_L_c Free
      pM_glu_L_e Free
      pM_glu_L_p Free
      pM_glucys_c Free
      pM_glutrna_c Free
      pM_glutrna_gln_c Free
      pM_glx_c Free
      pM_gly_c Free
      pM_gly_e Free
      pM_gly_p Free
      pM_glyald_c Free
      pM_glyc3p_c Free
      pM_glyc_R_c Free
      pM_glyc_c Free
      pM_glyclt_c Free
      pM_glycogen_c Free
      pM_glytrna_c Free
      pM_gmp_c Free
      pM_grdp_c Free
      pM_gthox_c Free
      pM_gthrd_c Free
      pM_gtocophe_c Free
      pM_gtp_c Free
      pM_h2_c Free
      pM_h2_e Free
      pM_h2_p Free
      pM_h2mb4p_c Free
      pM_h2o2_c Free
      pM_h2o_c Free
      pM_h2o_e Free
      pM_h2o_p Free
      pM_h2o_u Free
      pM_h2s_c Free
      pM_h_c Free
      pM_h_e Free
      pM_h_p Free
      pM_h_u Free
      pM_hco3_c Free
      pM_hco3_e Free
      pM_hco3_p Free
      pM_hcys_L_c Free
      pM_hdeACP_c Free
      pM_hemeA_c Free
      pM_hemeO_c Free
      pM_hexACP_c Free
      pM_hgbam_c Free
      pM_hgbyr_c Free
      pM_hgentis_c Free
      pM_his_L_c Free
      pM_his_L_e Free
      pM_his_L_p Free
      pM_hisp_c Free
      pM_histd_c Free
      pM_histrna_c Free
      pM_hmbil_c Free
      pM_hmppp9_c Free
      pM_hom_L_c Free
      pM_hpyr_c Free
      pM_iad_c Free
      pM_iasp_c Free
      pM_ichor_c Free
      pM_icit_c Free
      pM_id3acald_c Free
      pM_idp_c Free
      pM_ile_L_c Free
      pM_iletrna_c Free
      pM_im4ac_c Free
      pM_im4act_c Free
      pM_imacp_c Free
      pM_imp_c Free
      pM_ind3ac_c Free
      pM_indole_c Free
      pM_inost_c Free
      pM_ipdp_c Free
      pM_itp_c Free
      pM_k_c Free
      pM_k_e Free
      pM_k_p Free
      pM_lac_D_c Free
      pM_lald_D_c Free
      pM_lald_L_c Free
      pM_leu_L_c Free
      pM_leu_L_e Free
      pM_leu_L_p Free
      pM_leutrna_c Free
      pM_lgt_S_c Free
      pM_lipidAds_c Free
      pM_lipidX_c Free
      pM_lpam_c Free
      pM_lycop_c Free
      pM_lys_L_c Free
      pM_lys_L_e Free
      pM_lys_L_p Free
      pM_lystrna_c Free
      pM_malACP_c Free
      pM_mal_L_c Free
      pM_mal_L_e Free
      pM_mal_L_p Free
      pM_malcoa_c Free
      pM_malt_c Free
      pM_man1p_c Free
      pM_man6p_c Free
      pM_mercppyr_c Free
      pM_met_L_c Free
      pM_methf_c Free
      pM_mettrna_c Free
      pM_mg2_c Free
      pM_mg2_e Free
      pM_mg2_p Free
      pM_mgdg160_c Free
      pM_mgdg161_c Free
      pM_mgdg180_c Free
      pM_mgdg181_9_c Free
      pM_mgdg181_c Free
      pM_mgdg182_9_12_c Free
      pM_mgdg183_6_9_12_c Free
      pM_mgdg183_9_12_15_c Free
      pM_mgdg184_6_9_12_15_c Free
      pM_mi1p_D_c Free
      pM_mlthf_c Free
      pM_mmtsa_c Free
      pM_mn2_c Free
      pM_mn2_e Free
      pM_mn2_p Free
      pM_mobd_c Free
      pM_mobd_e Free
      pM_mobd_p Free
      pM_mppp9_c Free
      pM_mppp9om_c Free
      pM_mthgxl_c Free
      pM_myrsACP_c Free
      pM_n4abutn_c Free
      pM_na1_c Free
      pM_na1_e Free
      pM_na1_p Free
      pM_nad_c Free
      pM_nadh_c Free
      pM_nadp_c Free
      pM_nadph_c Free
      pM_ncam_c Free
      pM_neuspn_c Free
      pM_nh4_c Free
      pM_nh4_e Free
      pM_nh4_p Free
      pM_ni2_c Free
      pM_ni2_e Free
      pM_ni2_p Free
      pM_nicrnt_c Free
      pM_nmn_c Free
      pM_no2_c Free
      pM_no3_c Free
      pM_no3_e Free
      pM_no3_p Free
      pM_o2_c Free
      pM_o2_e Free
      pM_o2_p Free
      pM_o2_u Free
      pM_oaa_c Free
      pM_ocACP_c Free
      pM_ocdcaACP_c Free
      pM_octeACP_c Free
      pM_octe_6_9_12_15_ACP_c Free
      pM_octe_6_9_12_ACP_c Free
      pM_octe_9_12_15_ACP_c Free
      pM_octe_9_12_ACP_c Free
      pM_octe_9_ACP_c Free
      pM_omppp9_c Free
      pM_orn_c Free
      pM_orot5p_c Free
      pM_orot_c Free
      pM_oxa_c Free
      pM_pa160_c Free
      pM_pa161_c Free
      pM_pa180_c Free
      pM_pa181_9_c Free
      pM_pa181_c Free
      pM_pa182_9_12_c Free
      pM_pa183_6_9_12_c Free
      pM_pa183_9_12_15_c Free
      pM_pa184_6_9_12_15_c Free
      pM_pac_c Free
      pM_pad_c Free
      pM_palmACP_c Free
      pM_pan4p_c Free
      pM_pant_R_c Free
      pM_pap_c Free
      pM_paps_c Free
      pM_pchlld_c Free
      pM_pcox_p Free
      pM_pcox_u Free
      pM_pcrd_p Free
      pM_pcrd_u Free
      pM_pdx5p_c Free
      pM_pep_c Free
      pM_peptido_syn_c Free
      pM_peptido_syn_p Free
      pM_pg160_c Free
      pM_pg161_c Free
      pM_pg180_c Free
      pM_pg181_9_c Free
      pM_pg181_c Free
      pM_pg182_9_12_c Free
      pM_pg183_6_9_12_c Free
      pM_pg183_9_12_15_c Free
      pM_pg184_6_9_12_15_c Free
      pM_pgp160_c Free
      pM_pgp161_c Free
      pM_pgp180_c Free
      pM_pgp181_9_c Free
      pM_pgp181_c Free
      pM_pgp182_9_12_c Free
      pM_pgp183_6_9_12_c Free
      pM_pgp183_9_12_15_c Free
      pM_pgp184_6_9_12_15_c Free
      pM_phbg_c Free
      pM_phdp_c Free
      pM_phe_L_c Free
      pM_pheme_c Free
      pM_phetrna_c Free
      pM_phllqne_c Free
      pM_phom_c Free
      pM_photon_c Free
      pM_photon_e Free
      pM_photon_p Free
      pM_phpyr_c Free
      pM_phthr_c Free
      pM_phytoe_c Free
      pM_phytof_c Free
      pM_pi_c Free
      pM_pi_e Free
      pM_pi_p Free
      pM_pmcoa_c Free
      pM_pnto_R_c Free
      pM_ppbng_c Free
      pM_pphn_c Free
      pM_ppi_c Free
      pM_ppp9_c Free
      pM_pppg9_c Free
      pM_pppi_c Free
      pM_pq_p Free
      pM_pq_u Free
      pM_pqh2_p Free
      pM_pqh2_u Free
      pM_pram_c Free
      pM_pran_c Free
      pM_prbamp_c Free
      pM_prbatp_c Free
      pM_pre3a_c Free
      pM_pre3b_c Free
      pM_pre4_c Free
      pM_pre5_c Free
      pM_pre6a_c Free
      pM_pre6b_c Free
      pM_pre8_c Free
      pM_precyanphy_c Free
      pM_prephytedp_c Free
      pM_prfp_c Free
      pM_prlp_c Free
      pM_pro_L_c Free
      pM_pro_L_e Free
      pM_pro_L_p Free
      pM_protrna_c Free
      pM_prpp_c Free
      pM_ptrc_c Free
      pM_ptrc_e Free
      pM_ptrc_p Free
      pM_pyam5p_c Free
      pM_pydam_c Free
      pM_pydx5p_c Free
      pM_pydx_c Free
      pM_pydxn_c Free
      pM_pyr_c Free
      pM_pyr_e Free
      pM_pyr_p Free
      pM_q8_c Free
      pM_q8h2_c Free
      pM_quln_c Free
      pM_r1p_c Free
      pM_r5p_c Free
      pM_rb15bp_c Free
      pM_rdmbzi_c Free
      pM_rib_D_c Free
      pM_ribflv_c Free
      pM_ru5p_D_c Free
      pM_s17bp_c Free
      pM_s7p_c Free
      pM_s_c Free
      pM_sbzcoa_c Free
      pM_ser_L_c Free
      pM_ser_L_e Free
      pM_ser_L_p Free
      pM_sertrna_c Free
      pM_skm5p_c Free
      pM_skm_c Free
      pM_sl26da_c Free
      pM_sl2a6o_c Free
      pM_so3_c Free
      pM_so4_c Free
      pM_so4_e Free
      pM_so4_p Free
      pM_spmd_c Free
      pM_spmd_e Free
      pM_spmd_p Free
      pM_sqdg160_c Free
      pM_sqdg161_c Free
      pM_sqdg180_c Free
      pM_sqdg181_9_c Free
      pM_sqdg181_c Free
      pM_sqdg182_9_12_c Free
      pM_sqdg183_6_9_12_c Free
      pM_sqdg183_9_12_15_c Free
      pM_sqdg184_6_9_12_15_c Free
      pM_sql_c Free
      pM_ssaltpp_c Free
      pM_suc6p_c Free
      pM_sucbz_c Free
      pM_succ_c Free
      pM_succ_e Free
      pM_succ_p Free
      pM_succoa_c Free
      pM_sucr_c Free
      pM_sucr_e Free
      pM_sucr_p Free
      pM_sucsal_c Free
      pM_t3c11vaceACP_c Free
      pM_t3c5ddeceACP_c Free
      pM_t3c7mrseACP_c Free
      pM_t3c9palmeACP_c Free
      pM_tag6p_D_c Free
      pM_tagdp_D_c Free
      pM_tddec2eACP_c Free
      pM_tdeACP_c Free
      pM_tdec2eACP_c Free
      pM_thdp_c Free
      pM_thex2eACP_c Free
      pM_thf_c Free
      pM_thfglu_c Free
      pM_thm_c Free
      pM_thmmp_c Free
      pM_thmpp_c Free
      pM_thr_L_c Free
      pM_thrp_c Free
      pM_thrtrna_c Free
      pM_thymd_c Free
      pM_tmrs2eACP_c Free
      pM_toct2eACP_c Free
      pM_toctd2eACP_c Free
      pM_tpalm2eACP_c Free
      pM_trdox_c Free
      pM_trdrd_c Free
      pM_trnaala_c Free
      pM_trnaarg_c Free
      pM_trnaasn_c Free
      pM_trnaasp_c Free
      pM_trnacys_c Free
      pM_trnaglu_c Free
      pM_trnagly_c Free
      pM_trnahis_c Free
      pM_trnaile_c Free
      pM_trnaleu_c Free
      pM_trnalys_c Free
      pM_trnamet_c Free
      pM_trnaphe_c Free
      pM_trnapro_c Free
      pM_trnaser_c Free
      pM_trnathr_c Free
      pM_trnatrp_c Free
      pM_trnatyr_c Free
      pM_trnaval_c Free
      pM_trp_L_c Free
      pM_trptrna_c Free
      pM_tyr_L_c Free
      pM_tyrtrna_c Free
      pM_u23ga_c Free
      pM_u3aga_c Free
      pM_u3hga_c Free
      pM_uGgl_c Free
      pM_uGgla_c Free
      pM_uaagmda_c Free
      pM_uaccg_c Free
      pM_uacgam_c Free
      pM_uagmda_c Free
      pM_uama_c Free
      pM_uamag_c Free
      pM_uamr_c Free
      pM_udcpdp_c Free
      pM_udcpp_c Free
      pM_udp_c Free
      pM_udpg_c Free
      pM_udpgal_c Free
      pM_udpsq_c Free
      pM_ugmd_c Free
      pM_ugmda_c Free
      pM_ump_c Free
      pM_uppg1_c Free
      pM_uppg3_c Free
      pM_ura_c Free
      pM_urea_c Free
      pM_urea_e Free
      pM_urea_p Free
      pM_utp_c Free
      pM_val_L_c Free
      pM_valtrna_c Free
      pM_xmp_c Free
      pM_xu5p_D_c Free
      pM_zcarote_c Free
      pM_zeax_c Free
      pM_zn2_c Free
      pM_zn2_e Free
      pM_zn2_p Free
End
