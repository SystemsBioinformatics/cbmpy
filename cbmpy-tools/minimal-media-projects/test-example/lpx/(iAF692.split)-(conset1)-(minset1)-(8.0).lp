\ENCODING=ISO-8859-1
\Problem name: iAF692

Minimize
 MNAF: bin_R_EXi_4abz_e + bin_R_EXi_4hphac_e + bin_R_EXi_ac_e
       + bin_R_EXi_actn__R_e + bin_R_EXi_ala__L_e + bin_R_EXi_alac__S_e
       + bin_R_EXi_btn_e + bin_R_EXi_ca2_e + bin_R_EXi_cbi_e + bin_R_EXi_cbl1_e
       + bin_R_EXi_cbl1hbi_e + bin_R_EXi_cd2_e + bin_R_EXi_ch4_e
       + bin_R_EXi_ch4s_e + bin_R_EXi_cit_e + bin_R_EXi_cl_e
       + bin_R_EXi_cobalt2_e + bin_R_EXi_cu2_e + bin_R_EXi_cys__L_e
       + bin_R_EXi_dma_e + bin_R_EXi_dms_e + bin_R_EXi_etha_e + bin_R_EXi_fe2_e
       + bin_R_EXi_fe3_e + bin_R_EXi_fol_e + bin_R_EXi_gcald_e
       + bin_R_EXi_glcn_e + bin_R_EXi_glu__L_e + bin_R_EXi_gly_e
       + bin_R_EXi_glyald_e + bin_R_EXi_glyb_e + bin_R_EXi_glyc_e
       + bin_R_EXi_h2s_e + bin_R_EXi_ile__L_e + bin_R_EXi_ind3ac_e
       + bin_R_EXi_k_e + bin_R_EXi_leu__L_e + bin_R_EXi_lys__L_e
       + bin_R_EXi_meoh_e + bin_R_EXi_mg2_e + bin_R_EXi_mma_e + bin_R_EXi_mn2_e
       + bin_R_EXi_mobd_e + bin_R_EXi_n2_e + bin_R_EXi_na1_e + bin_R_EXi_nac_e
       + bin_R_EXi_nh4_e + bin_R_EXi_ni2_e + bin_R_EXi_pac_e + bin_R_EXi_pi_e
       + bin_R_EXi_pnto__R_e + bin_R_EXi_pro__L_e + bin_R_EXi_pyr_e
       + bin_R_EXi_ribflv_e + bin_R_EXi_s_e + bin_R_EXi_so3_e + bin_R_EXi_so4_e
       + bin_R_EXi_urea_e + bin_R_EXi_val__L_e
Subject To
 M_10fthf_c:          - R_FMETTRS - R_AICART + R_MTHFC - R_GARFT  = 0
 M_13dpg_c:           R_GAPD_nadp + R_PGK + R_GAPD - R_ACYP  = 0
 M_1pyr5c_c:          R_G5SADs - R_P5CRx  = 0
 M_23dhdp_c:          - R_DHDPRy + R_DHDPS  = 0
 M_23dhmb_c:          - R_KARA1 - R_DHAD1  = 0
 M_23dhmp_c:          - R_DHAD2 + R_KARA2  = 0
 M_23dpg_c:           - R_CDPGS + R_CDPGH - R_DPGase + R_2PGK  = 0
 M_25aics_c:          R_PRASCSi - R_ADSL2r  = 0
 M_25dhpp_c:          - R_DROPPRy + R_TAPTPP - R_DROPPRx  = 0
 M_25dthpp_c:         R_DROPPRy - R_DRTPPD + R_DROPPRx  = 0
 M_26dap__M_c:        - R_DAPDC + R_DAPDH  = 0
 M_2ahbut_c:          - R_KARA2 + R_ACHBS  = 0
 M_2c25dho_c:         - R_CMLDC  = 0
 M_2cpr5p_c:          - R_IGPS + R_PRAIi  = 0
 M_2dhp_c:            R_MOHMT - R_DPR  = 0
 M_2dr1p_c:           R_DURIPP + R_TMDPP  = 0
 M_2ins_c:            R_INS2D  = 0
 M_2ippm_c:           R_IPPMIa - R_IPPMIb  = 0
 M_2mahmp_c:          - R_TMPPP + R_PMPK  = 0
 M_2obut_c:           R_METGL - R_POR3 - R_ACHBS  = 0
 M_2ood_c:            R_TIH3CD - R_AKSDC - R_AKSCAL  = 0
 M_2pg_c:             - R_ENO - R_PGM - R_2PGK  = 0
 M_2pglyc_c:          - R_PGLYCP  = 0
 M_2plac__L_c:        - R_LPPGS + R_2PLS  = 0
 M_2ppoh_c:           - R_ALCD20y  = 0
 M_2saa_c:            - R_SMTZS + R_SPDC  = 0
 M_2tcc_c:            R_SMTZS - R_SECS  = 0
 M_34hpp_c:           R_TYRTA - R_IOR3 + R_PPND - R_DMHDRFS  = 0
 M_36dahx_c:          - R_BLAT + R_LYSAM  = 0
 M_3c2hmp_c:          - R_IPMD - R_IPPMIa  = 0
 M_3c3hmp_c:          R_IPPMIb + R_IPPS  = 0
 M_3c4mop_c:          R_IPMD - R_OMCDC  = 0
 M_3dhq_c:            - R_DHQTi + R_DHQS2  = 0
 M_3dhsk_c:           R_DHQTi - R_SHK3Dr  = 0
 M_3hcdgggp_c:        - R_CDGGGPP4 - R_CDGGGSAT2 - R_CDGGIPT2 + R_CDGGGS2  = 0
 M_3hdgggp_c:         R_DGGGPS2 - R_CDGGGS2  = 0
 M_3hdgggps_c:        R_CDGGGSAT2 - R_HATSH  = 0
 M_3hdggpg_c:         R_DGGPGP2 - R_HATGH  = 0
 M_3hdggpgp_c:        - R_DGGPGP2 + R_CDGGGPP4  = 0
 M_3hdggpi_c:         R_CDGGIPT2 - R_HATIH  = 0
 M_3hdpgpe_c:         R_HASD - 0.002698 R_BIOMASS_Mb_30  = 0
 M_3hdpgpg_c:         R_HATGH - 0.01012 R_BIOMASS_Mb_30  = 0
 M_3hdpgpi_c:         - 0.01349 R_BIOMASS_Mb_30 + R_HATIH  = 0
 M_3hdpgps_c:         - R_HASD - 0.01147 R_BIOMASS_Mb_30 + R_HATSH  = 0
 M_3hfrdp_c:          R_GRTT2 - R_FRTT2  = 0
 M_3hggdp_c:          - R_DGGGPS2 + R_FRTT2  = 0
 M_3hgrdp_c:          - R_GRTT2 + R_DMATT2  = 0
 M_3ig3p_c:           R_IGPS - R_TRPS3 - R_TRPS1  = 0
 M_3mob_c:            R_VALTA - R_MOHMT - R_VOR + R_DHAD1 - R_IPPS  = 0
 M_3mop_c:            R_DHAD2 + R_ILETA  = 0
 M_3pg_c:             R_PGM - R_PGK - R_PGCD + R_DPGase + R_ACYP + 2 R_RBPC
                       = 0
 M_3php_c:            - R_PSERT + R_PGCD  = 0
 M_3psme_c:           R_PSCVT - R_CHORS  = 0
 M_3spyr_c:           - R_SLDxi2 + R_SLDx + R_SLDy - R_SPDC  = 0
 M_3uib_c:            R_DHPM2  = 0
 M_4abut_c:           R_GLUDC - R_ABTA  = 0
 M_4abz_c:            R_4ABZt2r - R_DHPS2 - R_BRFAPS  = 0
 M_4abz_e:            R_EXi_4abz_e - R_EXo_4abz_e - R_4ABZt2r  = 0
 M_4ahmmp_c:          R_AHMMPS - R_HMPK1  = 0
 M_4ampm_c:           R_HMPK1 - R_PMPK  = 0
 M_4hba_c:            R_THZPSN  = 0
 M_4hphac_c:          - R_PACCOAL2 + R_HPACt2r  = 0
 M_4hphac_e:          R_EXi_4hphac_e - R_EXo_4hphac_e - R_HPACt2r  = 0
 M_4kfbp_c:           - R_4KFBPP + R_FBPO  = 0
 M_4mhetz_c:          - R_HETZK  = 0
 M_4mop_c:            R_OMCDC + R_LEUTA  = 0
 M_4mpetz_c:          R_HETZK + R_THZPSN - R_TMPPP  = 0
 M_4pasp_c:           - R_ASADi + R_ASPK  = 0
 M_4ppan_c:           R_PNTK - R_PPNCL2  = 0
 M_4ppcys_c:          R_PPNCL2 - R_PPCDC  = 0
 M_4r5au_c:           R_RBFSb - R_RBFSa + R_PMDPHT - R_DMHDRFS  = 0
 M_56dthm_c:          - R_DHPM2  = 0
 M_56dura_c:          - R_DHPM1  = 0
 M_5aizc_c:           - R_PRASCSi - R_AIRC3  = 0
 M_5aop_c:            - 2 R_PPBNGS + R_G1SAT  = 0
 M_5aprbu_c:          R_DRTPPD - R_PMDPHT  = 0
 M_5caiz_c:           R_AIRC3 + R_AIRC2  = 0
 M_5dpmev_c:          R_PMEVK - R_DPMVD - R_DPMVD2  = 0
 M_5hbzid_c:          - R_NN5HBPRT + R_5HBZIDS  = 0
 M_5mdr1p_c:          R_MTAP  = 0
 M_5mta_c:            - R_MTAP  = 0
 M_5mthf_c:           R_MTHFR3 - R_METS  = 0
 M_5odhf2a_c:         R_CMLDC  = 0
 M_5oxpro_c:          - R_OPAH  = 0
 M_5pmev_c:           - R_PMEVK + R_MEVK1  = 0
 M_5pr5hbz_c:         R_NN5HBPRT - R_RZ5PP2  = 0
 M_6ax6ax_c:          - R_AHXDH  = 0
 M_6ax_c:             2 R_AHXDH  = 0
 M_6hmhpt_c:          R_DHNPA2r - R_HPPK2  = 0
 M_6hmhptpp_c:        - R_DHPS3 - R_DHPS2 + R_HPPK2  = 0
 M_6pgl_c:            R_G6PDH3  = 0
 M_6pthp_c:           R_PTHPS  = 0
 M_7mhp_c:            R_7MHS - R_7MHTS  = 0
 M_7mht_c:            - R_COBS + R_7MHTS  = 0
 M_7ohp_c:            - R_7MHS + R_AKSDC  = 0
 M_8aonn_c:           - R_AMAOTr + R_AOXSr  = 0
 M_Brfap_c:           - R_DHPS3 + R_BRFAPS  = 0
 M_Rh2cit_c:          R_AKACAL - R_RH2CD  = 0
 M_Rh3cit_c:          R_AKPCAL - R_RH3CD  = 0
 M_S2hglut_c:         - R_H4MPTS14 - R_AHGDx  = 0
 M_Shcit_c:           - R_SHCD + R_THACH  = 0
 M_aacoa_c:           - R_HMGCOASi + R_ACACT1r  = 0
 M_ac_c:              R_ALDD2x + R_AGAID + R_ACYP_2 + R_ACt2r + R_CYSS - R_ACKr
                      + R_AHSERL2  = 0
 M_ac_e:              R_EXi_ac_e - R_EXo_ac_e - R_ACt2r  = 0
 M_acald_c:           - R_ALDD2x + R_PYRDC + R_ALCD2y + R_THRAi + R_CPC5MT  = 0
 M_accoa_c:           - R_BLAT - R_PTAr - R_HMGCOASi - R_HCITS - R_AKPCAL
                      - R_AKACAL - 0.00005 R_BIOMASS_Mb_30 - R_AKGCAL - R_CODHr
                      - R_ACGS + R_POR2 - 2 R_ACACT1r - R_HSERTA - R_SERAT
                      - R_CS - R_IPPS - R_AKSCAL - R_G1PACT  = 0
 M_acetone_c:         R_ALCD20y  = 0
 M_acg5p_c:           R_AGPR + R_ACGK  = 0
 M_acg5sa_c:          - R_AGPR + R_ACOTA  = 0
 M_acgam1p_c:         R_G1PACT - R_UAGDP  = 0
 M_acglu_c:           - R_ACGK + R_ACGS + R_ORNTAC  = 0
 M_achms_c:           R_HSERTA - R_AHSERL2  = 0
 M_acmana_c:          R_UAG2EMA  = 0
 M_acon_T_c:          - R_ACONMT  = 0
 M_aconm_c:           R_ACONMT  = 0
 M_acorn_c:           - R_ACOTA - R_ORNTAC  = 0
 M_acser_c:           - R_CYSS + R_SERAT  = 0
 M_actn__R_c:         R_ACTNt2r + R_ACLDC  = 0
 M_actn__R_e:         - R_EXo_actn__R_e + R_EXi_actn__R_e - R_ACTNt2r  = 0
 M_actp_c:            R_PTAr - R_ACYP_2 + R_ACKr  = 0
 M_adcobdam_c:        - R_ADCYRS + R_CYRDAAT  = 0
 M_adcobhex_c:        R_ADOCBIAH + R_ADCYRS - R_ADCPS2  = 0
 M_ade_c:             - R_ADD - R_ADPT + R_MTAP  = 0
 M_adn_c:             R_AHCi - R_ADNK1  = 0
 M_adocbi_c:          - R_ADOCBIAH + R_CBIAT  = 0
 M_adocbip_c:         - R_ACBIPGT + R_ADCPS2  = 0
 M_adocblhbi_c:       R_CBLAT2 - 0.00474 R_BIOMASS_Mb_30 + R_ADOCBLS2  = 0
 M_adp_c:             R_CTPS2 + R_CTPS1 + R_PMEVK + R_TMK + R_PPK2r + R_CA2abc
                      + R_PRAGSr + R_ATPM + R_PRFGS + R_H4MPTS12 + R_H4MPTS13
                      + R_GLNTRAT + R_CBL1abc + R_GLU5K + R_MEVK1 + R_PNTK
                      + R_PRASCSi + R_DPMVD + R_HSK + R_ACCOAL + R_NDPK9
                      + R_NDPK8 + R_NDPK7 + R_NDPK6 + R_NDPK5 + R_NDPK4
                      + R_NDPK3 + R_NDPK2 + R_NDPK1 + R_UMPK_copy2 + R_ZNabc
                      + R_TSULabc + R_NADK + R_HETZK + R_DPCOAK + R_PGK + R_PC
                      - R_PFK_adp + R_SHKK + R_FE2abc + R_PROabc + R_HEX7
                      + 29.998 R_BIOMASS_Mb_30 + R_HKt + 2 R_F430S1 + R_TMDK1
                      + R_COabc - R_RNDR1 + 32 R_NIT_n1p4 - R_PYK + R_ACGK
                      + 2 R_CYRDAS + 4 R_ADCYRS + R_GLNS + R_LTHRK + R_DTMPK
                      + R_CD2abc2 + R_CD2abc1 + R_NTP1 + R_ATPS1 + R_CYTK1
                      + R_CYTK2 + R_TMPK + R_CBL1HBIabc + R_AIRC2 + R_CDPGS
                      + R_DBTS + R_PIabc + R_GK1 + R_RBK + R_Cuabc + R_APPLP
                      + R_HMPK1 + R_DHFS + R_FEDCabc + R_NaKt_1 + R_UMPK_copy1
                      + R_URIDK2r + R_NIabc - R_ATPS4r + R_SULabc + R_CBIabc
                      + R_COCHL + R_PMPK - R_PPAKr + R_MG2abc + R_FE3abc
                      + R_COBS + 2 R_CBPS + 2 R_ADK1 + R_ADK3 + R_ADK2 + R_ADK4
                      + R_MNabc + R_ADCPS2 + R_OPAH + R_7MHTS + R_2PGK
                      + R_MOBDabc + R_ACKr + R_GLYBabc + R_Cut1 + R_ADNK1
                      + R_THMabc + R_H4MPTGL_atp - R_GK_adp + R_DADK + R_MPML
                      + R_PRAIS + R_ASPK  = 0
 M_agdpcbi_c:         R_ACBIPGT - R_ADOCBLS2  = 0
 M_agdpgpi_c:         R_AGAIAGT - R_AGAID  = 0
 M_agm_c:             R_ARGDC - R_AGMT  = 0
 M_ahcys_c:           R_CPC3MT + R_H4MPTS16 - R_AHCi + R_ACONMT + 2 R_UPP3MT
                      + R_H4MPTS15 + 2 R_CPC6MT + R_CPC2MT + R_CPC4MT
                      + R_CPC5MT + R_DIPS  = 0
 M_ahdt_c:            - R_PTHPS - R_AKP1  = 0
 M_aicar_c:           - R_AICART + R_IG3PS + R_ADSL2r  = 0
 M_air_c:             - R_AHMMPS - R_AIRC2 + R_PRAIS  = 0
 M_akg_c:             - R_GLUS_F420 + R_HSTPT - R_ASPTA - R_TYRTA - R_OORr
                      - R_ABTA - R_VALTA - R_PHETA1 - R_GLUDyi - R_HCITS
                      - R_TRPTA + R_PSERT - R_AKGCAL - R_ACOTA + R_GLUDxi
                      + R_ICDHyr - R_LEUTA + R_AHGDx - R_ALATA_L - R_ILETA  = 0
 M_ala_B_c:           - R_PANTS + R_ASP1DC + R_NBAHH_ir  = 0
 M_ala__L_c:          - 0.5621 R_BIOMASS_Mb_30 - R_ALATRS + R_THZPSN + R_ALAt4
                      - R_ALATA_L - R_AOXSr  = 0
 M_ala__L_e:          R_EXi_ala__L_e - R_EXo_ala__L_e - R_ALAt4  = 0
 M_alac__S_c:         R_KARA1 + R_ACLS + R_ALACt2r - R_ACLDC  = 0
 M_alac__S_e:         - R_EXo_alac__S_e + R_EXi_alac__S_e - R_ALACt2r  = 0
 M_alatrna_c:         R_ALATRS  = 0
 M_amet_c:            - R_CPC3MT - R_H4MPTS16 - R_AMAOTr + R_METAT - R_ACONMT
                      - 2 R_UPP3MT - R_H4MPTS15 - 2 R_CPC6MT - R_CPC2MT
                      - R_CPC4MT - R_CPC5MT - R_DIPS  = 0
 M_amob_c:            R_AMAOTr  = 0
 M_amp_c:             R_PRPPS + R_PANTS + R_CF3Ha + R_SERTRS + R_ACS2
                      + R_GLUTRS + R_PACCOAL + R_PROTRS + R_ADSL1r + R_GMPS
                      + R_ARGSS + R_PFK_adp + R_GMPS2 + R_ASPTRS + R_TRPTRS
                      + R_NADS1 + R_NADS2 - 0.001 R_BIOMASS_Mb_30 + R_PACCOAL2
                      + R_PACCOAL3 + R_PPS + R_ALATRS + R_ARGTRS + R_THZPSN
                      + R_ILETRS + R_LEUTRS + R_GLYTRS + R_GLUTRS_Gln
                      + R_VALTRS + R_METTRS - R_ADK1 - R_ADK3 - R_ADK2 - R_ADK4
                      + R_ASNS1 + R_CYSTRS + R_ADPT + R_ADNK1 + R_TYRTRS
                      + R_LYSTRS + R_GK_adp + R_HPPK2 + R_HISTRS + R_THRTRS
                      + R_PPDK + R_PHETRS  = 0
 M_anth_c:            R_ANS - R_ANPRT  = 0
 M_appl_c:            R_ADOCBIAH - R_APPLP  = 0
 M_applp_c:           R_APPLP + R_THRPD - R_ADCPS2  = 0
 M_arg__L_c:          - R_ARGDC - R_ARGDr - 0.3237 R_BIOMASS_Mb_30 - R_ARGTRS
                      + R_ARGSL  = 0
 M_argsuc_c:          R_ARGSS - R_ARGSL  = 0
 M_argtrna_c:         R_ARGTRS  = 0
 M_asn__L_c:          - R_ASNN - 0.2638 R_BIOMASS_Mb_30 + R_ASNS1  = 0
 M_asp__L_c:          R_ASNN - R_ASP1DC - R_ASPCT - R_PRASCSi - R_ASPTA
                      - R_ARGSS - R_ASPTRS - 0.2638 R_BIOMASS_Mb_30 - R_ADSS
                      - R_ASNS1 - R_ASPK  = 0
 M_aspsa_c:           R_ASADi - R_DKMPASPL - R_DHDPS + R_HSDy  = 0
 M_asptrna_c:         R_ASPTRS  = 0
 M_atp_c:             - R_CTPS2 - R_CTPS1 - R_PRPPS - R_PANTS - R_PMEVK - R_TMK
                      - R_SERTRS - R_NNATr - R_PPK2r - R_ACS2 - R_CA2abc
                      - R_PRAGSr - R_ATPM - R_PRFGS - R_H4MPTS12 - R_H4MPTS13
                      - R_GLNTRAT - R_GLUTRS - R_CBL1abc - R_GLU5K - R_MEVK1
                      - R_PNTK - R_PRASCSi - R_PACCOAL - R_DPMVD - R_ADNCYC
                      - R_PROTRS - R_HSK - R_ACCOAL - R_NDPK9 - R_NDPK8
                      - R_NDPK7 - R_NDPK6 - R_NDPK5 - R_NDPK4 - R_NDPK3
                      - R_NDPK2 - R_NDPK1 - R_UMPK_copy2 - R_ZNabc - R_RNTR1
                      - R_NMNAT - R_GMPS - R_TSULabc - R_NADK - R_ARGSS
                      - R_HETZK - R_DPCOAK - R_PGK - R_PC - R_CBLAT2 - R_GMPS2
                      - R_SHKK - R_FE2abc - R_ASPTRS - R_TRPTRS - R_PROabc
                      - R_NADS1 - R_NADS2 - R_HEX7 - 30.1846 R_BIOMASS_Mb_30
                      - R_PACCOAL2 - R_PACCOAL3 - R_HKt - 2 R_F430S1 - R_TMDK1
                      - R_METAT - R_COabc - 32 R_NIT_n1p4 + R_PYK - R_ACGK
                      - R_PPS - 2 R_CYRDAS - 4 R_ADCYRS - R_GLNS - R_LTHRK
                      - R_ALATRS - R_BACCL - R_DTMPK - R_CD2abc2 - R_CD2abc1
                      - R_NTP1 - R_ARGTRS - R_ATPS1 - R_THZPSN - R_CYTK1
                      - R_CYTK2 - R_TMPK - R_CBL1HBIabc - R_AIRC2 - R_CDPGS
                      - R_DBTS - R_ILETRS - R_PIabc - R_LEUTRS - R_GK1
                      - R_GLYTRS - R_GLUTRS_Gln - R_RBK - R_Cuabc - R_APPLP
                      - R_HMPK1 - R_DHFS - R_FEDCabc - R_NaKt_1 - R_VALTRS
                      - R_UMPK_copy1 - R_URIDK2r - R_CYRDAAT - R_METTRS
                      - R_NIabc + R_ATPS4r - R_SULabc - R_CBIabc - R_COCHL
                      - R_PMPK + R_PPAKr - R_MG2abc - R_FE3abc - R_COBS
                      - 2 R_CBPS - R_ADK1 - R_MNabc - R_ADCPS2 - R_ASNS1
                      - R_OPAH - R_7MHTS - R_2PGK - R_MOBDabc - R_ACKr
                      - R_CYSTRS - R_GLYBabc - R_Cut1 - R_ADNK1 - R_THMabc
                      - R_TYRTRS - R_ATPPRT - R_LYSTRS - R_H4MPTGL_atp
                      - R_CBIAT - R_DADK - R_HPPK2 - R_MPML - R_PTPATi
                      - R_HISTRS - R_THRTRS - R_PRAIS - R_CF3Sa - R_PPDK
                      - R_PHETRS - R_ASPK  = 0
 M_atrz_c:            - R_ATRZCH  = 0
 M_btamp_c:           R_BACCL  = 0
 M_btn_c:             - R_BACCL + R_BTNt2i + R_BTS3r  = 0
 M_btn_e:             - R_EXo_btn_e + R_EXi_btn_e - R_BTNt2i  = 0
 M_ca2_c:             R_CA2abc - R_CAt6  = 0
 M_ca2_e:             - R_CA2abc + R_CAt6 + R_EXi_ca2_e - R_EXo_ca2_e  = 0
 M_cala_c:            R_DHPM1  = 0
 M_camp_c:            R_ADNCYC  = 0
 M_caphis_c:          - R_DIPS  = 0
 M_carn_c:            - R_NBAHH_ir  = 0
 M_cbasp_c:           R_ASPCT + R_DHORTS  = 0
 M_cbi_c:             R_CBIabc - R_CBIAT + R_CBLD  = 0
 M_cbi_e:             - R_CBIabc + R_EXi_cbi_e - R_EXo_cbi_e  = 0
 M_cbl1_c:            R_CBL1abc - R_CBLD  = 0
 M_cbl1_e:            - R_CBL1abc - R_EXo_cbl1_e + R_EXi_cbl1_e  = 0
 M_cbl1hbi_c:         - R_CBLAT2 + R_CBL1HBIabc  = 0
 M_cbl1hbi_e:         - R_EXo_cbl1hbi_e - R_CBL1HBIabc + R_EXi_cbl1hbi_e  = 0
 M_cbp_c:             - R_ASPCT - R_OCBT + R_CBPS  = 0
 M_cd2_c:             R_CD2abc2 - R_CD2abc1  = 0
 M_cd2_e:             - R_EXo_cd2_e + R_EXi_cd2_e - R_CD2abc2 + R_CD2abc1  = 0
 M_cdgggp_c:          - R_CDGGIPT - R_CDGGGPP3 + R_CDGGGS - R_CDGGGSAT  = 0
 M_cdp_c:             - R_NDPK3 - R_RNDR3 + R_NTP5 + R_CYTK1  = 0
 M_cdpg_c:            R_CDPGS - R_CDPGH  = 0
 M_ch4_c:             R_CH4t + R_MCR  = 0
 M_ch4_e:             - R_CH4t + R_EXi_ch4_e - R_EXo_ch4_e  = 0
 M_ch4s_c:            R_METGL + R_CH4St + R_MTCMMT  = 0
 M_ch4s_e:            R_EXi_ch4s_e - R_CH4St - R_EXo_ch4s_e  = 0
 M_chor_c:            - R_ANS - R_CHORM + R_CHORS  = 0
 M_cit_c:             - R_ACONT + 2 R_FEDCabc + R_CS  = 0
 M_cit_e:             - 2 R_FEDCabc + R_EXi_cit_e - R_EXo_cit_e  = 0
 M_citr__L_c:         R_OCBT - R_ARGSS + R_ARGDr  = 0
 M_cl_c:              R_ATRZCH + R_Clt  = 0
 M_cl_e:              - R_EXo_cl_e + R_EXi_cl_e - R_Clt  = 0
 M_cmaphis_c:         R_DIPS  = 0
 M_cmp_c:             R_CDGGIPT + R_CDGGGPP4 + R_CDGGGPP3 + R_PPNCL2 - R_CYTK1
                      + R_CDGGGSAT2 + R_CDGGIPT2 + R_CDGGGSAT  = 0
 M_co1dam_c:          2 R_CYRDAR - R_CYRDAAT  = 0
 M_co2_c:             R_POR3 + R_GLUDC + R_ASP1DC + R_PYRDC + R_DPMVD + R_ORNDC
                      + R_TIH2CD + R_HASD + R_CO2t + R_OORr + R_FDH + R_TIH3CD
                      + R_IGPS + R_ARGDC + R_ASD + R_HISDC + R_BRFAPS + R_CMLDC
                      + R_DPMVD2 - R_FMFD_b + R_PPNDH + R_CODHr - R_HCO3E
                      + R_IOR + R_NNDPR + R_PPCDC + R_IOR3 + R_IOR2 + R_POR2
                      + R_ACLS + R_THZPSN + R_DAPDC + R_CODH2r + R_VOR - R_DBTS
                      + R_THRPD + R_ICDHyr + R_AKSDC + R_PPND + R_OMCDC + R_ME1
                      + R_CPC6MT + R_OMPDC + R_TIHCD + R_OXADC - R_RBPC
                      + R_AOXSr + R_ACHBS + R_SPDC + R_ACLDC  = 0
 M_co2_e:             - R_CO2t - R_EXo_co2_e + R_EXi_co2_e  = 0
 M_co2dam_c:          - 2 R_CYRDAR + R_CYRDAS + R_COCHL  = 0
 M_co_c:              - R_CODH2r - R_COt  = 0
 M_co_e:              R_EXi_co_e - R_EXo_co_e + R_COt  = 0
 M_coa_c:             R_HMGCOARi + R_BLAT - R_ACS2 - R_POR3 - R_PACCOAL
                      - R_ACCOAL + R_PTAr - R_OORr + R_DPCOAK + R_HMGCOASi
                      + R_HCITS + R_AKPCAL + R_AKACAL
                      - 0.000006 R_BIOMASS_Mb_30 - R_PACCOAL2 - R_PACCOAL3
                      + R_AKGCAL + R_CODHr + R_ACGS - R_IOR - R_IOR3 - R_IOR2
                      - R_POR2 - R_VOR + R_ACACT1r + R_HSERTA + R_SERAT + R_CS
                      + R_IPPS + R_AOXSr + R_G1PACT  = 0
 M_cob_c:             - R_MCR - 0.0005 R_BIOMASS_Mb_30 + R_HDR + R_COBS  = 0
 M_cobalt2_c:         R_COabc - R_COCHL - R_COBALTt5 - R_SHCHCC  = 0
 M_cobalt2_e:         - R_COabc + R_EXi_cobalt2_e - R_EXo_cobalt2_e
                      + R_COBALTt5  = 0
 M_cobya_c:           - R_CYRDAS + R_CPC8MM  = 0
 M_codhpre6_c:        R_CPC6R - R_CPC6MT  = 0
 M_com_c:             R_COMS - R_DMAMT - 0.02056 R_BIOMASS_Mb_30 - R_MTCMMT
                      + R_HDR - R_MCMMT - R_TMAMT - R_MMAMT - R_MTSPCMMT  = 0
 M_copre2_c:          - R_CPC2MT + R_SHCHCC  = 0
 M_copre3_c:          - R_CPC3MT + R_CPC2MT  = 0
 M_copre4_c:          R_CPC3MT - R_CPC4MT  = 0
 M_copre5_c:          R_CPC4MT - R_CPC5MT  = 0
 M_copre6_c:          - R_CPC6R + R_CPC5MT  = 0
 M_copre8_c:          R_CPC6MT - R_CPC8MM  = 0
 M_csn_c:             - R_CSND  = 0
 M_ctp_c:             R_CTPS2 + R_CTPS1 + R_NDPK3 - R_RNTR3
                      - 0.1379 R_BIOMASS_Mb_30 - R_PPNCL2 - R_NTP5 - R_CDGGGS
                      - R_CDGGGS2  = 0
 M_cu2_c:             R_Cuabc - R_Cut1  = 0
 M_cu2_e:             R_EXi_cu2_e - R_EXo_cu2_e - R_Cuabc + R_Cut1  = 0
 M_cys__L_c:          - R_7MHS + R_CYSt2r - R_SMTZS - 0.1002 R_BIOMASS_Mb_30
                      - R_PPNCL2 - R_CYSDS - R_THZPSN + R_CYSS - R_CYSTRS  = 0
 M_cys__L_e:          - R_CYSt2r - R_EXo_cys__L_e + R_EXi_cys__L_e  = 0
 M_cystrna_c:         R_CYSTRS  = 0
 M_dadp_c:            - R_NDPK8 + R_RNDR1 + R_NTP2 + 2 R_ADKd + R_DADK  = 0
 M_damp_c:            - R_ADKd - R_DADK  = 0
 M_dann_c:            R_AMAOTr - R_DBTS  = 0
 M_datp_c:            R_NDPK8 + R_RNTR1 - 0.03312 R_BIOMASS_Mb_30 - R_NTP2
                      - R_ADKd  = 0
 M_db4p_c:            R_DB4PS - R_RBFSa  = 0
 M_dcamp_c:           - R_ADSL1r + R_ADSS  = 0
 M_dcdp_c:            - R_NDPK7 + R_RNDR3 + R_NTP6 + R_CYTK2  = 0
 M_dcmp_c:            - R_CYTK2 - R_DCMPDA  = 0
 M_dctp_c:            R_NDPK7 + R_RNTR3 - 0.02153 R_BIOMASS_Mb_30 - R_NTP6
                      - R_DCTPD  = 0
 M_ddhrb_c:           - R_LPFPLT + R_DMHDRFS  = 0
 M_dgdp_c:            - R_NDPK5 + R_RNDR2 + R_NTP4  = 0
 M_dgggp_c:           - R_CDGGGS + R_DGGGPS  = 0
 M_dgggps_c:          - R_ATSH + R_CDGGGSAT  = 0
 M_dggpg_c:           - R_ATGH + R_DGGPGP  = 0
 M_dggpgp_c:          - R_DGGPGP + R_CDGGGPP3  = 0
 M_dggpi_c:           R_CDGGIPT - R_ATIH  = 0
 M_dgtp_c:            R_NDPK5 + R_RNTR2 - 0.02153 R_BIOMASS_Mb_30 - R_NTP4  = 0
 M_dha_c:             R_ALKP  = 0
 M_dhadr_c:           R_H4MPTS10 - R_H4MPTS11  = 0
 M_dhadrdpr_c:        R_H4MPTS12 - R_H4MPTS13  = 0
 M_dhadrp_c:          R_H4MPTS9 - R_H4MPTS10  = 0
 M_dhadrpr_c:         R_H4MPTS11 - R_H4MPTS12  = 0
 M_dhadrtpr_c:        R_H4MPTS13 - R_H4MPTS14  = 0
 M_dhap_c:            R_FCLPA - R_G1PDH + R_FBA2 - R_TPI - R_ALKP - R_QULNS
                      + R_FBA  = 0
 M_dhf_c:             R_TMDS + R_FOLR2 + R_DHFS - R_DHFR  = 0
 M_dhnpt_c:           - R_DHNPA2r + R_DNMPPA + R_AKP1  = 0
 M_dhor__S_c:         - R_DHORTS - R_DHORD7  = 0
 M_dhp23cp_c:         - R_DHPCPH + R_TAPTPPM  = 0
 M_dhpmp_c:           - R_DNMPPA + R_DHPCPH  = 0
 M_dhpt_c:            R_DHPS2 - R_DHFS  = 0
 M_dhrfap_c:          - R_H4MPTS9 + R_DHPS3  = 0
 M_dkdofp_c:          - R_DKFPS + R_4KFBPP  = 0
 M_dkfp_c:            R_DKFPS - R_DKFPR  = 0
 M_dkmp_c:            R_DKFPR - R_DKMPASPL  = 0
 M_dma_c:             - R_DMAMT + R_TMAMT + R_DMAt  = 0
 M_dma_e:             - R_EXo_dma_e + R_EXi_dma_e - R_DMAt  = 0
 M_dmh2mpt_c:         R_H4MPTS14 - R_H4MPTS15  = 0
 M_dmlz_c:            - 2 R_RBFSb + R_RBFSa  = 0
 M_dmpp_c:            R_IPDDI3y + R_IPDDI3x - R_DMATT  = 0
 M_dms_c:             R_DMSt - R_MTCMMT  = 0
 M_dms_e:             - R_DMSt - R_EXo_dms_e + R_EXi_dms_e  = 0
 M_dnad_c:            R_NNATr - R_NADS1 - R_NADS2  = 0
 M_dohau_c:           - R_DOHDUS + R_DKMPASPL  = 0
 M_dohdu_c:           R_DOHDUS - R_DHQS2  = 0
 M_dpcoa_c:           - R_DPCOAK + R_PTPATi  = 0
 M_dpgpe_c:           R_ASD - 0.000674 R_BIOMASS_Mb_30  = 0
 M_dpgpg_c:           R_ATGH - 0.000674 R_BIOMASS_Mb_30  = 0
 M_dpgpi_c:           - R_AGAIAGT - 0.002698 R_BIOMASS_Mb_30 + R_ATIH  = 0
 M_dpgps_c:           - R_ASD - 0.001349 R_BIOMASS_Mb_30 + R_ATSH  = 0
 M_dscl_c:            - R_SHCHD2 - R_F430S1 + R_UPP3MT - R_SHCHCC  = 0
 M_dtbt_c:            R_DBTS - R_BTS3r  = 0
 M_dtdp4d6dg_c:       - R_TDPDRE + R_TDPGDH  = 0
 M_dtdp4d6dm_c:       - R_TDPDRR + R_TDPDRE  = 0
 M_dtdp_c:            - R_NDPK4 + R_DTMPK + R_NTP9  = 0
 M_dtdpglu_c:         R_G1PTT - R_TDPGDH  = 0
 M_dtdprmn_c:         R_TDPDRR  = 0
 M_dtmp_c:            R_TMDS + R_TMDK1 - R_DTMPK  = 0
 M_dttp_c:            R_NDPK4 - 0.03312 R_BIOMASS_Mb_30 - R_NTP9 - R_G1PTT  = 0
 M_dudp_c:            - R_NDPK6 + R_RNDR4 + R_NTP8 + R_URIDK2r  = 0
 M_dump_c:            - R_TMDS - R_URIDK2r + R_DCMPDA  = 0
 M_duri_c:            - R_DURIPP  = 0
 M_dutp_c:            R_NDPK6 + R_RNTR4 - R_NTP8 + R_DCTPD  = 0
 M_dxyl5p_c:          - R_THZPSN  = 0
 M_e4p_c:             R_TALA - R_TKT2  = 0
 M_eig3p_c:           - R_IGPDH + R_IG3PS  = 0
 M_etha_c:            R_ETHAt6  = 0
 M_etha_e:            - R_ETHAt6 - R_EXo_etha_e + R_EXi_etha_e  = 0
 M_etoh_c:            - R_ALCD2y  = 0
 M_f1p_c:             - R_FBA2  = 0
 M_f390a_c:           - R_CF3Ha - 0.000005 R_BIOMASS_Mb_30 + R_CF3Sa  = 0
 M_f390g_c:           - R_CF3Hg - 0.000005 R_BIOMASS_Mb_30 + R_CF3Sg  = 0
 M_f420_0_c:          R_LPFPLT - R_GF4GL_0  = 0
 M_f420_1_c:          - R_GF4GL_1 + R_GF4GL_0  = 0
 M_f420_2_c:          - R_G6PDH3 + R_CF3Hg + R_CF3Ha - R_F4RHi + R_H4MPTS9
                      + R_F4MTSPD + R_F4MTSPR + R_GLUS_F420 - R_GF4GL_2
                      + R_GF4GL_1 - R_F4NPR - 0.000043 R_BIOMASS_Mb_30
                      - R_DHORD7 + R_F4D - R_IPDDI3y - R_IPDDI3x + R_H2MPTR
                      + 2 R_F4H2O - R_CF3Sg - R_CF3Sa  = 0
 M_f420_2h2_c:        R_G6PDH3 + R_F4RHi - R_H4MPTS9 - R_F4MTSPD - R_F4MTSPR
                      - R_GLUS_F420 + R_F4NPR + R_DHORD7 - R_F4D + R_IPDDI3y
                      + R_IPDDI3x - R_H2MPTR - 2 R_F4H2O  = 0
 M_f420_3_c:          - R_GF4GL_3 + R_GF4GL_2 - 0.000026 R_BIOMASS_Mb_30  = 0
 M_f420_4_c:          - R_GF4GL_4 + R_GF4GL_3 - 0.000374 R_BIOMASS_Mb_30  = 0
 M_f420_5_c:          - R_GF4GL_5 + R_GF4GL_4 - 0.000299 R_BIOMASS_Mb_30  = 0
 M_f420_6_c:          - R_GF4GL_6 + R_GF4GL_5 - 0.000011 R_BIOMASS_Mb_30  = 0
 M_f420_7_c:          R_GF4GL_6 - 0.000001 R_BIOMASS_Mb_30  = 0
 M_f430_c:            - 0.002 R_BIOMASS_Mb_30 + R_F430S5  = 0
 M_f430p1_c:          - R_F430S2 + R_F430S1  = 0
 M_f430p2_c:          - R_F430S3 + R_F430S2  = 0
 M_f430p3_c:          R_F430S3 - R_F430S4  = 0
 M_f6p_c:             - R_GF6PTA + R_PGI - R_PFK_adp + R_HEX7 + R_TALA + R_TKT2
                      + R_FBP + R_MAN6PI  = 0
 M_fald_c:            - R_ALDD1  = 0
 M_fc1p_c:            - R_FCLPA  = 0
 M_fdox_c:            - 2 R_POR3 - 2 R_OORr + 2 R_FMFD_b - 2 R_CODHr
                      + 8 R_NIT_n1p4 - 2 R_IOR - 2 R_IOR3 - 2 R_IOR2 - 2 R_POR2
                      - 2 R_CODH2r - 2 R_VOR + 2 R_ECHH_10  = 0
 M_fdp_c:             R_PFK_adp - R_FBPO - R_FBP - R_FBA  = 0
 M_fdred_c:           2 R_POR3 + 2 R_OORr - 2 R_FMFD_b + 2 R_CODHr
                      - 8 R_NIT_n1p4 + 2 R_IOR + 2 R_IOR3 + 2 R_IOR2 + 2 R_POR2
                      + 2 R_CODH2r + 2 R_VOR - 2 R_ECHH_10  = 0
 M_fe2_c:             R_FE2abc - R_SHCHF  = 0
 M_fe2_e:             - R_FE2abc - R_EXo_fe2_e + R_EXi_fe2_e  = 0
 M_fe3_c:             R_FEDCabc + R_FE3abc  = 0
 M_fe3_e:             - R_FEDCabc + R_EXi_fe3_e - R_FE3abc - R_EXo_fe3_e  = 0
 M_fgam_c:            - R_PRFGS + R_GARFT  = 0
 M_fmettrna_c:        R_FMETTRS  = 0
 M_fmn_c:             - R_ACP1_FMN  = 0
 M_fol_c:             R_FOLt - R_FOLR2  = 0
 M_fol_e:             - R_FOLt - R_EXo_fol_e + R_EXi_fol_e  = 0
 M_for_c:             - R_FDH + R_DB4PS + R_ALDD1 + R_OXADC + R_GTPCIII  = 0
 M_formh4spt_c:       R_FMFTSPFT_b - R_MTSPC  = 0
 M_formmfr_b_c:       - R_FMFTSPFT_b + R_FMFD_b  = 0
 M_fpram_c:           R_PRFGS - R_PRAIS  = 0
 M_fprica_c:          R_IMPC + R_AICART  = 0
 M_frdp_c:            R_GRTT - R_FRTT  = 0
 M_fru_c:             - R_HEX7  = 0
 M_fum_c:             R_ADSL1r - R_FUM + R_ADSL2r + R_ARGSL  = 0
 M_g1p_c:             R_GLCP - R_GALUi - R_PGMT - R_G1PTT  = 0
 M_g3p_c:             - R_GAPD_nadp + R_TPI + R_TRPS3 + R_TRPS1 - R_TALA
                      - R_GAPD + R_DKMPASPL + R_TKT1 + R_TKT2 + R_FBA  = 0
 M_g6p_c:             - R_G6PDH3 - R_PGI + R_PGMT - R_MI1PS + R_GK_adp  = 0
 M_gal1p_c:           - R_GALT  = 0
 M_gam1p_c:           - R_PGAMT - R_G1PACT  = 0
 M_gam6p_c:           R_GF6PTA + R_PGAMT  = 0
 M_gar_c:             R_PRAGSr - R_GARFT  = 0
 M_gcald_c:           R_DHNPA2r + R_AHMMPS + R_GCALDt  = 0
 M_gcald_e:           R_EXi_gcald_e - R_EXo_gcald_e - R_GCALDt  = 0
 M_gdp_c:             R_H4MPTGL_gtp - R_NDPK1 + R_GF4GL_6 + R_GF4GL_5
                      + R_GF4GL_4 + R_GF4GL_3 + R_GF4GL_2 + R_GF4GL_1
                      + R_GF4GL_0 - R_RNDR2 + R_ADSS + R_NTP3 + R_GK1
                      - R_MAN1PT2 + R_ADK3 + R_2PLS  = 0
 M_gdpddman_c:        - R_GDMANE + R_GMAND  = 0
 M_gdpfuc_c:          R_GOFUCR  = 0
 M_gdpgpi_c:          - 0.01619 R_BIOMASS_Mb_30 + R_AGAID  = 0
 M_gdpmann_c:         - R_GMAND + R_MAN1PT2  = 0
 M_gdpofuc_c:         R_GDMANE - R_GOFUCR  = 0
 M_ggdp_c:            - R_GGGPS + R_FRTT - R_DGGGPS  = 0
 M_gggp_c:            R_GGGPS - R_DGGGPS2 - R_DGGGPS  = 0
 M_glc__D_c:          R_GLCGSD - R_GK_adp  = 0
 M_glcn_c:            R_GLCNt2r  = 0
 M_glcn_e:            R_EXi_glcn_e - R_GLCNt2r - R_EXo_glcn_e  = 0
 M_gln__L_c:          - R_CTPS2 - R_PRFGS - R_GLNTRAT - R_GLUS_F420 - R_ANS
                      - R_GF6PTA - R_GMPS2 - R_NADS2 - 0.288 R_BIOMASS_Mb_30
                      - 2 R_F430S1 - 2 R_CYRDAS - 4 R_ADCYRS + R_GLNS - R_CBPS
                      - R_ASNS1 - R_IG3PS - R_GLUPRT  = 0
 M_glntrna_c:         R_GLNTRAT  = 0
 M_glu1sa_c:          R_GLUTRR - R_G1SAT  = 0
 M_glu5p_c:           - R_G5SD2 + R_GLU5K  = 0
 M_glu5sa_c:          R_G5SD2 - R_G5SADs  = 0
 M_glu__L_c:          R_CTPS2 - R_H4MPTGL_gtp + R_PRFGS - R_GLUDC + R_GLNTRAT
                      + R_GLUt2r + 2 R_GLUS_F420 - R_GLUTRS + R_ANS - R_HSTPT
                      - R_GLU5K + R_GF6PTA + R_ASPTA - R_GF4GL_6 - R_GF4GL_5
                      - R_GF4GL_4 - R_GF4GL_3 - R_GF4GL_2 - R_GF4GL_1
                      - R_GF4GL_0 + R_TYRTA + R_ABTA + R_VALTA + R_PHETA1
                      + R_GLUDyi + R_GMPS2 + R_TRPTA + R_NADS2 - R_PSERT
                      - 0.288 R_BIOMASS_Mb_30 + 2 R_F430S1 - R_ACGS
                      + 2 R_CYRDAS + 4 R_ADCYRS - R_GLNS + R_ACOTA - R_GLUDxi
                      - R_GLUTRS_Gln - R_DHFS + R_LEUTA + R_CBPS + R_ASNS1
                      + R_OPAH + R_ALATA_L + R_IG3PS - R_H4MPTGL_atp + R_GLUPRT
                      + R_ILETA - R_ORNTAC  = 0
 M_glu__L_e:          - R_GLUt2r - R_EXo_glu__L_e + R_EXi_glu__L_e  = 0
 M_glutrna_c:         R_GLUTRS - R_GLUTRR  = 0
 M_glutrna_gln_c:     - R_GLNTRAT + R_GLUTRS_Gln  = 0
 M_gly_c:             - R_PRAGSr - 0.6704 R_BIOMASS_Mb_30 + R_GLYt4 - R_GLYTRS
                      + R_GHMT2r + R_THRAi  = 0
 M_gly_e:             R_EXi_gly_e - R_EXo_gly_e - R_GLYt4  = 0
 M_glyald_c:          R_FBA2 + R_GLYALDt  = 0
 M_glyald_e:          - R_EXo_glyald_e - R_GLYALDt + R_EXi_glyald_e  = 0
 M_glyb_c:            R_GLYBabc  = 0
 M_glyb_e:            R_EXi_glyb_e - R_EXo_glyb_e - R_GLYBabc  = 0
 M_glyc1p_c:          R_G1PDH - R_GGGPS - R_CDGGGPP4 - R_CDGGGPP3  = 0
 M_glyc__R_c:         R_AKSCAL  = 0
 M_glyc_c:            - R_GLYCt  = 0
 M_glyc_e:            R_GLYCt - R_EXo_glyc_e + R_EXi_glyc_e  = 0
 M_glyclt_c:          R_PGLYCP  = 0
 M_glycogen_c:        - R_GLCP - 0.01543 R_BIOMASS_Mb_30 + R_GLCS2 - R_GLCGSD
                       = 0
 M_glytrna_c:         R_GLYTRS  = 0
 M_gmp_c:             R_LPFPLT + R_CF3Hg + R_GUAPRT + R_GMPS + R_GMPS2 - R_GK1
                      + R_ADOCBLS2  = 0
 M_grdp_c:            - R_GRTT + R_DMATT  = 0
 M_gthox_c:           - R_GTHOr  = 0
 M_gthrd_c:           2 R_GTHOr  = 0
 M_gtp_c:             - R_H4MPTGL_gtp - R_LPPGS + R_NDPK1 - R_GF4GL_6
                      - R_GF4GL_5 - R_GF4GL_4 - R_GF4GL_3 - R_GF4GL_2
                      - R_GF4GL_1 - R_GF4GL_0 - R_RNTR2 - R_ACBIPGT
                      - 0.2222 R_BIOMASS_Mb_30 - R_ADSS - R_NTP3 - R_ADK3
                      - R_GTPCIII - R_2PLS - R_CF3Sg  = 0
 M_gua_c:             - R_GUAPRT - R_GUAD  = 0
 M_h2_c:              - R_F4RHi - R_H2td + R_NIT_n1p4 - R_TSR + R_ECHH_10  = 0
 M_h2_e:              - R_EXo_h2_e + R_H2td + R_EXi_h2_e - R_F4NH  = 0
 M_h2acon_C_c:        - R_CH2ACH + R_RH2CD  = 0
 M_h2mpt_c:           R_H4MPTS16 - R_H2MPTR  = 0
 M_h2o2_c:            - 2 R_CAT + R_SPODM  = 0
 M_h2o_c:             - R_CTPS2 - R_ALDD2x - R_PRPPP - R_ASNN - R_METGL
                      - R_DOHDUS - R_CF3Hg - R_CF3Ha + R_DHQTi - R_ACP1_FMN
                      + 2 R_CAT - R_PPA_1 - R_CHACH - R_CA2abc - R_ATPM
                      + R_DHAD2 - R_PRFGS - R_NTP7 - R_ADOCBIAH - R_H4MPTS10
                      - R_GLNTRAT - R_PRAMPC + R_G5SADs + 2 R_NTRIR2y
                      + 2 R_NTRIR2x - R_CBL1abc - R_CH2ACH - R_DGGPGP2
                      - R_ATRZCH - R_AHXDH + R_SHCD - R_DGGPGP - R_COMS
                      - R_IMPD - R_ZNabc + R_ENO + R_RNTR4 + R_RNTR2 + R_RNTR3
                      + R_RNTR1 - R_AHCi - R_TSULabc + R_H2Ot - R_SSALy
                      - R_SSALx + R_UPP3S + R_IGPS - R_HMBS - R_ARGDr + R_TRPS2
                      + R_TRPS1 - R_THACH + R_GLUDyi - R_HMGCOASi - R_DRTPPD
                      - R_GMPS2 - R_HCITS - R_FE2abc + R_SMTZS - R_UDPGD
                      - R_THRS - R_PROabc - R_NADS2 - R_AGMT - R_AKPCAL
                      - R_AKACAL - R_IMPC - R_DHORTS - 30 R_BIOMASS_Mb_30
                      + R_FMFD_b - R_DHPM1 - R_HKt + R_PPNDH - R_CYSDS
                      - 2 R_F430S1 - R_CODHr + R_F430S5 - R_METAT + R_RH3CD
                      + R_RNDR4 + R_RNDR2 + R_RNDR3 - R_COabc + R_RNDR1
                      - 32 R_NIT_n1p4 - R_ALDD1 - R_HCO3E + 2 R_RBFSa - R_PPS
                      - R_CYRDAS - R_FUM - 4 R_ADCYRS - R_PPA - R_CH3ACH
                      - R_CD2abc2 - R_CD2abc1 - R_NTP5 - R_NTP4 - R_NTP6
                      - R_NTP1 - R_NTP3 - R_NTP2 - R_NTP9 - R_NTP8 - R_MOHMT
                      - R_DNMPPA - R_TAPTPP - R_ATPS1 + R_THZPSN - R_ALKP
                      + 2 R_PPBNGS - R_CODH2r - R_AGAID + R_GMAND - R_GLUDxi
                      - 3 R_AKP1 - R_CBL1HBIabc - R_DHPM2 - R_CDPGH - R_PIabc
                      + R_IGPDH - R_PMDPHT - R_Cuabc - R_GLCGSD - R_FEDCabc
                      - R_NaKt_1 - R_PGLYCP - R_GUAD - R_NIabc + R_MCMMT
                      - R_RZ5PP2 + R_ATPS4r - R_SULabc - R_CBIabc - R_MTHFC
                      - R_COCHL - R_ACYP_2 - R_MG2abc + R_IPPMIa - R_IPPMIb
                      + R_GHMT2r + R_DHAD1 + 2 R_QULNS - R_FE3abc - R_DMHDRFS
                      + R_MTSPC - R_MI1PP + 2 R_F4H2O - R_CS - R_CBPS - R_ADD
                      - R_MNabc - R_DHPCPH - R_IPPS + 2 R_DHDPS + R_TAPTPPM
                      - R_DPGase - R_ASNS1 - 2 R_OPAH - R_AKSCAL - R_FBP
                      - R_MOBDabc + R_TDPGDH - R_DCTPD - R_GLYBabc - R_Cut1
                      - R_HISTP - R_THMabc - R_HISTD - R_ACYP - R_RBPC
                      - 2 R_GTPCIII - R_UAG2EMA + R_YUMPS - R_GLUPRT + R_RH2CD
                      - R_PSP_L - R_MPML - R_CSND - R_NBAHH_ir - R_PSLP
                      - R_DCMPDA - R_CPC5MT - R_PRATPP  = 0
 M_h2o_e:             R_EXi_h2o_e - R_H2Ot - R_EXo_h2o_e  = 0
 M_h2s_c:             R_H2St + R_CYSDS + R_TSR + R_BTS3r - R_CYSS - R_AHSERL2
                       = 0
 M_h2s_e:             R_EXi_h2s_e - R_H2St - R_EXo_h2s_e  = 0
 M_h3acon_C_c:        R_RH3CD - R_CH3ACH  = 0
 M_h4mpt_c:           - R_H4MPTGL_gtp + R_H2MPTR - R_H4MPTGL_atp  = 0
 M_h4spt_c:           R_H4MPTGL_gtp - R_FMFTSPFT_b - 0.02362 R_BIOMASS_Mb_30
                      - R_CODHr + R_MTSPCMMT + R_H4MPTGL_atp  = 0
 M_h_c:               - 2 R_HMGCOARi + 2 R_CTPS2 + 2 R_CTPS1 + 2 R_ALDD2x
                      + R_PRPPS + R_PRPPP + R_IPMD + R_H4MPTGL_gtp + R_PROt2r
                      + R_PANTS + R_LPFPLT + R_DOHDUS + R_CF3Hg + R_BLAT
                      + R_CF3Ha + R_TMK - R_NNATr + R_CDGGIPT - R_LPPGS
                      - R_7MHS + R_CA2abc - R_G1PDH - R_SHK3Dr + R_PRAGSr
                      - R_G5SD2 + R_POR3 + R_CYSt2r + R_ATPM + R_FOLt + R_PRFGS
                      + 2 R_PAPSR + R_NTP7 - R_GLUDC - R_DKFPS - R_DKFPR
                      + R_H4MPTS16 + R_GLNTRAT - R_NAPRT + R_G5SADs + R_GLUt2r
                      - 5 R_NTRIR2y - 5 R_NTRIR2x + R_F4MTSPD - R_ASP1DC
                      + R_4ABZt2r - R_PYRDC + R_FMETTRS - R_DROPPRy - R_GLUTRR
                      + R_INS2D + R_CBL1abc + R_ASPCT + R_NN5HBPRT - 8 R_ATGH
                      + R_ANS + R_GAPD_nadp + R_MEVK1 + R_ATRZCH + R_PNTK
                      - R_GALUi + R_PRASCSi - 7 R_HATGH - R_ORNDC - R_DHDPRy
                      - R_FMFTSPFT_b + R_HSK - R_SLDxi2 + R_IMPD + R_GF4GL_6
                      + R_GF4GL_5 + R_GF4GL_4 + R_GF4GL_3 + R_GF4GL_2
                      + R_GF4GL_1 + R_GF4GL_0 + R_ZNabc + R_SHCHD2 - R_NMNAT
                      + R_CDGGGPP4 + R_CDGGGPP3 + R_OCBT - R_HASD - R_F4NPR
                      + R_ETHAt6 - R_ACBIPGT + 2 R_GMPS + R_OORr + R_TSULabc
                      + R_TMAt2 + R_NADK + 2 R_SSALy + 2 R_SSALx + R_SO4t2
                      + R_ARGSS + R_PACt2r - R_IGPS + R_HETZK + R_MDH
                      + R_DPCOAK - R_ARGDC + R_PC + R_PFK_adp - R_ASD
                      - R_GLUDyi - R_HISDC + R_HMGCOASi - R_CBLAT2 - R_DRTPPD
                      + 2 R_GMPS2 + R_AGPR + R_HCITS + R_SHKK - R_BRFAPS
                      - R_TRDR + R_FE2abc + R_MDHy + 3 R_UDPGD + R_CMLDC
                      + R_LDH_L - R_DPMVD2 + R_PROabc + R_NADS1 + R_NADS2
                      + R_AKPCAL + R_AKACAL + R_DB4PS + R_AGAIAGT + R_DHORTS
                      + R_HEX7 + R_MMAt - 2 R_SPODM - R_ASADi
                      + 29.99755 R_BIOMASS_Mb_30 + R_KARA1 - R_FMFD_b
                      + R_ACTNt2r + R_DHPM1 + R_PPNCL2 + R_AKGCAL - 2 R_MTHFR3
                      - R_PPNDH - 2 R_AHMMPS - R_F430S3 - 2 R_F430S2
                      + 5 R_F430S1 + R_TMDK1 + 2 R_CODHr - R_F430S5 + R_GAPD
                      + R_COabc + 6 R_NIT_n1p4 - R_PYK + 2 R_ALDD1 + R_HCO3E
                      + R_RIBFLVt2 + R_ACGS + R_IOR - R_GOFUCR + 2 R_PPS
                      + R_CYRDAR + R_CYRDAS - 2 R_NNDPR + 4 R_ADCYRS + R_GLNS
                      + R_PPA - R_PPCDC - R_PSLS + R_IOR3 + R_IOR2 + R_LTHRK
                      - R_SECS + R_ALCD20y + 2 R_ADSS - R_BACCL + R_POR2
                      + R_CD2abc2 + R_CD2abc1 + R_NTP5 + R_NTP4 + R_NTP6
                      + R_NTP1 + R_NTP3 + R_NTP2 + R_NTP9 + R_NTP8 - R_CPC6R
                      - R_ACLS + R_TAPTPP + R_GLCS2 - R_ATPS1 + R_THZPSN
                      + R_ILEt2r - R_DAPDC + R_PPBNGS - R_IND3ACt2r
                      + 2 R_CODH2r + R_UPP3MT - 1.8 R_HDR + R_GLUDxi
                      - 1.8 R_F4D + 2 R_AKP1 + R_CBL1HBIabc + R_VOR - R_GTHOr
                      - R_KARA2 + R_DHPM2 + R_AIRC2 - R_CDPGS - R_GALT
                      - R_TDPDRR + 3 R_DBTS - 2 R_DAPDH - 8 R_ATIH + 2 R_CDPGH
                      + R_PIabc + R_RBK + R_Cuabc + R_APPLP - R_THRPD + R_HMPK1
                      + R_DHFS - R_DHFR + R_FEDCabc + R_LEUt2r + R_ALCD2y
                      + R_NaKt_1 + R_METS - R_MAN1PT2 + R_BTNt2i + R_H4MPTS15
                      - R_AKSDC - R_GUAD + R_NAt3_1 - R_UAPGR - R_CYRDAAT
                      - R_TMPPP + R_NIabc + 3 R_ATPS4r + R_SULabc - R_IPDDI3y
                      - R_IPDDI3x + R_CBIabc + R_ALACt2r + R_MTHFC + R_PGCD
                      + 2 R_COCHL - R_OMCDC + R_ACYP_2 + R_PYRt2 - R_DPR
                      + R_MG2abc + R_VALt2r + R_ACt2r - R_G1PTT - R_CDGGGS
                      + R_CPC6MT + R_PNTOt2 + R_FE3abc + 2 R_DMHDRFS - R_OMPDC
                      + R_COBS - R_MTSPC + R_CDGGGSAT2 + R_CYSS - 7 R_HATSH
                      + R_CS + 2 R_CBPS + R_CPC2MT + R_FBPO + R_AHGDx
                      - 3 R_ECHH_10 - R_ADD + R_LYSt3r + R_MNabc + R_ADCPS2
                      + R_DHPCPH + 3 R_SHCHF + R_IPPS + R_DHDPS + R_ASNS1
                      + R_OPAH - 1.8 R_F4NH + R_7MHTS + R_SLDx + R_SLDy
                      + R_2PGK - R_MDHi2 + R_MOBDabc - R_DCTPD + R_GLYBabc
                      - R_DROPPRx + R_Cut1 + R_IG3PS + R_ADNK1 + R_Kt2r_copy1
                      + R_Kt2r_copy2 + R_THMabc + 3 R_HISTD - 8 R_ATSH
                      + R_CDGGIPT2 - R_CDGGGS2 + R_ACYP - 2 R_P5CRx + R_HPACt2r
                      + R_H4MPTGL_atp + R_GK_adp - R_OXADC + R_SO3t + R_HSDy
                      + 2 R_RBPC - R_CBIAT + R_GTPCIII + R_UAG2EMA - R_AOXSr
                      + R_GLCNt2r + R_G1PACT + R_CPC4MT + R_2PLS - R_UAGDP
                      - R_ACHBS + R_DMAt + R_HPPK2 - 7 R_HATIH + 3 R_MPML
                      - R_PTPATi + R_ADOCBLS2 - R_SPDC - R_CSND + 2 R_PRAIS
                      - R_DCMPDA + 3 R_SHCHCC + R_PPDK + R_CDGGGSAT - R_ACLDC
                      + R_CPC5MT + R_GARFT + R_DIPS + R_PRATPP + R_AHSERL2  = 0
 M_h_e:               - R_PROt2r + R_PPA_1 - R_CYSt2r - R_FOLt - R_GLUt2r
                      - R_4ABZt2r - R_ETHAt6 - R_TMAt2 - R_SO4t2 - R_PACt2r
                      + R_EXi_h_e - R_MMAt - R_ACTNt2r + R_HKt - R_RIBFLVt2
                      + 2 R_ATPS1 - R_ILEt2r + R_IND3ACt2r + 1.8 R_HDR
                      + 1.8 R_F4D - R_LEUt2r - R_BTNt2i - R_NAt3_1 - 4 R_ATPS4r
                      - R_ALACt2r - R_PYRt2 - R_VALt2r - R_ACt2r - R_PNTOt2
                      + R_ECHH_10 - R_LYSt3r + 1.8 R_F4NH - R_Kt2r_copy1
                      - R_Kt2r_copy2 - R_HPACt2r - R_SO3t - R_GLCNt2r - R_DMAt
                      - R_EXo_h_e  = 0
 M_hacon_C_c:         - R_CHACH + R_SHCD  = 0
 M_hacon_T_c:         - R_THACH + R_AKGCAL  = 0
 M_hatrz_c:           R_ATRZCH  = 0
 M_hcit_c:            R_HCITS  = 0
 M_hco3_c:            - R_PC + R_HCO3E - R_AIRC2 - R_CBPS  = 0
 M_hcys__L_c:         R_AHCi - R_METS + R_AHSERL2  = 0
 M_hgbam_c:           - R_COCHL  = 0
 M_his__L_c:          - R_HISDC - 0.1037 R_BIOMASS_Mb_30 + R_HISTD - R_HISTRS
                      + R_NBAHH_ir  = 0
 M_hisp_c:            R_HSTPT - R_HISTP  = 0
 M_hista_c:           R_HISDC  = 0
 M_histd_c:           R_HISTP - R_HISTD  = 0
 M_histrna_c:         R_HISTRS  = 0
 M_hmbil_c:           - R_UPP3S + R_HMBS  = 0
 M_hmgcoa_c:          - R_HMGCOARi + R_HMGCOASi  = 0
 M_hom__L_c:          - R_HSK - R_HSERTA - R_HSDy  = 0
 M_hphaccoa_c:        R_PACCOAL2 + R_IOR3  = 0
 M_hsfd_c:            R_MCR - R_HDR  = 0
 M_hspmd_c:           - 0.0047 R_BIOMASS_Mb_30 + R_HSPMS  = 0
 M_hxan_c:            - R_HXPRT + R_ADD  = 0
 M_iasp_c:            - R_QULNS  = 0
 M_ibcoa_c:           R_VOR  = 0
 M_icit_c:            R_ACONT - R_ICDHyr  = 0
 M_idp_c:             - R_NDPK9 + R_ADK4  = 0
 M_ihcit_T_c:         R_CHACH - R_TIHCD  = 0
 M_ile__L_c:          - 0.3179 R_BIOMASS_Mb_30 + R_ILEt2r - R_ILETRS - R_ILETA
                       = 0
 M_ile__L_e:          R_EXi_ile__L_e - R_ILEt2r - R_EXo_ile__L_e  = 0
 M_iletrna_c:         R_ILETRS  = 0
 M_imacp_c:           - R_HSTPT + R_IGPDH  = 0
 M_imp_c:             - R_IMPD - R_IMPC - R_ADSS + R_HXPRT  = 0
 M_ind3ac_c:          - R_PACCOAL3 - R_IND3ACt2r  = 0
 M_ind3ac_e:          - R_EXo_ind3ac_e + R_IND3ACt2r + R_EXi_ind3ac_e  = 0
 M_indaccoa_c:        R_PACCOAL3 + R_IOR  = 0
 M_indole_c:          R_TRPS3 - R_TRPS2  = 0
 M_indpyr_c:          R_TRPTA - R_IOR  = 0
 M_inost_c:           - R_CDGGIPT - R_INS2D + R_MI1PP - R_CDGGIPT2  = 0
 M_ipdp_c:            R_DPMVD - R_GRTT - R_GRTT2 - R_FRTT - R_FRTT2 - R_IPDDI3y
                      - R_IPDDI3x - R_DMATT2 - R_DMATT  = 0
 M_itp_c:             R_NDPK9 - R_ADK4  = 0
 M_k_c:               R_HKt + R_NaKt_1 + R_Kt2r_copy1 + R_Kt2r_copy2  = 0
 M_k_e:               - R_HKt + R_EXi_k_e - R_NaKt_1 - R_EXo_k_e - R_Kt2r_copy1
                      - R_Kt2r_copy2  = 0
 M_lac__L_c:          - R_LDH_L - R_2PLS  = 0
 M_lald__L_c:         R_FCLPA  = 0
 M_leu__L_c:          - 0.493 R_BIOMASS_Mb_30 - R_LEUTRS + R_LEUt2r - R_LEUTA
                       = 0
 M_leu__L_e:          R_EXi_leu__L_e - R_LEUt2r - R_EXo_leu__L_e  = 0
 M_leutrna_c:         R_LEUTRS  = 0
 M_lppg_c:            - R_LPFPLT + R_LPPGS  = 0
 M_lys__L_c:          - 0.3755 R_BIOMASS_Mb_30 + R_DAPDC - R_LYSt3r - R_LYSTRS
                      - R_LYSAM  = 0
 M_lys__L_e:          R_LYSt3r + R_EXi_lys__L_e - R_EXo_lys__L_e  = 0
 M_lystrna_c:         R_LYSTRS  = 0
 M_m3hdp_c:           R_DPMVD2 - R_DMATT2  = 0
 M_mal__L_c:          - R_MDH - R_MDHy + R_FUM - R_ME1 + R_MDHi2  = 0
 M_man1p_c:           - R_MAN1PT2 - R_PMANM  = 0
 M_man6p_c:           R_PMANM - R_MAN6PI  = 0
 M_mcom_c:            - R_MCR + R_DMAMT + R_MTCMMT + R_MCMMT + R_TMAMT
                      + R_MMAMT + R_MTSPCMMT  = 0
 M_menylh4spt_c:      - R_F4MTSPD + R_MTSPC  = 0
 M_meoh_c:            - R_MCMMT + R_MEOHt2  = 0
 M_meoh_e:            - R_EXo_meoh_e + R_EXi_meoh_e - R_MEOHt2  = 0
 M_met__L_c:          - R_METGL - 0.1682 R_BIOMASS_Mb_30 - R_METAT + R_METS
                      - R_METTRS  = 0
 M_methf_c:           - R_MTHFC + R_MTHFD2i  = 0
 M_mettrna_c:         - R_FMETTRS + R_METTRS  = 0
 M_mev__R_c:          R_HMGCOARi - R_MEVK1  = 0
 M_mfr_b_c:           R_FMFTSPFT_b - R_FMFD_b  = 0
 M_mg2_c:             - R_MGt5 + R_MG2abc - R_MPML  = 0
 M_mg2_e:             R_EXi_mg2_e + R_MGt5 - R_EXo_mg2_e - R_MG2abc  = 0
 M_mh4spt_c:          R_F4MTSPR + R_CODHr - R_MTSPCMMT  = 0
 M_mi1p__D_c:         R_MI1PS - R_MI1PP  = 0
 M_mleneh4spt_c:      R_F4MTSPD - R_F4MTSPR  = 0
 M_mlthf_c:           - R_TMDS - R_MTHFR3 - R_MOHMT + R_GHMT2r - R_MTHFD2i  = 0
 M_mma_c:             R_DMAMT + R_MMAt - R_MMAMT  = 0
 M_mma_e:             - R_MMAt + R_EXi_mma_e - R_EXo_mma_e  = 0
 M_mmh2mpt_c:         - R_H4MPTS16 + R_H4MPTS15  = 0
 M_mn2_c:             R_MNabc  = 0
 M_mn2_e:             R_EXi_mn2_e - R_EXo_mn2_e - R_MNabc  = 0
 M_mobd_c:            R_MOBDabc  = 0
 M_mobd_e:            R_EXi_mobd_e - R_MOBDabc - R_EXo_mobd_e  = 0
 M_mphen_c:           R_HDR - R_F4D - R_F4NH  = 0
 M_mphenh2_c:         - R_HDR + R_F4D + R_F4NH  = 0
 M_mppp9_c:           R_MPML  = 0
 M_n2_c:              R_N2tr - R_NIT_n1p4  = 0
 M_n2_e:              - R_EXo_n2_e - R_N2tr + R_EXi_n2_e  = 0
 M_na1_c:             R_CAt6 + R_GLYt4 + R_ALAt4 - R_NaKt_1 - R_NAt3_1
                      - 1.7 R_MTSPCMMT  = 0
 M_na1_e:             - R_EXo_na1_e - R_CAt6 - R_GLYt4 + R_EXi_na1_e - R_ALAt4
                      + R_NaKt_1 + R_NAt3_1 + 1.7 R_MTSPCMMT  = 0
 M_nabl_c:            R_BLAT  = 0
 M_nac_c:             - R_NAPRT + R_NN5HBPRT + R_NACt  = 0
 M_nac_e:             - R_EXo_nac_e + R_EXi_nac_e - R_NACt  = 0
 M_nad_c:             - R_ALDD2x - R_IPMD - R_DOHDUS + R_7MHS + R_G1PDH
                      + R_G5SD2 + R_DKFPS + R_DKFPR + 3 R_NTRIR2x - R_INS2D
                      + R_SLDxi2 - R_IMPD - R_SHCHD2 - R_TIH2CD + R_NMNAT
                      - R_FDH - R_NADK - R_SSALx - R_TIH3CD - R_MDH - 2 R_UDPGD
                      - R_LDH_L + R_NADS1 + R_NADS2 - 0.00215 R_BIOMASS_Mb_30
                      + R_F430S3 + R_F430S2 - R_GAPD - R_ALDD1 + R_CYRDAR
                      + R_SECS - R_GLUDxi - R_PPND + R_IPDDI3x - R_PGCD - R_ME1
                      - R_FBPO - R_AHGDx - R_MTHFD2i - R_SLDx + R_MDHi2
                      + R_DROPPRx - R_TIHCD - 2 R_HISTD + R_P5CRx  = 0
 M_nadh_c:            R_ALDD2x + R_IPMD + R_DOHDUS - R_7MHS - R_G1PDH - R_G5SD2
                      - R_DKFPS - R_DKFPR - 3 R_NTRIR2x + R_INS2D - R_SLDxi2
                      + R_IMPD + R_SHCHD2 + R_TIH2CD + R_FDH + R_SSALx
                      + R_TIH3CD + R_MDH + 2 R_UDPGD + R_LDH_L
                      - 0.00005 R_BIOMASS_Mb_30 - R_F430S3 - R_F430S2 + R_GAPD
                      + R_ALDD1 - R_CYRDAR - R_SECS + R_GLUDxi + R_PPND
                      - R_IPDDI3x + R_PGCD + R_ME1 + R_FBPO + R_AHGDx
                      + R_MTHFD2i + R_SLDx - R_MDHi2 - R_DROPPRx + R_TIHCD
                      + 2 R_HISTD - R_P5CRx  = 0
 M_nadp_c:            2 R_HMGCOARi + R_SHK3Dr + 3 R_NTRIR2y + R_DROPPRy
                      + R_GLUTRR + 8 R_ATGH - R_GAPD_nadp + 7 R_HATGH
                      + R_DHDPRy + R_F4NPR + R_NADK - R_SSALy + R_GLUDyi
                      - R_AGPR + R_TRDR - R_MDHy + R_ASADi
                      - 0.00013 R_BIOMASS_Mb_30 - R_KARA1 + R_MTHFR3 + R_GOFUCR
                      - R_ALCD20y + R_FOLR2 + R_CPC6R + R_GTHOr + R_KARA2
                      + R_TDPDRR + R_DAPDH + 8 R_ATIH + R_DHFR - R_ALCD2y
                      - R_ICDHyr + R_UAPGR + R_IPDDI3y + R_DPR - 2 R_DMHDRFS
                      + 7 R_HATSH - R_SLDy + 8 R_ATSH - R_HSDy + 7 R_HATIH  = 0
 M_nadph_c:           - 2 R_HMGCOARi - R_SHK3Dr - 3 R_NTRIR2y - R_DROPPRy
                      - R_GLUTRR - 8 R_ATGH + R_GAPD_nadp - 7 R_HATGH
                      - R_DHDPRy - R_F4NPR + R_SSALy - R_GLUDyi + R_AGPR
                      - R_TRDR + R_MDHy - R_ASADi - 0.0004 R_BIOMASS_Mb_30
                      + R_KARA1 - R_MTHFR3 - R_GOFUCR + R_ALCD20y - R_FOLR2
                      - R_CPC6R - R_GTHOr - R_KARA2 - R_TDPDRR - R_DAPDH
                      - 8 R_ATIH - R_DHFR + R_ALCD2y + R_ICDHyr - R_UAPGR
                      - R_IPDDI3y - R_DPR + 2 R_DMHDRFS - 7 R_HATSH + R_SLDy
                      - 8 R_ATSH + R_HSDy - 7 R_HATIH  = 0
 M_nh4_c:             - R_CTPS1 + R_ASNN + R_METGL + R_DOHDUS + R_NTRIR2y
                      + R_NTRIR2x + R_ORNCD + R_COMS - R_GMPS + 4 R_HMBS
                      + R_ARGDr - R_GLUDyi + R_DRTPPD - R_NADS1 + R_CYSDS
                      + 2 R_NIT_n1p4 - R_GLNS + R_GLUDxi + R_HSPMS - R_DAPDH
                      + R_GUAD + R_NH4t + R_MMAMT + R_DMHDRFS + R_ADD + R_DCTPD
                      + R_CSND + R_DCMPDA  = 0
 M_nh4_e:             - R_NH4t - R_EXo_nh4_e + R_EXi_nh4_e  = 0
 M_ni2_c:             - R_NIt5 - R_F430S1 + R_NIabc  = 0
 M_ni2_e:             R_NIt5 + R_EXi_ni2_e - R_EXo_ni2_e - R_NIabc  = 0
 M_nicrnt_c:          - R_NNATr + R_NAPRT - R_NN5HBPRT + R_NNDPR  = 0
 M_nmn_c:             - R_NMNAT  = 0
 M_no2_c:             - R_NTRIR2y - R_NTRIR2x  = 0
 M_o2_c:              R_CAT + R_SPODM - R_F4H2O  = 0
 M_o2s_c:             - 2 R_SPODM  = 0
 M_oaa_c:             R_ASPTA + R_MDH + R_PC + R_MDHy - R_CS - R_MDHi2  = 0
 M_ohepa_c:           R_TIH2CD - R_AKPCAL  = 0
 M_ohexa_c:           - R_AKACAL + R_TIHCD  = 0
 M_orn_c:             - R_ORNCD - R_ORNDC - R_OCBT + R_ORNTAC  = 0
 M_orot5p_c:          - R_OMPDC - R_ORPT  = 0
 M_orot_c:            R_DHORD7 + R_ORPT  = 0
 M_oxa_c:             R_DMHDRFS - R_OXADC  = 0
 M_pac_c:             - R_PACCOAL + R_PACt2r  = 0
 M_pac_e:             R_EXi_pac_e - R_PACt2r - R_EXo_pac_e  = 0
 M_pan4p_c:           R_PPCDC - R_PTPATi  = 0
 M_pant__R_c:         - R_PANTS + R_DPR  = 0
 M_pap_c:             R_PAPSR  = 0
 M_paps_c:            - R_PAPSR  = 0
 M_pep_c:             - R_PSCVT + R_ENO - R_PYK + R_PPS - R_PSLS - R_UAGCVT
                      + R_PPDK  = 0
 M_phaccoa_c:         R_PACCOAL + R_IOR2  = 0
 M_phe__L_c:          - R_PHETA1 - 0.2027 R_BIOMASS_Mb_30 - R_PHETRS  = 0
 M_phetrna_c:         R_PHETRS  = 0
 M_phom_c:            R_HSK - R_THRS  = 0
 M_phpyr_c:           R_PHETA1 + R_PPNDH - R_IOR2  = 0
 M_pi_c:              R_CTPS2 + R_CTPS1 + R_PRPPP + R_H4MPTGL_gtp + R_PSCVT
                      + R_ACP1_FMN + 2 R_PPA_1 + R_CA2abc + R_PRAGSr + R_G5SD2
                      + R_ATPM + R_PRFGS + R_NTP7 + R_H4MPTS10 + R_GLNTRAT
                      + R_CBL1abc + R_ASPCT + R_DGGPGP2 - R_GAPD_nadp - R_GLCP
                      + R_PRASCSi + R_DPMVD + R_DGGPGP + R_ACCOAL + R_GF4GL_6
                      + R_GF4GL_5 + R_GF4GL_4 + R_GF4GL_3 + R_GF4GL_2
                      + R_GF4GL_1 + R_GF4GL_0 + R_ZNabc - R_PTAr + R_OCBT
                      + R_TSULabc + R_PC - R_DURIPP - R_AGPR + R_FE2abc
                      + R_THRS + R_PROabc + R_ASADi + 30 R_BIOMASS_Mb_30
                      + R_4KFBPP + R_HKt - R_TMDPP + R_AHMMPS + 2 R_F430S1
                      - R_GAPD + R_METAT + R_COabc + 32 R_NIT_n1p4 + R_RBFSa
                      + R_PPS + 4 R_ADCYRS + R_GLNS + 2 R_PPA + R_ADSS
                      + R_CD2abc2 + R_CD2abc1 + R_NTP5 + R_NTP4 + R_NTP6
                      + R_NTP1 + R_NTP3 + R_NTP2 + R_NTP9 + R_NTP8 + R_DNMPPA
                      + R_ATPS1 + R_ALKP + 3 R_AKP1 + R_CBL1HBIabc + R_AIRC2
                      + R_CDPGS + R_DBTS + 2 R_PIabc + R_PMDPHT + R_Cuabc
                      + R_DHFS + R_FEDCabc + R_NaKt_1 + R_MAN1PT2 + R_PGLYCP
                      + R_NIabc + R_RZ5PP2 - R_ATPS4r + R_SULabc + R_UAGCVT
                      + R_CBIabc + R_COCHL + R_ACYP_2 + R_MG2abc + R_QULNS
                      + R_FE3abc + R_MI1PP + R_CBPS + R_MNabc + R_ADCPS2
                      + R_DPGase + R_OPAH + R_7MHTS + R_FBP + R_MOBDabc
                      + R_GLYBabc + R_Cut1 + R_HISTP + R_THMabc + R_ACYP
                      + R_H4MPTGL_atp - R_MTAP + R_PSP_L + R_MPML + R_CHORS
                      + R_PRAIS + R_PSLP - R_PPDK  = 0
 M_pi_e:              R_EXi_pi_e - R_PIabc - R_EXo_pi_e  = 0
 M_pmcoa_c:           R_AKSCAL - R_AOXSr  = 0
 M_pnto__R_c:         R_PANTS - R_PNTK + R_PNTOt2  = 0
 M_pnto__R_e:         - R_EXo_pnto__R_e + R_EXi_pnto__R_e - R_PNTOt2  = 0
 M_ppa_c:             - R_ACS2 - R_ACCOAL + R_PPAKr  = 0
 M_ppap_c:            - R_PPAKr  = 0
 M_ppbng_c:           - 4 R_HMBS + R_PPBNGS  = 0
 M_ppcoa_c:           R_ACS2 + R_POR3 + R_ACCOAL  = 0
 M_pphn_c:            - R_PPNDH + R_CHORM - R_PPND  = 0
 M_ppi_c:             R_PANTS - R_PPA_1 + R_SERTRS + R_NNATr - R_PPK2r
                      + R_LPPGS + R_ACS2 + R_GGGPS + R_GUAPRT + R_H4MPTS11
                      + R_NAPRT + R_GLUTRS + R_GALUi + R_PACCOAL + R_ADNCYC
                      + R_PROTRS + R_GRTT + R_DGGGPS2 + R_DHPS3 + R_DHPS2
                      + R_NMNAT + R_ACBIPGT + R_GMPS + R_GRTT2 + R_ARGSS
                      + R_GMPS2 + R_BRFAPS + R_ASPTRS + R_TRPTRS + R_FRTT
                      + R_NADS1 + R_NADS2 + 0.799 R_BIOMASS_Mb_30 + R_PACCOAL2
                      + R_PACCOAL3 + R_PPNCL2 + R_METAT + R_CYRDAS + R_NNDPR
                      - R_PPA + R_ALATRS + R_BACCL + R_TAPTPP + R_ARGTRS
                      + R_THZPSN + R_FRTT2 + R_GALT + R_ILETRS + R_LEUTRS
                      + R_GLYTRS + R_GLUTRS_Gln + R_VALTRS + R_H4MPTS14
                      + R_METTRS + R_TMPPP + R_ANPRT + R_XPPT + R_G1PTT
                      + R_CDGGGS + R_HXPRT - R_ORPT + R_ADK2 + R_TAPTPPM
                      + R_ASNS1 + R_CYSTRS + R_ADPT + R_DMATT2 + R_DMATT
                      + R_TYRTRS + R_CDGGGS2 + R_ATPPRT + R_LYSTRS + R_GLUPRT
                      + R_UAGDP + R_PTPATi + R_HISTRS + R_THRTRS + R_CF3Sg
                      + R_CF3Sa + R_PPDK + R_PHETRS + R_PRATPP + R_DGGGPS  = 0
 M_ppp9_c:            - R_MPML  = 0
 M_pppi_c:            R_PTHPS + R_PPK2r + R_CBLAT2 + R_CYRDAAT - R_ADK2
                      + R_CBIAT  = 0
 M_pram_c:            - R_PRAGSr + R_GLUPRT  = 0
 M_pran_c:            R_ANPRT - R_PRAIi  = 0
 M_prbamp_c:          - R_PRAMPC + R_PRATPP  = 0
 M_prbatp_c:          R_ATPPRT - R_PRATPP  = 0
 M_prfp_c:            R_PRAMPC - R_PRMICI  = 0
 M_prlp_c:            R_PRMICI - R_IG3PS  = 0
 M_pro__L_c:          R_PROt2r + R_ORNCD - R_PROTRS + R_PROabc
                      - 0.2419 R_BIOMASS_Mb_30 + R_P5CRx  = 0
 M_pro__L_e:          - R_PROt2r - R_EXo_pro__L_e - R_PROabc + R_EXi_pro__L_e
                       = 0
 M_protrna_c:         R_PROTRS  = 0
 M_prpp_c:            R_PRPPS - R_PRPPP - R_GUAPRT - R_H4MPTS11 - R_NAPRT
                      - R_BRFAPS - R_NNDPR - R_ANPRT - R_XPPT - R_HXPRT
                      + R_ORPT - R_ADPT - R_ATPPRT - R_GLUPRT  = 0
 M_psd5p_c:           R_YUMPS  = 0
 M_pser__L_c:         R_PSERT - R_PSP_L  = 0
 M_ptp_c:             - R_TAPTPP - R_TAPTPPM + R_GTPCIII  = 0
 M_ptrc_c:            R_ORNDC + R_AGMT - 0.0262 R_BIOMASS_Mb_30 - 2 R_HSPMS
                       = 0
 M_pyr_c:             - R_PYRDC + R_ANS + R_COMS - R_PC + R_LDH_L + R_CYSDS
                      + R_PYK - R_PPS - R_POR2 - 2 R_ACLS + R_PYRt2 + R_ME1
                      - R_DHDPS + R_ALATA_L - R_ACHBS - R_PPDK  = 0
 M_pyr_e:             R_EXi_pyr_e - R_EXo_pyr_e - R_PYRt2  = 0
 M_quln_c:            - R_NNDPR + R_QULNS  = 0
 M_r15bp_c:           R_PRPPP - R_RBPI  = 0
 M_r5hbzi_c:          R_RZ5PP2 - R_ADOCBLS2  = 0
 M_r5p_c:             - R_PRPPS + R_RBK - R_TKT1 - R_RPI - R_YUMPS  = 0
 M_rb15bp_c:          - R_RBPC + R_RBPI  = 0
 M_rib__D_c:          - R_RBK  = 0
 M_ribflv_c:          R_ACP1_FMN + R_RIBFLVt2 + R_RBFSb - R_5HBZIDS  = 0
 M_ribflv_e:          R_EXi_ribflv_e - R_RIBFLVt2 - R_EXo_ribflv_e  = 0
 M_ru5p__D_c:         - R_DB4PS - R_RPE + R_RPI  = 0
 M_s7p_c:             - R_TALA + R_TKT1  = 0
 M_s_c:               R_St - 2 R_BTS3r  = 0
 M_s_e:               - R_EXo_s_e - R_St + R_EXi_s_e  = 0
 M_scl_c:             R_SHCHD2 - R_SHCHF  = 0
 M_sec_c:             - R_COMS + R_SECS  = 0
 M_ser__L_c:          - R_SERTRS + R_7MHS - R_TRPS2 - R_TRPS1
                      - 0.2361 R_BIOMASS_Mb_30 - R_GHMT2r - R_CDGGGSAT2
                      - R_SERAT + R_PSP_L - R_CDGGGSAT  = 0
 M_sertrna_c:         R_SERTRS  = 0
 M_sf430a_c:          - R_F430S5 + R_F430S4  = 0
 M_sheme_c:           R_SHCHF  = 0
 M_skm5p_c:           - R_PSCVT + R_SHKK  = 0
 M_skm_c:             R_SHK3Dr - R_SHKK  = 0
 M_sl__L_c:           R_SLDxi2 - R_SLDx - R_SLDy + R_PSLP  = 0
 M_slp__L_c:          R_PSLS - R_PSLP  = 0
 M_so3_c:             R_PAPSR - R_PSLS + R_TSR + R_SO3t  = 0
 M_so3_e:             - R_EXo_so3_e + R_EXi_so3_e - R_SO3t  = 0
 M_so4_c:             R_SO4t2 + R_SULabc  = 0
 M_so4_e:             - R_SO4t2 - R_SULabc + R_EXi_so4_e - R_EXo_so4_e  = 0
 M_succ_c:            R_SSALy + R_SSALx  = 0
 M_succoa_c:          R_OORr - 0.000003 R_BIOMASS_Mb_30  = 0
 M_sucsal_c:          R_ABTA - R_SSALy - R_SSALx  = 0
 M_thdp_c:            R_DHDPRy - R_DAPDH  = 0
 M_thf_c:             R_FMETTRS - 0.000061 R_BIOMASS_Mb_30 + R_MOHMT + R_AICART
                      + R_DHFR + R_METS - R_GHMT2r + R_GARFT  = 0
 M_thm_c:             - R_TMK + R_THMabc  = 0
 M_thm_e:             - R_EXo_thm_e + R_EXi_thm_e - R_THMabc  = 0
 M_thmmp_c:           R_TMK - R_TMPK + R_TMPPP  = 0
 M_thmpp_c:           R_TMPK  = 0
 M_thr__L_c:          R_THRS - 0.2776 R_BIOMASS_Mb_30 - R_LTHRK - R_7MHTS
                      - R_THRAi - R_THRTRS  = 0
 M_thrp_c:            R_LTHRK - R_THRPD  = 0
 M_thrtrna_c:         R_THRTRS  = 0
 M_thym_c:            R_TMDPP  = 0
 M_thymd_c:           - R_TMDPP - R_TMDK1  = 0
 M_tih2cit_c:         R_CH2ACH - R_TIH2CD  = 0
 M_tih3cit_c:         - R_TIH3CD + R_CH3ACH  = 0
 M_tma_c:             R_TMAt2 - R_TMAMT  = 0
 M_tma_e:             - R_EXo_tma_e - R_TMAt2 + R_EXi_tma_e  = 0
 M_trdox_c:           R_PAPSR + R_RNTR4 + R_RNTR2 + R_RNTR3 + R_RNTR1 - R_TRDR
                      + R_RNDR4 + R_RNDR2 + R_RNDR3 + R_RNDR1  = 0
 M_trdrd_c:           - R_PAPSR - R_RNTR4 - R_RNTR2 - R_RNTR3 - R_RNTR1
                      + R_TRDR - R_RNDR4 - R_RNDR2 - R_RNDR3 - R_RNDR1  = 0
 M_trnaala_c:         - R_ALATRS  = 0
 M_trnaarg_c:         - R_ARGTRS  = 0
 M_trnaasp_c:         - R_ASPTRS  = 0
 M_trnacys_c:         - R_CYSTRS  = 0
 M_trnagln_c:         - R_GLUTRS_Gln  = 0
 M_trnaglu_c:         - R_GLUTRS + R_GLUTRR  = 0
 M_trnagly_c:         - R_GLYTRS  = 0
 M_trnahis_c:         - R_HISTRS  = 0
 M_trnaile_c:         - R_ILETRS  = 0
 M_trnaleu_c:         - R_LEUTRS  = 0
 M_trnalys_c:         - R_LYSTRS  = 0
 M_trnamet_c:         - R_METTRS  = 0
 M_trnaphe_c:         - R_PHETRS  = 0
 M_trnapro_c:         - R_PROTRS  = 0
 M_trnaser_c:         - R_SERTRS  = 0
 M_trnathr_c:         - R_THRTRS  = 0
 M_trnatrp_c:         - R_TRPTRS  = 0
 M_trnatyr_c:         - R_TYRTRS  = 0
 M_trnaval_c:         - R_VALTRS  = 0
 M_trp__L_c:          R_TRPS2 + R_TRPS1 - R_TRPTRS - R_TRPTA
                      - 0.0622 R_BIOMASS_Mb_30  = 0
 M_trptrna_c:         R_TRPTRS  = 0
 M_tsul_c:            R_TSULabc - R_TSR  = 0
 M_tsul_e:            - R_TSULabc + R_EXi_tsul_e - R_EXo_tsul_e  = 0
 M_tyr__L_c:          - R_TYRTA - 0.1509 R_BIOMASS_Mb_30 - R_THZPSN - R_TYRTRS
                       = 0
 M_tyrtrna_c:         R_TYRTRS  = 0
 M_uaccg_c:           - R_UAPGR + R_UAGCVT  = 0
 M_uacgam_c:          - R_AGAIAGT - R_UAG2E - R_UAGCVT - R_UAG2EMA + R_UAGDP
                       = 0
 M_uacmam_c:          R_UAG2E  = 0
 M_uamr_c:            R_UAPGR  = 0
 M_udp_c:             R_NTP7 - R_NDPK2 + R_UMPK_copy2 + R_AGAIAGT - R_RNDR4
                      + R_GLCS2 + R_UMPK_copy1 + R_UAG2EMA  = 0
 M_udpg_c:            R_GALUi - R_UDPGD - R_GLCS2 - R_UDPG4E  = 0
 M_udpgal_c:          R_GALT + R_UDPG4E  = 0
 M_udpglcur_c:        R_UDPGD  = 0
 M_ump_c:             - R_UMPK_copy2 - R_UMPK_copy1 + R_OMPDC  = 0
 M_unknown_cbl1deg_c: R_CBLD + R_UNK_CBL1DEGt  = 0
 M_unknown_cbl1deg_e: - R_EXo_unknown_cbl1deg_e + R_EXi_unknown_cbl1deg_e
                      - R_UNK_CBL1DEGt  = 0
 M_unknown_rbfdeg_c:  R_UNK_RBFDEGt + R_5HBZIDS  = 0
 M_unknown_rbfdeg_e:  R_EXi_unknown_rbfdeg_e - R_EXo_unknown_rbfdeg_e
                      - R_UNK_RBFDEGt  = 0
 M_uppg3_c:           R_UPP3S - R_UPP3MT  = 0
 M_ura_c:             R_DURIPP - R_YUMPS + R_CSND  = 0
 M_urea_c:            R_UREAt + R_AGMT  = 0
 M_urea_e:            R_EXi_urea_e - R_UREAt - R_EXo_urea_e  = 0
 M_utp_c:             - R_CTPS2 - R_CTPS1 - R_NTP7 - R_GALUi + R_NDPK2
                      - R_RNTR4 - 0.1489 R_BIOMASS_Mb_30 - R_GALT - R_UAGDP
                       = 0
 M_val__L_c:          - R_VALTA - 0.4631 R_BIOMASS_Mb_30 - R_VALTRS + R_VALt2r
                       = 0
 M_val__L_e:          R_EXi_val__L_e - R_VALt2r - R_EXo_val__L_e  = 0
 M_valtrna_c:         R_VALTRS  = 0
 M_xan_c:             R_GUAD - R_XPPT  = 0
 M_xmp_c:             R_IMPD - R_GMPS - R_GMPS2 + R_XPPT  = 0
 M_xu5p__D_c:         - R_TKT1 - R_TKT2 + R_RPE  = 0
 M_zn2_c:             R_ZNabc  = 0
 M_zn2_e:             R_EXi_zn2_e - R_ZNabc - R_EXo_zn2_e  = 0
 ind_R_EXi_4abz_e:    bin_R_EXi_4abz_e = 0 -> R_EXi_4abz_e  = 0
 ind_R_EXi_4hphac_e:  bin_R_EXi_4hphac_e = 0 -> R_EXi_4hphac_e  = 0
 ind_R_EXi_ac_e:      bin_R_EXi_ac_e = 0 -> R_EXi_ac_e  = 0
 ind_R_EXi_actn__R_e: bin_R_EXi_actn__R_e = 0 -> R_EXi_actn__R_e  = 0
 ind_R_EXi_ala__L_e:  bin_R_EXi_ala__L_e = 0 -> R_EXi_ala__L_e  = 0
 ind_R_EXi_alac__S_e: bin_R_EXi_alac__S_e = 0 -> R_EXi_alac__S_e  = 0
 ind_R_EXi_btn_e:     bin_R_EXi_btn_e = 0 -> R_EXi_btn_e  = 0
 ind_R_EXi_ca2_e:     bin_R_EXi_ca2_e = 0 -> R_EXi_ca2_e  = 0
 ind_R_EXi_cbi_e:     bin_R_EXi_cbi_e = 0 -> R_EXi_cbi_e  = 0
 ind_R_EXi_cbl1_e:    bin_R_EXi_cbl1_e = 0 -> R_EXi_cbl1_e  = 0
 ind_R_EXi_cbl1hbi_e: bin_R_EXi_cbl1hbi_e = 0 -> R_EXi_cbl1hbi_e  = 0
 ind_R_EXi_cd2_e:     bin_R_EXi_cd2_e = 0 -> R_EXi_cd2_e  = 0
 ind_R_EXi_ch4_e:     bin_R_EXi_ch4_e = 0 -> R_EXi_ch4_e  = 0
 ind_R_EXi_ch4s_e:    bin_R_EXi_ch4s_e = 0 -> R_EXi_ch4s_e  = 0
 ind_R_EXi_cit_e:     bin_R_EXi_cit_e = 0 -> R_EXi_cit_e  = 0
 ind_R_EXi_cl_e:      bin_R_EXi_cl_e = 0 -> R_EXi_cl_e  = 0
 ind_R_EXi_cobalt2_e: bin_R_EXi_cobalt2_e = 0 -> R_EXi_cobalt2_e  = 0
 ind_R_EXi_cu2_e:     bin_R_EXi_cu2_e = 0 -> R_EXi_cu2_e  = 0
 ind_R_EXi_cys__L_e:  bin_R_EXi_cys__L_e = 0 -> R_EXi_cys__L_e  = 0
 ind_R_EXi_dma_e:     bin_R_EXi_dma_e = 0 -> R_EXi_dma_e  = 0
 ind_R_EXi_dms_e:     bin_R_EXi_dms_e = 0 -> R_EXi_dms_e  = 0
 ind_R_EXi_etha_e:    bin_R_EXi_etha_e = 0 -> R_EXi_etha_e  = 0
 ind_R_EXi_fe2_e:     bin_R_EXi_fe2_e = 0 -> R_EXi_fe2_e  = 0
 ind_R_EXi_fe3_e:     bin_R_EXi_fe3_e = 0 -> R_EXi_fe3_e  = 0
 ind_R_EXi_fol_e:     bin_R_EXi_fol_e = 0 -> R_EXi_fol_e  = 0
 ind_R_EXi_gcald_e:   bin_R_EXi_gcald_e = 0 -> R_EXi_gcald_e  = 0
 ind_R_EXi_glcn_e:    bin_R_EXi_glcn_e = 0 -> R_EXi_glcn_e  = 0
 ind_R_EXi_glu__L_e:  bin_R_EXi_glu__L_e = 0 -> R_EXi_glu__L_e  = 0
 ind_R_EXi_gly_e:     bin_R_EXi_gly_e = 0 -> R_EXi_gly_e  = 0
 ind_R_EXi_glyald_e:  bin_R_EXi_glyald_e = 0 -> R_EXi_glyald_e  = 0
 ind_R_EXi_glyb_e:    bin_R_EXi_glyb_e = 0 -> R_EXi_glyb_e  = 0
 ind_R_EXi_glyc_e:    bin_R_EXi_glyc_e = 0 -> R_EXi_glyc_e  = 0
 ind_R_EXi_h2s_e:     bin_R_EXi_h2s_e = 0 -> R_EXi_h2s_e  = 0
 ind_R_EXi_ile__L_e:  bin_R_EXi_ile__L_e = 0 -> R_EXi_ile__L_e  = 0
 ind_R_EXi_ind3ac_e:  bin_R_EXi_ind3ac_e = 0 -> R_EXi_ind3ac_e  = 0
 ind_R_EXi_k_e:       bin_R_EXi_k_e = 0 -> R_EXi_k_e  = 0
 ind_R_EXi_leu__L_e:  bin_R_EXi_leu__L_e = 0 -> R_EXi_leu__L_e  = 0
 ind_R_EXi_lys__L_e:  bin_R_EXi_lys__L_e = 0 -> R_EXi_lys__L_e  = 0
 ind_R_EXi_meoh_e:    bin_R_EXi_meoh_e = 0 -> R_EXi_meoh_e  = 0
 ind_R_EXi_mg2_e:     bin_R_EXi_mg2_e = 0 -> R_EXi_mg2_e  = 0
 ind_R_EXi_mma_e:     bin_R_EXi_mma_e = 0 -> R_EXi_mma_e  = 0
 ind_R_EXi_mn2_e:     bin_R_EXi_mn2_e = 0 -> R_EXi_mn2_e  = 0
 ind_R_EXi_mobd_e:    bin_R_EXi_mobd_e = 0 -> R_EXi_mobd_e  = 0
 ind_R_EXi_n2_e:      bin_R_EXi_n2_e = 0 -> R_EXi_n2_e  = 0
 ind_R_EXi_na1_e:     bin_R_EXi_na1_e = 0 -> R_EXi_na1_e  = 0
 ind_R_EXi_nac_e:     bin_R_EXi_nac_e = 0 -> R_EXi_nac_e  = 0
 ind_R_EXi_nh4_e:     bin_R_EXi_nh4_e = 0 -> R_EXi_nh4_e  = 0
 ind_R_EXi_ni2_e:     bin_R_EXi_ni2_e = 0 -> R_EXi_ni2_e  = 0
 ind_R_EXi_pac_e:     bin_R_EXi_pac_e = 0 -> R_EXi_pac_e  = 0
 ind_R_EXi_pi_e:      bin_R_EXi_pi_e = 0 -> R_EXi_pi_e  = 0
 ind_R_EXi_pnto__R_e: bin_R_EXi_pnto__R_e = 0 -> R_EXi_pnto__R_e  = 0
 ind_R_EXi_pro__L_e:  bin_R_EXi_pro__L_e = 0 -> R_EXi_pro__L_e  = 0
 ind_R_EXi_pyr_e:     bin_R_EXi_pyr_e = 0 -> R_EXi_pyr_e  = 0
 ind_R_EXi_ribflv_e:  bin_R_EXi_ribflv_e = 0 -> R_EXi_ribflv_e  = 0
 ind_R_EXi_s_e:       bin_R_EXi_s_e = 0 -> R_EXi_s_e  = 0
 ind_R_EXi_so3_e:     bin_R_EXi_so3_e = 0 -> R_EXi_so3_e  = 0
 ind_R_EXi_so4_e:     bin_R_EXi_so4_e = 0 -> R_EXi_so4_e  = 0
 ind_R_EXi_urea_e:    bin_R_EXi_urea_e = 0 -> R_EXi_urea_e  = 0
 ind_R_EXi_val__L_e:  bin_R_EXi_val__L_e = 0 -> R_EXi_val__L_e  = 0
Bounds
 0 <= R_EXi_ribflv_e <= 999999
-999999 <= R_GLYCt <= 999999
 0 <= R_HMGCOARi <= 999999
 0 <= R_CTPS2 <= 999999
 0 <= R_CTPS1 <= 999999
 0 <= R_ALDD2x <= 999999
-999999 <= R_CH4t <= 999999
-999999 <= R_PRPPS <= 999999
 0 <= R_PRPPP <= 999999
 0 <= R_EXo_glyald_e <= 999999
 0 <= R_EXo_so3_e <= 999999
 0 <= R_IPMD <= 999999
 0 <= R_EXo_h2_e <= 999999
 0 <= R_EXo_cd2_e <= 999999
 0 <= R_ASNN <= 999999
 0 <= R_H4MPTGL_gtp <= 999999
 0 <= R_METGL <= 999999
 0 <= R_PROt2r <= 999999
 0 <= R_PANTS <= 999999
 0 <= R_EXi_4abz_e <= 999999
-999999 <= R_PSCVT <= 999999
 0 <= R_G6PDH3 <= 999999
-999999 <= R_DMSt <= 999999
 0 <= R_LPFPLT <= 999999
 0 <= R_DOHDUS <= 999999
 0 <= R_CF3Hg <= 999999
 0 <= R_EXo_s_e <= 999999
 0 <= R_PTHPS <= 999999
 0 <= R_BLAT <= 999999
 0 <= R_CF3Ha <= 999999
 0 <= R_PMEVK <= 999999
-999999 <= R_DHQTi <= 999999
-999999 <= R_TMK <= 999999
 0 <= R_EXi_glyb_e <= 999999
 0 <= R_EXi_pyr_e <= 999999
 0 <= R_ACP1_FMN <= 999999
-999999 <= R_FCLPA <= 999999
 0 <= R_CAT <= 999999
 0 <= R_PPA_1 <= 999999
-999999 <= R_CHACH <= 999999
 0 <= R_SERTRS <= 999999
 0 <= R_NNATr <= 999999
 0 <= R_CDGGIPT <= 999999
 0 <= R_EXo_pro__L_e <= 999999
 0 <= R_PPK2r <= 999999
 0 <= R_LPPGS <= 999999
 0 <= R_EXo_n2_e <= 999999
 0 <= R_ACS2 <= 999999
 0 <= R_7MHS <= 999999
 0 <= R_CPC3MT <= 999999
 0 <= R_CA2abc <= 999999
 0 <= R_EXi_mg2_e <= 999999
 0 <= R_EXi_mn2_e <= 999999
-999999 <= R_G1PDH <= 999999
-999999 <= R_SHK3Dr <= 999999
-999999 <= R_PRAGSr <= 999999
 0 <= R_F4RHi <= 999999
 0 <= R_G5SD2 <= 999999
 0 <= R_EXo_btn_e <= 999999
 0 <= R_EXo_4abz_e <= 999999
-999999 <= R_POR3 <= 999999
-999999 <= R_CYSt2r <= 999999
      R_ATPM = 0.8
 0 <= R_EXi_h2s_e <= 999999
-999999 <= R_FOLt <= 999999
 0 <= R_EXo_na1_e <= 999999
 0 <= R_DHAD2 <= 999999
 0 <= R_PRFGS <= 999999
      R_EXi_zn2_e = 0
 0 <= R_PAPSR <= 999999
 0 <= R_GGGPS <= 999999
 0 <= R_EXi_pi_e <= 999999
 0 <= R_H4MPTS9 <= 999999
 0 <= R_GUAPRT <= 999999
 0 <= R_NTP7 <= 999999
 0 <= R_EXo_mn2_e <= 999999
 0 <= R_ADOCBIAH <= 999999
 0 <= R_GLUDC <= 999999
 0 <= R_DKFPS <= 999999
-999999 <= R_DKFPR <= 999999
 0 <= R_H4MPTS16 <= 999999
 0 <= R_H4MPTS10 <= 999999
 0 <= R_H4MPTS11 <= 999999
 0 <= R_H4MPTS12 <= 999999
 0 <= R_H4MPTS13 <= 999999
 0 <= R_GLNTRAT <= 999999
 0 <= R_MCR <= 999999
 0 <= R_PRAMPC <= 999999
-999999 <= R_NAPRT <= 999999
-999999 <= R_G5SADs <= 999999
 0 <= R_EXi_val__L_e <= 999999
-999999 <= R_GLUt2r <= 999999
 0 <= R_NTRIR2y <= 999999
 0 <= R_NTRIR2x <= 999999
-999999 <= R_F4MTSPD <= 999999
 0 <= R_EXi_urea_e <= 999999
 0 <= R_ASP1DC <= 999999
 0 <= R_EXi_leu__L_e <= 999999
-999999 <= R_FBA2 <= 999999
 0 <= R_EXo_actn__R_e <= 999999
-999999 <= R_F4MTSPR <= 999999
-999999 <= R_4ABZt2r <= 999999
 0 <= R_PYRDC <= 999999
 0 <= R_GLUS_F420 <= 999999
 0 <= R_FMETTRS <= 999999
 0 <= R_DROPPRy <= 999999
 0 <= R_GLUTRS <= 999999
 0 <= R_GLUTRR <= 999999
 0 <= R_INS2D <= 999999
 0 <= R_EXi_pac_e <= 999999
 0 <= R_CBL1abc <= 999999
 0 <= R_ASPCT <= 999999
-999999 <= R_CH2ACH <= 999999
 0 <= R_NN5HBPRT <= 999999
 0 <= R_ATGH <= 999999
 0 <= R_ANS <= 999999
 0 <= R_HSTPT <= 999999
 0 <= R_DGGPGP2 <= 999999
 0 <= R_EXi_h2o_e <= 999999
-999999 <= R_GLYALDt <= 999999
 0 <= R_EXo_cys__L_e <= 999999
-999999 <= R_GAPD_nadp <= 999999
 0 <= R_GLU5K <= 999999
      R_EXi_unknown_rbfdeg_e = 0
 0 <= R_ORNCD <= 999999
 0 <= R_MEVK1 <= 999999
 0 <= R_GLCP <= 999999
 0 <= R_EXi_4hphac_e <= 999999
 0 <= R_TMDS <= 999999
 0 <= R_ATRZCH <= 999999
 0 <= R_AHXDH <= 999999
-999999 <= R_SHCD <= 999999
 0 <= R_PNTK <= 999999
 0 <= R_EXo_thm_e <= 999999
 0 <= R_GALUi <= 999999
-999999 <= R_PRASCSi <= 999999
 0 <= R_PACCOAL <= 999999
 0 <= R_HATGH <= 999999
 0 <= R_GDMANE <= 999999
 0 <= R_DPMVD <= 999999
 0 <= R_ORNDC <= 999999
 0 <= R_DHDPRy <= 999999
 0 <= R_ADNCYC <= 999999
 0 <= R_GF6PTA <= 999999
-999999 <= R_FMFTSPFT_b <= 999999
 0 <= R_DGGPGP <= 999999
 0 <= R_PROTRS <= 999999
 0 <= R_HSK <= 999999
 0 <= R_SLDxi2 <= 999999
-999999 <= R_TPI <= 999999
-999999 <= R_ACCOAL <= 999999
 0 <= R_EXi_ala__L_e <= 999999
 0 <= R_EXi_actn__R_e <= 999999
 0 <= R_GRTT <= 999999
-999999 <= R_ASPTA <= 999999
-999999 <= R_NDPK9 <= 999999
-999999 <= R_NDPK8 <= 999999
-999999 <= R_NDPK7 <= 999999
-999999 <= R_NDPK6 <= 999999
-999999 <= R_NDPK5 <= 999999
-999999 <= R_NDPK4 <= 999999
-999999 <= R_NDPK3 <= 999999
-999999 <= R_NDPK2 <= 999999
-999999 <= R_NDPK1 <= 999999
 0 <= R_UMPK_copy2 <= 999999
 0 <= R_COMS <= 999999
 0 <= R_IMPD <= 999999
 0 <= R_DGGGPS2 <= 999999
 0 <= R_GF4GL_6 <= 999999
 0 <= R_GF4GL_5 <= 999999
 0 <= R_GF4GL_4 <= 999999
 0 <= R_GF4GL_3 <= 999999
 0 <= R_GF4GL_2 <= 999999
 0 <= R_GF4GL_1 <= 999999
 0 <= R_GF4GL_0 <= 999999
 0 <= R_EXo_tma_e <= 999999
 0 <= R_EXo_alac__S_e <= 999999
 0 <= R_ZNabc <= 999999
 0 <= R_DHPS3 <= 999999
-999999 <= R_ENO <= 999999
 0 <= R_DHPS2 <= 999999
-999999 <= R_TYRTA <= 999999
 0 <= R_DHNPA2r <= 999999
 0 <= R_SHCHD2 <= 999999
 0 <= R_TIH2CD <= 999999
-999999 <= R_AMAOTr <= 999999
 0 <= R_RNTR4 <= 999999
 0 <= R_RNTR2 <= 999999
 0 <= R_RNTR3 <= 999999
 0 <= R_RNTR1 <= 999999
 0 <= R_NMNAT <= 999999
 0 <= R_EXo_cbl1hbi_e <= 999999
 0 <= R_CDGGGPP4 <= 999999
 0 <= R_CDGGGPP3 <= 999999
 0 <= R_EXo_ala__L_e <= 999999
 0 <= R_ADSL1r <= 999999
      R_EXi_co_e = 0
-999999 <= R_PTAr <= 999999
-999999 <= R_OCBT <= 999999
 0 <= R_HASD <= 999999
 0 <= R_EXi_gcald_e <= 999999
-999999 <= R_F4NPR <= 999999
-999999 <= R_ETHAt6 <= 999999
 0 <= R_ACBIPGT <= 999999
 0 <= R_EXo_dms_e <= 999999
 0 <= R_GMPS <= 999999
-999999 <= R_CO2t <= 999999
 0 <= R_EXi_ch4s_e <= 999999
-999999 <= R_OORr <= 999999
 0 <= R_ABTA <= 999999
-999999 <= R_AHCi <= 999999
-999999 <= R_H2td <= 999999
-999999 <= R_CH4St <= 999999
 0 <= R_TSULabc <= 999999
 0 <= R_GRTT2 <= 999999
-999999 <= R_TMAt2 <= 999999
-999999 <= R_NIt5 <= 999999
 0 <= R_FDH <= 999999
 0 <= R_NADK <= 999999
-999999 <= R_H2Ot <= 999999
 0 <= R_SSALy <= 999999
 0 <= R_SSALx <= 999999
      R_EXi_cu2_e = 0
 0 <= R_SO4t2 <= 999999
-999999 <= R_MGt5 <= 999999
 0 <= R_EXi_ch4_e <= 999999
-999999 <= R_VALTA <= 999999
 0 <= R_ARGSS <= 999999
-999999 <= R_PACt2r <= 999999
 0 <= R_TIH3CD <= 999999
-999999 <= R_PHETA1 <= 999999
 0 <= R_UPP3S <= 999999
 0 <= R_IGPS <= 999999
-999999 <= R_UREAt <= 999999
 0 <= R_HETZK <= 999999
-999999 <= R_MDH <= 999999
 0 <= R_EXo_co2_e <= 999999
 0 <= R_DPCOAK <= 999999
 0 <= R_ARGDC <= 999999
 0 <= R_EXi_gly_e <= 999999
 0 <= R_EXo_unknown_rbfdeg_e <= 999999
-999999 <= R_PGM <= 999999
 0 <= R_HMBS <= 999999
-999999 <= R_PGK <= 999999
-999999 <= R_CAt6 <= 999999
-999999 <= R_PGI <= 999999
 0 <= R_EXo_urea_e <= 999999
-999999 <= R_ARGDr <= 999999
 0 <= R_PC <= 999999
 0 <= R_TRPS3 <= 999999
 0 <= R_TRPS2 <= 999999
 0 <= R_TRPS1 <= 999999
 0 <= R_EXi_cys__L_e <= 999999
      R_EXi_cd2_e = 0
 0 <= R_PFK_adp <= 999999
 0 <= R_ASD <= 999999
-999999 <= R_THACH <= 999999
 0 <= R_GLUDyi <= 999999
 0 <= R_HISDC <= 999999
 0 <= R_HMGCOASi <= 999999
 0 <= R_EXo_glu__L_e <= 999999
-999999 <= R_DURIPP <= 999999
-999999 <= R_St <= 999999
 0 <= R_CBLAT2 <= 999999
 0 <= R_EXo_4hphac_e <= 999999
 0 <= R_DRTPPD <= 999999
 0 <= R_GMPS2 <= 999999
-999999 <= R_AGPR <= 999999
 0 <= R_HCITS <= 999999
 0 <= R_SHKK <= 999999
 0 <= R_BRFAPS <= 999999
 0 <= R_TRDR <= 999999
 0 <= R_FE2abc <= 999999
 0 <= R_SMTZS <= 999999
-999999 <= R_MDHy <= 999999
 0 <= R_DMAMT <= 999999
 0 <= R_EXo_meoh_e <= 999999
      R_EXi_tsul_e = 0
 0 <= R_PRMICI <= 999999
-999999 <= R_UDPGD <= 999999
 0 <= R_EXi_s_e <= 999999
 0 <= R_EXo_ch4s_e <= 999999
 0 <= R_CMLDC <= 999999
 0 <= R_ASPTRS <= 999999
 0 <= R_EXi_ni2_e <= 999999
 0 <= R_THRS <= 999999
-999999 <= R_LDH_L <= 999999
 0 <= R_DPMVD2 <= 999999
 0 <= R_TRPTRS <= 999999
 0 <= R_EXi_alac__S_e <= 999999
 0 <= R_FRTT <= 999999
-999999 <= R_TRPTA <= 999999
 0 <= R_PROabc <= 999999
 0 <= R_NADS1 <= 999999
 0 <= R_NADS2 <= 999999
 0 <= R_AGMT <= 999999
 0 <= R_PSERT <= 999999
 0 <= R_AKPCAL <= 999999
 0 <= R_AKACAL <= 999999
 0 <= R_DB4PS <= 999999
 0 <= R_EXo_gly_e <= 999999
      R_EXi_h_e = 0
-999999 <= R_H2St <= 999999
-999999 <= R_IMPC <= 999999
 0 <= R_AGAIAGT <= 999999
-999999 <= R_DHORTS <= 999999
-999999 <= R_ACONT <= 999999
 0 <= R_EXo_cbl1_e <= 999999
 0 <= R_HEX7 <= 999999
-999999 <= R_MMAt <= 999999
 0 <= R_SPODM <= 999999
 0 <= R_ASADi <= 999999
 0 <= R_EXo_cu2_e <= 999999
 1 <= R_BIOMASS_Mb_30 <= 1000
-999999 <= R_KARA1 <= 999999
 0 <= R_EXo_fol_e <= 999999
-999999 <= R_FMFD_b <= 999999
 0 <= R_PACCOAL2 <= 999999
 0 <= R_PACCOAL3 <= 999999
-999999 <= R_ACTNt2r <= 999999
-999999 <= R_TALA <= 999999
 0 <= R_EXo_mg2_e <= 999999
 0 <= R_DHPM1 <= 999999
 0 <= R_4KFBPP <= 999999
 0 <= R_HKt <= 999999
 0 <= R_PPNCL2 <= 999999
 0 <= R_EXo_pyr_e <= 999999
-999999 <= R_TMDPP <= 999999
 0 <= R_AKGCAL <= 999999
 0 <= R_MTHFR3 <= 999999
 0 <= R_PPNDH <= 999999
 0 <= R_CYSDS <= 999999
 0 <= R_EXo_nac_e <= 999999
 0 <= R_AHMMPS <= 999999
 0 <= R_F430S3 <= 999999
 0 <= R_F430S2 <= 999999
 0 <= R_F430S1 <= 999999
 0 <= R_TMDK1 <= 999999
-999999 <= R_CODHr <= 999999
 0 <= R_F430S5 <= 999999
 0 <= R_F430S4 <= 999999
 0 <= R_EXo_ni2_e <= 999999
-999999 <= R_DHORD7 <= 999999
-999999 <= R_N2tr <= 999999
-999999 <= R_GAPD <= 999999
 0 <= R_METAT <= 999999
 0 <= R_EXo_pnto__R_e <= 999999
-999999 <= R_RH3CD <= 999999
 0 <= R_RNDR4 <= 999999
 0 <= R_RNDR2 <= 999999
 0 <= R_RNDR3 <= 999999
 0 <= R_COabc <= 999999
 0 <= R_RNDR1 <= 999999
 0 <= R_NIT_n1p4 <= 999999
 0 <= R_PYK <= 999999
 0 <= R_MTCMMT <= 999999
 0 <= R_EXi_nac_e <= 999999
 0 <= R_ACGK <= 999999
-999999 <= R_PGAMT <= 999999
 0 <= R_ALDD1 <= 999999
-999999 <= R_HCO3E <= 999999
 0 <= R_RIBFLVt2 <= 999999
 0 <= R_ACGS <= 999999
 0 <= R_CHORM <= 999999
 0 <= R_EXi_mobd_e <= 999999
 0 <= R_RBFSb <= 999999
 0 <= R_RBFSa <= 999999
-999999 <= R_IOR <= 999999
 0 <= R_GOFUCR <= 999999
 0 <= R_PPS <= 999999
 0 <= R_CYRDAR <= 999999
 0 <= R_CYRDAS <= 999999
 0 <= R_NNDPR <= 999999
-999999 <= R_FUM <= 999999
      R_EXi_thm_e = 0
 0 <= R_ACONMT <= 999999
 0 <= R_ADCYRS <= 999999
 0 <= R_GLNS <= 999999
 0 <= R_PPA <= 999999
 0 <= R_PPCDC <= 999999
 0 <= R_PSLS <= 999999
-999999 <= R_IOR3 <= 999999
-999999 <= R_IOR2 <= 999999
-999999 <= R_CH3ACH <= 999999
 0 <= R_LTHRK <= 999999
 0 <= R_ALATRS <= 999999
 0 <= R_SECS <= 999999
-999999 <= R_ALCD20y <= 999999
 0 <= R_FOLR2 <= 999999
 0 <= R_ADSS <= 999999
 0 <= R_BACCL <= 999999
-999999 <= R_UNK_RBFDEGt <= 999999
-999999 <= R_DTMPK <= 999999
-999999 <= R_POR2 <= 999999
 0 <= R_CD2abc2 <= 999999
 0 <= R_CD2abc1 <= 999999
 0 <= R_NTP5 <= 999999
 0 <= R_NTP4 <= 999999
-999999 <= R_ACOTA <= 999999
 0 <= R_NTP6 <= 999999
 0 <= R_NTP1 <= 999999
 0 <= R_NTP3 <= 999999
 0 <= R_NTP2 <= 999999
 0 <= R_NTP9 <= 999999
 0 <= R_NTP8 <= 999999
 0 <= R_CPC6R <= 999999
 0 <= R_EXi_ac_e <= 999999
 0 <= R_MOHMT <= 999999
 0 <= R_ACLS <= 999999
 0 <= R_DNMPPA <= 999999
 0 <= R_EXo_ind3ac_e <= 999999
 0 <= R_TAPTPP <= 999999
 0 <= R_EXi_ile__L_e <= 999999
 0 <= R_GLCS2 <= 999999
 0 <= R_ARGTRS <= 999999
-999999 <= R_GLYt4 <= 999999
 0 <= R_TSR <= 999999
 0 <= R_ATPS1 <= 999999
 0 <= R_EXo_cl_e <= 999999
 0 <= R_THZPSN <= 999999
-999999 <= R_ILEt2r <= 999999
 0 <= R_EXo_dma_e <= 999999
 0 <= R_EXo_ribflv_e <= 999999
 0 <= R_ALKP <= 999999
 0 <= R_DAPDC <= 999999
 0 <= R_PPBNGS <= 999999
-999999 <= R_CYTK1 <= 999999
-999999 <= R_IND3ACt2r <= 999999
-999999 <= R_CYTK2 <= 999999
 0 <= R_CODH2r <= 999999
 0 <= R_EXo_tsul_e <= 999999
 0 <= R_AGAID <= 999999
 0 <= R_EXi_btn_e <= 999999
 0 <= R_GMAND <= 999999
-999999 <= R_TMPK <= 999999
 0 <= R_UPP3MT <= 999999
 0 <= R_HDR <= 999999
 0 <= R_GLUDxi <= 999999
 0 <= R_F4D <= 999999
 0 <= R_AKP1 <= 999999
 0 <= R_EXi_na1_e <= 999999
 0 <= R_CBL1HBIabc <= 999999
-999999 <= R_VOR <= 999999
-999999 <= R_GTHOr <= 999999
 0 <= R_EXi_ind3ac_e <= 999999
-999999 <= R_KARA2 <= 999999
 0 <= R_DHPM2 <= 999999
-999999 <= R_AIRC3 <= 999999
 0 <= R_AIRC2 <= 999999
 0 <= R_FRTT2 <= 999999
 0 <= R_HSPMS <= 999999
 0 <= R_CDPGS <= 999999
-999999 <= R_GALT <= 999999
 0 <= R_TDPDRR <= 999999
-999999 <= R_DBTS <= 999999
-999999 <= R_ALAt4 <= 999999
 0 <= R_DAPDH <= 999999
 0 <= R_TDPDRE <= 999999
 0 <= R_ATIH <= 999999
 0 <= R_ILETRS <= 999999
      R_EXi_cl_e = 0
 0 <= R_CDPGH <= 999999
 0 <= R_PIabc <= 999999
 0 <= R_LEUTRS <= 999999
 0 <= R_IGPDH <= 999999
-999999 <= R_GK1 <= 999999
 0 <= R_EXo_gcald_e <= 999999
 0 <= R_GLYTRS <= 999999
 0 <= R_PMDPHT <= 999999
 0 <= R_GLUTRS_Gln <= 999999
 0 <= R_DHQS2 <= 999999
 0 <= R_RBK <= 999999
 0 <= R_Cuabc <= 999999
 0 <= R_EXi_fol_e <= 999999
 0 <= R_EXo_glyb_e <= 999999
-999999 <= R_AICART <= 999999
 0 <= R_EXi_k_e <= 999999
 0 <= R_EXi_so3_e <= 999999
-999999 <= R_Clt <= 999999
 0 <= R_APPLP <= 999999
 0 <= R_EXi_cobalt2_e <= 999999
 0 <= R_GLCGSD <= 999999
 0 <= R_5HBZIDS <= 999999
 0 <= R_THRPD <= 999999
 0 <= R_HMPK1 <= 999999
 0 <= R_DHFS <= 999999
-999999 <= R_DHFR <= 999999
 0 <= R_FEDCabc <= 999999
-999999 <= R_LEUt2r <= 999999
-999999 <= R_ALCD2y <= 999999
 0 <= R_EXo_unknown_cbl1deg_e <= 999999
 0 <= R_NaKt_1 <= 999999
 0 <= R_VALTRS <= 999999
 0 <= R_METS <= 999999
 0 <= R_ACACT1r <= 999999
-999999 <= R_MAN1PT2 <= 999999
-999999 <= R_UMPK_copy1 <= 999999
 0 <= R_BTNt2i <= 999999
 0 <= R_H4MPTS14 <= 999999
-999999 <= R_URIDK2r <= 999999
 0 <= R_PGLYCP <= 999999
 0 <= R_H4MPTS15 <= 999999
-999999 <= R_ICDHyr <= 999999
 0 <= R_AKSDC <= 999999
 0 <= R_EXo_etha_e <= 999999
 0 <= R_EXi_cit_e <= 999999
 0 <= R_GUAD <= 999999
-999999 <= R_BTS3r <= 999999
-999999 <= R_NAt3_1 <= 999999
-999999 <= R_NH4t <= 999999
 0 <= R_EXi_pnto__R_e <= 999999
 0 <= R_UAPGR <= 999999
 0 <= R_CYRDAAT <= 999999
 0 <= R_METTRS <= 999999
 0 <= R_PPND <= 999999
 0 <= R_EXo_ac_e <= 999999
 0 <= R_TMPPP <= 999999
 0 <= R_ANPRT <= 999999
-999999 <= R_UAG2E <= 999999
 0 <= R_NIabc <= 999999
 0 <= R_MCMMT <= 999999
 0 <= R_RZ5PP2 <= 999999
-999999 <= R_ATPS4r <= 999999
 0 <= R_SULabc <= 999999
 0 <= R_EXi_meoh_e <= 999999
-999999 <= R_GCALDt <= 999999
 0 <= R_EXo_k_e <= 999999
 0 <= R_UAGCVT <= 999999
 0 <= R_EXo_glyc_e <= 999999
 0 <= R_IPDDI3y <= 999999
 0 <= R_IPDDI3x <= 999999
 0 <= R_CBIabc <= 999999
 0 <= R_XPPT <= 999999
-999999 <= R_PMANM <= 999999
-999999 <= R_ALACt2r <= 999999
-999999 <= R_MTHFC <= 999999
 0 <= R_PGCD <= 999999
 0 <= R_COCHL <= 999999
 0 <= R_OMCDC <= 999999
 0 <= R_ACYP_2 <= 999999
      R_EXi_h2_e = 0
 0 <= R_EXi_glcn_e <= 999999
 0 <= R_H2MPTR <= 999999
 0 <= R_EXo_co_e <= 999999
 0 <= R_EXi_pro__L_e <= 999999
-999999 <= R_PGMT <= 999999
 0 <= R_EXi_glu__L_e <= 999999
 0 <= R_PMPK <= 999999
      R_EXi_cbl1_e = 0
 0 <= R_PYRt2 <= 999999
 0 <= R_DPR <= 999999
-999999 <= R_PPAKr <= 999999
-999999 <= R_LEUTA <= 999999
 0 <= R_MG2abc <= 999999
 0 <= R_EXo_ile__L_e <= 999999
-999999 <= R_VALt2r <= 999999
 0 <= R_EXo_h2o_e <= 999999
-999999 <= R_ACt2r <= 999999
 0 <= R_EXi_etha_e <= 999999
 0 <= R_TMAMT <= 999999
 0 <= R_EXi_fe3_e <= 999999
-999999 <= R_ME1 <= 999999
 0 <= R_G1PTT <= 999999
 0 <= R_EXi_dms_e <= 999999
 0 <= R_EXo_pac_e <= 999999
-999999 <= R_IPPMIa <= 999999
-999999 <= R_IPPMIb <= 999999
 0 <= R_GHMT2r <= 999999
 0 <= R_DKMPASPL <= 999999
 0 <= R_CDGGGS <= 999999
-999999 <= R_HSERTA <= 999999
 0 <= R_CPC6MT <= 999999
-999999 <= R_TKT1 <= 999999
-999999 <= R_TKT2 <= 999999
 0 <= R_MMAMT <= 999999
 0 <= R_EXo_cit_e <= 999999
 0 <= R_DHAD1 <= 999999
-999999 <= R_PNTOt2 <= 999999
 0 <= R_EXi_mma_e <= 999999
 0 <= R_QULNS <= 999999
 0 <= R_FE3abc <= 999999
 0 <= R_DMHDRFS <= 999999
-999999 <= R_ADKd <= 999999
-999999 <= R_MTSPCMMT <= 999999
 0 <= R_OMPDC <= 999999
      R_EXi_unknown_cbl1deg_e = 0
 0 <= R_COBS <= 999999
 0 <= R_EXo_fe2_e <= 999999
-999999 <= R_MTSPC <= 999999
 0 <= R_CDGGGSAT2 <= 999999
 0 <= R_HXPRT <= 999999
 0 <= R_MI1PS <= 999999
-999999 <= R_CYSS <= 999999
 0 <= R_MI1PP <= 999999
-999999 <= R_SERAT <= 999999
 0 <= R_HATSH <= 999999
-999999 <= R_F4H2O <= 999999
 0 <= R_CS <= 999999
 0 <= R_CBPS <= 999999
 0 <= R_CPC2MT <= 999999
 0 <= R_FBPO <= 999999
      R_EXi_so4_e = 0
-999999 <= R_ADK1 <= 999999
-999999 <= R_ORPT <= 999999
-999999 <= R_ADK3 <= 999999
-999999 <= R_ADK2 <= 999999
-999999 <= R_ADK4 <= 999999
-999999 <= R_AHGDx <= 999999
 0 <= R_MTHFD2i <= 999999
-999999 <= R_ECHH_10 <= 999999
 0 <= R_EXo_zn2_e <= 999999
 0 <= R_ADD <= 999999
-999999 <= R_LYSt3r <= 999999
-999999 <= R_RPE <= 999999
 0 <= R_MNabc <= 999999
 0 <= R_ADCPS2 <= 999999
-999999 <= R_RPI <= 999999
 0 <= R_DHPCPH <= 999999
 0 <= R_SHCHF <= 999999
 0 <= R_IPPS <= 999999
 0 <= R_DHDPS <= 999999
-999999 <= R_COt <= 999999
 0 <= R_TAPTPPM <= 999999
 0 <= R_DPGase <= 999999
 0 <= R_ASNS1 <= 999999
-999999 <= R_OPAH <= 999999
 0 <= R_F4NH <= 999999
 0 <= R_AKSCAL <= 999999
-999999 <= R_ALATA_L <= 999999
 0 <= R_7MHTS <= 999999
      R_EXi_cbl1hbi_e = 0
 0 <= R_EXo_leu__L_e <= 999999
-999999 <= R_SLDx <= 999999
-999999 <= R_SLDy <= 999999
 0 <= R_NACt <= 999999
 0 <= R_2PGK <= 999999
 0 <= R_MDHi2 <= 999999
      R_EXo_h2s_e = 0
 0 <= R_EXo_ch4_e <= 999999
      R_EXi_co2_e = 0
 0 <= R_FBP <= 999999
 0 <= R_MOBDabc <= 999999
 0 <= R_TDPGDH <= 999999
 0 <= R_DCTPD <= 999999
-999999 <= R_ACKr <= 999999
 0 <= R_CYSTRS <= 999999
 0 <= R_GLYBabc <= 999999
 0 <= R_ADPT <= 999999
 0 <= R_DMATT2 <= 999999
-999999 <= R_FBA <= 999999
 0 <= R_EXi_n2_e <= 999999
 0 <= R_DROPPRx <= 999999
 0 <= R_Cut1 <= 999999
 0 <= R_EXi_glyc_e <= 999999
 0 <= R_HISTP <= 999999
 0 <= R_EXi_lys__L_e <= 999999
 0 <= R_EXi_dma_e <= 999999
 0 <= R_IG3PS <= 999999
 0 <= R_DMATT <= 999999
 0 <= R_ADNK1 <= 999999
 0 <= R_Kt2r_copy1 <= 999999
-999999 <= R_Kt2r_copy2 <= 999999
 0 <= R_THMabc <= 999999
 0 <= R_G1SAT <= 999999
 0 <= R_TIHCD <= 999999
 0 <= R_HISTD <= 999999
 0 <= R_ATSH <= 999999
 0 <= R_TYRTRS <= 999999
 0 <= R_CDGGIPT2 <= 999999
 0 <= R_CDGGGS2 <= 999999
 0 <= R_ACYP <= 999999
 0 <= R_ATPPRT <= 999999
 0 <= R_P5CRx <= 999999
 0 <= R_EXo_mma_e <= 999999
 0 <= R_LYSTRS <= 999999
-999999 <= R_HPACt2r <= 999999
 0 <= R_H4MPTGL_atp <= 999999
-999999 <= R_UDPG4E <= 999999
 0 <= R_MTAP <= 999999
 0 <= R_GK_adp <= 999999
-999999 <= R_THRAi <= 999999
 0 <= R_OXADC <= 999999
-999999 <= R_SO3t <= 999999
-999999 <= R_HSDy <= 999999
 0 <= R_RBPC <= 999999
-999999 <= R_CPC8MM <= 999999
      R_EXi_ca2_e = 0
 0 <= R_RBPI <= 999999
-999999 <= R_CBIAT <= 999999
 0 <= R_EXo_nh4_e <= 999999
 0 <= R_GTPCIII <= 999999
-999999 <= R_UAG2EMA <= 999999
-999999 <= R_YUMPS <= 999999
-999999 <= R_AOXSr <= 999999
 0 <= R_CBLD <= 999999
-999999 <= R_GLCNt2r <= 999999
 0 <= R_G1PACT <= 999999
-999999 <= R_DADK <= 999999
 0 <= R_EXo_mobd_e <= 999999
 0 <= R_EXo_cobalt2_e <= 999999
 0 <= R_EXi_glyald_e <= 999999
 0 <= R_EXo_so4_e <= 999999
 0 <= R_EXi_fe2_e <= 999999
 0 <= R_EXi_nh4_e <= 999999
 0 <= R_GLUPRT <= 999999
 0 <= R_CPC4MT <= 999999
 0 <= R_2PLS <= 999999
-999999 <= R_COBALTt5 <= 999999
 0 <= R_UAGDP <= 999999
 0 <= R_ACHBS <= 999999
 0 <= R_EXo_lys__L_e <= 999999
-999999 <= R_DMAt <= 999999
-999999 <= R_RH2CD <= 999999
-999999 <= R_UNK_CBL1DEGt <= 999999
-999999 <= R_MAN6PI <= 999999
 0 <= R_HPPK2 <= 999999
 0 <= R_PSP_L <= 999999
 0 <= R_ADSL2r <= 999999
-999999 <= R_ARGSL <= 999999
-999999 <= R_LYSAM <= 999999
-999999 <= R_ILETA <= 999999
 0 <= R_HATIH <= 999999
 0 <= R_EXo_ca2_e <= 999999
-999999 <= R_ORNTAC <= 999999
-999999 <= R_MPML <= 999999
 0 <= R_PTPATi <= 999999
-999999 <= R_MEOHt2 <= 999999
 0 <= R_ADOCBLS2 <= 999999
 0 <= R_SPDC <= 999999
 0 <= R_CHORS <= 999999
      R_EXi_tma_e = 0
 0 <= R_HISTRS <= 999999
 0 <= R_EXo_pi_e <= 999999
 0 <= R_THRTRS <= 999999
 0 <= R_CSND <= 999999
-999999 <= R_NBAHH_ir <= 999999
 0 <= R_EXo_fe3_e <= 999999
 0 <= R_CF3Sg <= 999999
 0 <= R_PRAIS <= 999999
 0 <= R_PRAIi <= 999999
 0 <= R_CF3Sa <= 999999
 0 <= R_PSLP <= 999999
 0 <= R_EXo_val__L_e <= 999999
-999999 <= R_DCMPDA <= 999999
      R_EXi_cbi_e = 0
 0 <= R_EXo_glcn_e <= 999999
 0 <= R_SHCHCC <= 999999
 0 <= R_PPDK <= 999999
 0 <= R_CDGGGSAT <= 999999
 0 <= R_ACLDC <= 999999
 0 <= R_CPC5MT <= 999999
-999999 <= R_GARFT <= 999999
 0 <= R_DIPS <= 999999
 0 <= R_PHETRS <= 999999
 0 <= R_PRATPP <= 999999
 0 <= R_EXo_h_e <= 999999
 0 <= R_EXo_cbi_e <= 999999
 0 <= R_ASPK <= 999999
 0 <= R_DGGGPS <= 999999
 0 <= R_AHSERL2 <= 999999
 0 <= bin_R_EXi_4abz_e <= 1
 0 <= bin_R_EXi_4hphac_e <= 1
 0 <= bin_R_EXi_ac_e <= 1
 0 <= bin_R_EXi_actn__R_e <= 1
 0 <= bin_R_EXi_ala__L_e <= 1
 0 <= bin_R_EXi_alac__S_e <= 1
 0 <= bin_R_EXi_btn_e <= 1
 0 <= bin_R_EXi_ca2_e <= 1
 0 <= bin_R_EXi_cbi_e <= 1
 0 <= bin_R_EXi_cbl1_e <= 1
 0 <= bin_R_EXi_cbl1hbi_e <= 1
 0 <= bin_R_EXi_cd2_e <= 1
 0 <= bin_R_EXi_ch4_e <= 1
 0 <= bin_R_EXi_ch4s_e <= 1
 0 <= bin_R_EXi_cit_e <= 1
 0 <= bin_R_EXi_cl_e <= 1
 0 <= bin_R_EXi_cobalt2_e <= 1
 0 <= bin_R_EXi_cu2_e <= 1
 0 <= bin_R_EXi_cys__L_e <= 1
 0 <= bin_R_EXi_dma_e <= 1
 0 <= bin_R_EXi_dms_e <= 1
 0 <= bin_R_EXi_etha_e <= 1
 0 <= bin_R_EXi_fe2_e <= 1
 0 <= bin_R_EXi_fe3_e <= 1
 0 <= bin_R_EXi_fol_e <= 1
 0 <= bin_R_EXi_gcald_e <= 1
 0 <= bin_R_EXi_glcn_e <= 1
 0 <= bin_R_EXi_glu__L_e <= 1
 0 <= bin_R_EXi_gly_e <= 1
 0 <= bin_R_EXi_glyald_e <= 1
 0 <= bin_R_EXi_glyb_e <= 1
 0 <= bin_R_EXi_glyc_e <= 1
 0 <= bin_R_EXi_h2s_e <= 1
 0 <= bin_R_EXi_ile__L_e <= 1
 0 <= bin_R_EXi_ind3ac_e <= 1
 0 <= bin_R_EXi_k_e <= 1
 0 <= bin_R_EXi_leu__L_e <= 1
 0 <= bin_R_EXi_lys__L_e <= 1
 0 <= bin_R_EXi_meoh_e <= 1
 0 <= bin_R_EXi_mg2_e <= 1
 0 <= bin_R_EXi_mma_e <= 1
 0 <= bin_R_EXi_mn2_e <= 1
 0 <= bin_R_EXi_mobd_e <= 1
 0 <= bin_R_EXi_n2_e <= 1
 0 <= bin_R_EXi_na1_e <= 1
 0 <= bin_R_EXi_nac_e <= 1
 0 <= bin_R_EXi_nh4_e <= 1
 0 <= bin_R_EXi_ni2_e <= 1
 0 <= bin_R_EXi_pac_e <= 1
 0 <= bin_R_EXi_pi_e <= 1
 0 <= bin_R_EXi_pnto__R_e <= 1
 0 <= bin_R_EXi_pro__L_e <= 1
 0 <= bin_R_EXi_pyr_e <= 1
 0 <= bin_R_EXi_ribflv_e <= 1
 0 <= bin_R_EXi_s_e <= 1
 0 <= bin_R_EXi_so3_e <= 1
 0 <= bin_R_EXi_so4_e <= 1
 0 <= bin_R_EXi_urea_e <= 1
 0 <= bin_R_EXi_val__L_e <= 1
Binaries
 bin_R_EXi_4abz_e  bin_R_EXi_4hphac_e  bin_R_EXi_ac_e  bin_R_EXi_actn__R_e 
 bin_R_EXi_ala__L_e  bin_R_EXi_alac__S_e  bin_R_EXi_btn_e  bin_R_EXi_ca2_e 
 bin_R_EXi_cbi_e  bin_R_EXi_cbl1_e  bin_R_EXi_cbl1hbi_e  bin_R_EXi_cd2_e 
 bin_R_EXi_ch4_e  bin_R_EXi_ch4s_e  bin_R_EXi_cit_e  bin_R_EXi_cl_e 
 bin_R_EXi_cobalt2_e  bin_R_EXi_cu2_e  bin_R_EXi_cys__L_e  bin_R_EXi_dma_e 
 bin_R_EXi_dms_e  bin_R_EXi_etha_e  bin_R_EXi_fe2_e  bin_R_EXi_fe3_e 
 bin_R_EXi_fol_e  bin_R_EXi_gcald_e  bin_R_EXi_glcn_e  bin_R_EXi_glu__L_e 
 bin_R_EXi_gly_e  bin_R_EXi_glyald_e  bin_R_EXi_glyb_e  bin_R_EXi_glyc_e 
 bin_R_EXi_h2s_e  bin_R_EXi_ile__L_e  bin_R_EXi_ind3ac_e  bin_R_EXi_k_e 
 bin_R_EXi_leu__L_e  bin_R_EXi_lys__L_e  bin_R_EXi_meoh_e  bin_R_EXi_mg2_e 
 bin_R_EXi_mma_e  bin_R_EXi_mn2_e  bin_R_EXi_mobd_e  bin_R_EXi_n2_e 
 bin_R_EXi_na1_e  bin_R_EXi_nac_e  bin_R_EXi_nh4_e  bin_R_EXi_ni2_e 
 bin_R_EXi_pac_e  bin_R_EXi_pi_e  bin_R_EXi_pnto__R_e  bin_R_EXi_pro__L_e 
 bin_R_EXi_pyr_e  bin_R_EXi_ribflv_e  bin_R_EXi_s_e  bin_R_EXi_so3_e 
 bin_R_EXi_so4_e  bin_R_EXi_urea_e  bin_R_EXi_val__L_e 
End
