\ENCODING=ISO-8859-1
\Problem name: 

Maximize
 pd_R09: R09 + 0 pADP + 0 pATP + 0 pM1 + 0 pM10 + 0 pM11 + 0 pM12 + 0 pM13
         + 0 pM2 + 0 pM3 + 0 pM4 + 0 pM5 + 0 pM8 + 0 pM9 + 0 bin_G25
         + 0 bin_G16 + 0 bin_G25 + 0 bin_G37 + 0 bin_G25 + 0 bin_G26
         + 0 bin_G21 + 0 bin_G02 + 0 bin_G20 + 0 bin_G21 + 0 bin_G23
         + 0 bin_G03 + 0 bin_G08 + 0 bin_G08
Subject To
 pADP:                      R01 - R03 - R13  = 0
 pATP:                      - R01 + R03 + R13  = 0
 pM1:                       R01 - R02 - R10  = 0
 pM10:                      R12  = 0
 pM11:                      R13 - R14  = 0
 pM12:                      0.7 R14 - R16  = 0
 pM13:                      0.3 R14 - R15  = 0
 pM2:                       - R03 + R02 + R11  = 0
 pM3:                       R03 - R04 - R06 - R18  = 0
 pM4:                       R04 - R05  = 0
 pM5:                       - R17 + R16 + R05 + 0.5 R08  = 0
 pM6:                       R06 - R07  = 0
 pM7:                       - R08 + R07  = 0
 pM8:                       2 R18 + 0.5 R08 - R09  = 0
 pM9:                       - R13 + R10 - R12 - R11  = 0
 dR07:                      - pM6 + pM7 + u17 >= 0
 objLinkC:                  R17 + 10 l1 - 1000 l5 - 1000 l6 - 1000 l11
                            - 1000 l12 - 1000 l15 - 1000 u1 - 100 u2 - 1000 u3
                            - 1000 u4 - 1000 u5 - 1000 u6 - 1000 u7 - 1000 u8
                            - 1000 u9 - 1000 u10 - 1000 u11 - 1000 u12
                            - 1000 u13 - 1000 u14 - 1000 u15 - 1000 u16
                            - 1000 u17 - 1000 u18  = 0
 MinActG:                   bin_G01 + bin_G19 + bin_G16 + bin_G25 + bin_G37
                            + bin_G23 + bin_G03 + bin_G02 + bin_G20 + bin_G21
                            + bin_G05 + bin_G26 + bin_G04 + bin_G24 + bin_G10
                            + bin_G31 + bin_G32 + bin_G06 + bin_G27 + bin_G28
                            + bin_G09 + bin_G30 + bin_G08 + bin_G11 + bin_G22
                            + bin_G14 + bin_G34 + bin_G35 + bin_G18 + bin_G15
                            + bin_G36 + bin_G12 + bin_G33 + bin_G17 + bin_G38
                            + bin_G07 >= 34
 ON_bin_R17:                bin_R17  = 1
 bin_E1_R12:                bin_R12 - bin_E1_R12 >= 0
 bin_R12_E:                 bin_R12 - bin_E1_R12 <= 0
 bin_E1_R15:                bin_R15 - bin_E1_R15 >= 0
 bin_R15_E:                 bin_R15 - bin_E1_R15 <= 0
 bin_E1_R11:                bin_R11 - bin_E1_R11 >= 0
 bin_E2_R11:                bin_R11 - bin_E2_R11 >= 0
 bin_R11_E:                 bin_R11 - bin_E1_R11 - bin_E2_R11 <= 0
 bin_E1_R05:                bin_R05 - bin_E1_R05 >= 0
 bin_E2_R05:                bin_R05 - bin_E2_R05 >= 0
 bin_R05_E:                 bin_R05 - bin_E1_R05 - bin_E2_R05 <= 0
 bin_E1_R14:                bin_R14 - bin_E1_R14 >= 0
 bin_R14_E:                 bin_R14 - bin_E1_R14 <= 0
 bin_E1_R06:                bin_R06 - bin_E1_R06 >= 0
 bin_R06_E:                 bin_R06 - bin_E1_R06 <= 0
 bin_E1_R16:                bin_R16 - bin_E1_R16 >= 0
 bin_R16_E:                 bin_R16 - bin_E1_R16 <= 0
 bin_E1_R04:                bin_R04 - bin_E1_R04 >= 0
 bin_R04_E:                 bin_R04 - bin_E1_R04 <= 0
 bin_E1_R03:                bin_R03 - bin_E1_R03 >= 0
 bin_E2_R03:                bin_R03 - bin_E2_R03 >= 0
 bin_R03_E:                 bin_R03 - bin_E1_R03 - bin_E2_R03 <= 0
 bin_E1_R02:                bin_R02 - bin_E1_R02 >= 0
 bin_R02_E:                 bin_R02 - bin_E1_R02 <= 0
 bin_E1_R01:                bin_R01 - bin_E1_R01 >= 0
 bin_R01_E:                 bin_R01 - bin_E1_R01 <= 0
 bin_E1_R13:                bin_R13 - bin_E1_R13 >= 0
 bin_R13_E:                 bin_R13 - bin_E1_R13 <= 0
 bin_E1_R10:                bin_R10 - bin_E1_R10 >= 0
 bin_R10_E:                 bin_R10 - bin_E1_R10 <= 0
 bin_E1_R18:                bin_R18 - bin_E1_R18 >= 0
 bin_R18_E:                 bin_R18 - bin_E1_R18 <= 0
 bin_E1_R17:                bin_R17 - bin_E1_R17 >= 0
 bin_R17_E:                 bin_R17 - bin_E1_R17 <= 0
 bin_E1_R09:                bin_R09 - bin_E1_R09 >= 0
 bin_R09_E:                 bin_R09 - bin_E1_R09 <= 0
 bin_E1_R08:                bin_R08 - bin_E1_R08 >= 0
 bin_R08_E:                 bin_R08 - bin_E1_R08 <= 0
 bin_G09_bin_E1_R09:        bin_E1_R09 - bin_G09 <= 0
 bin_G30_bin_E1_R09:        bin_E1_R09 - bin_G30 <= 0
 bin_G30_E:                 bin_E1_R09 - bin_G09 - bin_G30 >= -1
 bin_G08_bin_E1_R08:        bin_E1_R08 - bin_G08 <= 0
 bin_G25_bin_E1_R08:        - bin_G25 + bin_E1_R08 <= 0
 bin_G26_bin_E1_R08:        - bin_G26 + bin_E1_R08 <= 0
 bin_G26_E:                 - bin_G25 - bin_G26 + bin_E1_R08 - bin_G08 >= -2
 bin_G18_bin_E1_R18:        bin_E1_R18 - bin_G18 <= 0
 bin_G08_or_G07_bin_E1_R18: bin_E1_R18 - bin_G08_or_G07 <= 0
 bin_G08_or_G07_E:          bin_E1_R18 - bin_G18 - bin_G08_or_G07 >= -1
 bin_G16_bin_E2_R05:        - bin_G16 + bin_E2_R05 <= 0
 bin_G25_bin_E2_R05:        - bin_G25 + bin_E2_R05 <= 0
 bin_G37_bin_E2_R05:        - bin_G37 + bin_E2_R05 <= 0
 bin_G37_E:                 - bin_G16 - bin_G25 - bin_G37 + bin_E2_R05 >= -2
 bin_G02_bin_E2_R11:        - bin_G02 + bin_E2_R11 <= 0
 bin_G20_bin_E2_R11:        - bin_G20 + bin_E2_R11 <= 0
 bin_G21_bin_E2_R11:        - bin_G21 + bin_E2_R11 <= 0
 bin_G21_E:                 - bin_G02 - bin_G20 - bin_G21 + bin_E2_R11 >= -2
 bin_G03_bin_E2_R03:        bin_E2_R03 - bin_G03 <= 0
 bin_G03_E:                 bin_E2_R03 - bin_G03 >= 0
 bin_G01_bin_E1_R01:        bin_E1_R01 - bin_G01 <= 0
 bin_G19_bin_E1_R01:        bin_E1_R01 - bin_G19 <= 0
 bin_G19_E:                 bin_E1_R01 - bin_G01 - bin_G19 >= -1
 bin_G16_bin_E1_R16:        bin_E1_R16 - bin_G16 <= 0
 bin_G25_bin_E1_R16:        bin_E1_R16 - bin_G25 <= 0
 bin_G37_bin_E1_R16:        bin_E1_R16 - bin_G37 <= 0
 bin_G37_E:                 bin_E1_R16 - bin_G16 - bin_G25 - bin_G37 >= -2
 bin_G23_bin_E1_R03:        bin_E1_R03 - bin_G23 <= 0
 bin_G23_E:                 bin_E1_R03 - bin_G23 >= 0
 bin_G02_bin_E1_R02:        bin_E1_R02 - bin_G02 <= 0
 bin_G20_bin_E1_R02:        bin_E1_R02 - bin_G20 <= 0
 bin_G21_bin_E1_R02:        bin_E1_R02 - bin_G21 <= 0
 bin_G21_E:                 bin_E1_R02 - bin_G02 - bin_G20 - bin_G21 >= -2
 bin_G05_bin_E1_R05:        bin_E1_R05 - bin_G05 <= 0
 bin_G25_bin_E1_R05:        - bin_G25 + bin_E1_R05 <= 0
 bin_G26_bin_E1_R05:        bin_E1_R05 - bin_G26 <= 0
 bin_G26_E:                 - bin_G25 + bin_E1_R05 - bin_G05 - bin_G26 >= -2
 bin_G04_bin_E1_R04:        bin_E1_R04 - bin_G04 <= 0
 bin_G24_bin_E1_R04:        bin_E1_R04 - bin_G24 <= 0
 bin_G24_E:                 bin_E1_R04 - bin_G04 - bin_G24 >= -1
 bin_G14_bin_E1_R14:        bin_E1_R14 - bin_G14 <= 0
 bin_G34_bin_E1_R14:        bin_E1_R14 - bin_G34 <= 0
 bin_G35_bin_E1_R14:        bin_E1_R14 - bin_G35 <= 0
 bin_G35_E:                 bin_E1_R14 - bin_G14 - bin_G34 - bin_G35 >= -2
 bin_G06_bin_E1_R06:        bin_E1_R06 - bin_G06 <= 0
 bin_G27_bin_E1_R06:        bin_E1_R06 - bin_G27 <= 0
 bin_G28_bin_E1_R06:        bin_E1_R06 - bin_G28 <= 0
 bin_G28_E:                 bin_E1_R06 - bin_G06 - bin_G27 - bin_G28 >= -2
 bin_G12_bin_E1_R12:        bin_E1_R12 - bin_G12 <= 0
 bin_G33_bin_E1_R12:        bin_E1_R12 - bin_G33 <= 0
 bin_G33_E:                 bin_E1_R12 - bin_G12 - bin_G33 >= -1
 bin_G03_or_G08_bin_E1_R13: bin_E1_R13 - bin_G03_or_G08 <= 0
 bin_G23_bin_E1_R13:        - bin_G23 + bin_E1_R13 <= 0
 bin_G23_E:                 - bin_G23 + bin_E1_R13 - bin_G03_or_G08 >= -1
 bin_G17_bin_E1_R17:        bin_E1_R17 - bin_G17 <= 0
 bin_G38_bin_E1_R17:        bin_E1_R17 - bin_G38 <= 0
 bin_G38_E:                 bin_E1_R17 - bin_G17 - bin_G38 >= -1
 bin_G15_bin_E1_R15:        bin_E1_R15 - bin_G15 <= 0
 bin_G36_bin_E1_R15:        bin_E1_R15 - bin_G36 <= 0
 bin_G36_E:                 bin_E1_R15 - bin_G15 - bin_G36 >= -1
 bin_G10_bin_E1_R10:        bin_E1_R10 - bin_G10 <= 0
 bin_G31_bin_E1_R10:        bin_E1_R10 - bin_G31 <= 0
 bin_G32_bin_E1_R10:        bin_E1_R10 - bin_G32 <= 0
 bin_G32_E:                 bin_E1_R10 - bin_G10 - bin_G31 - bin_G32 >= -2
 bin_G11_bin_E1_R11:        bin_E1_R11 - bin_G11 <= 0
 bin_G22_bin_E1_R11:        bin_E1_R11 - bin_G22 <= 0
 bin_G21_bin_E1_R11:        - bin_G21 + bin_E1_R11 <= 0
 bin_G21_E:                 - bin_G21 + bin_E1_R11 - bin_G11 - bin_G22 >= -2
 G03_or_G08_G03:            - bin_G03 + bin_G03_or_G08 >= 0
 G03_or_G08_G08:            - bin_G08 + bin_G03_or_G08 >= 0
 G03_or_G08:                - bin_G03 - bin_G08 + bin_G03_or_G08 <= 0
 G08_or_G07_G08:            - bin_G08 + bin_G08_or_G07 >= 0
 G08_or_G07_G07:            bin_G08_or_G07 - bin_G07 >= 0
 G08_or_G07:                - bin_G08 + bin_G08_or_G07 - bin_G07 <= 0
 ind_R01:  bin_R01 = 0 -> R01  = 0
 ind_R16:  bin_R16 = 0 -> R16  = 0
 ind_R03:  bin_R03 = 0 -> R03  = 0
 ind_R02:  bin_R02 = 0 -> R02  = 0
 ind_R05:  bin_R05 = 0 -> R05  = 0
 ind_R04:  bin_R04 = 0 -> R04  = 0
 ind_R10:  bin_R10 = 0 -> R10  = 0
 ind_R06:  bin_R06 = 0 -> R06  = 0
 ind_R09:  bin_R09 = 0 -> R09  = 0
 ind_R08:  bin_R08 = 0 -> R08  = 0
 ind_R11:  bin_R11 = 0 -> R11  = 0
 ind_R14:  bin_R14 = 0 -> R14  = 0
 ind_R18:  bin_R18 = 0 -> R18  = 0
 ind_R15:  bin_R15 = 0 -> R15  = 0
 ind_R12:  bin_R12 = 0 -> R12  = 0
 ind_R17:  bin_R17 = 0 -> R17  = 0
 ind_R13:  bin_R13 = 0 -> R13  = 0
 ind_dR17: bin_R17 = 1 -> - pM5 - l1 + u1  = 1
 ind_dR01: bin_R01 = 1 -> pADP - pATP + pM1 + u2 >= 0
 ind_dR03: bin_R03 = 1 -> - pADP + pATP - pM2 + pM3 + u3 >= 0
 ind_dR13: bin_R13 = 1 -> - pADP + pATP + pM11 - pM9 + u4 >= 0
 ind_dR02: bin_R02 = 1 -> - pM1 + pM2 - l5 + u5  = 0
 ind_dR10: bin_R10 = 1 -> - pM1 + pM9 - l6 + u6  = 0
 ind_dR12: bin_R12 = 1 -> pM10 - pM9 + u7 >= 0
 ind_dR14: bin_R14 = 1 -> - pM11 + 0.7 pM12 + 0.3 pM13 + u8 >= 0
 ind_dR16: bin_R16 = 1 -> - pM12 + pM5 + u9 >= 0
 ind_dR15: bin_R15 = 1 -> - pM13 + u10 >= 0
 ind_dR11: bin_R11 = 1 -> pM2 - pM9 - l11 + u11  = 0
 ind_dR04: bin_R04 = 1 -> - pM3 + pM4 - l12 + u12  = 0
 ind_dR06: bin_R06 = 1 -> - pM3 + pM6 + u13 >= 0
 ind_dR18: bin_R18 = 1 -> - pM3 + 2 pM8 + u14 >= 0
 ind_dR05: bin_R05 = 1 -> - pM4 + pM5 - l15 + u15  = 0
 ind_dR08: bin_R08 = 1 -> 0.5 pM5 - pM7 + 0.5 pM8 + u16 >= 0
 ind_dR09: bin_R09 = 1 -> - pM8 + u18 >= 0
Bounds
 10 <= R17 <= 1000
 0 <= R01 <= 100
 0 <= R03 <= 1000
 0 <= R13 <= 1000
-1000 <= R02 <= 1000
-1000 <= R10 <= 1000
 0 <= R12 <= 1000
 0 <= R14 <= 1000
 0 <= R16 <= 1000
 0 <= R15 <= 1000
-1000 <= R11 <= 1000
-1000 <= R04 <= 1000
 0 <= R06 <= 1000
 0 <= R18 <= 1000
-1000 <= R05 <= 1000
 0 <= R08 <= 1000
 0 <= R07 <= 1000
 0 <= R09 <= 1000
      pADP Free
      pATP Free
      pM1 Free
      pM10 Free
      pM11 Free
      pM12 Free
      pM13 Free
      pM2 Free
      pM3 Free
      pM4 Free
      pM5 Free
      pM6 Free
      pM7 Free
      pM8 Free
      pM9 Free
 0 <= bin_R01 <= 1
 0 <= bin_E1_R01 <= 1
 0 <= bin_G01 <= 1
 0 <= bin_G19 <= 1
 0 <= bin_R16 <= 1
 0 <= bin_E1_R16 <= 1
 0 <= bin_G16 <= 1
 0 <= bin_G25 <= 1
 0 <= bin_G37 <= 1
 0 <= bin_R03 <= 1
 0 <= bin_E1_R03 <= 1
 0 <= bin_G23 <= 1
 0 <= bin_E2_R03 <= 1
 0 <= bin_G03 <= 1
 0 <= bin_R02 <= 1
 0 <= bin_E1_R02 <= 1
 0 <= bin_G02 <= 1
 0 <= bin_G20 <= 1
 0 <= bin_G21 <= 1
 0 <= bin_R05 <= 1
 0 <= bin_E1_R05 <= 1
 0 <= bin_G05 <= 1
 0 <= bin_G25 <= 1
 0 <= bin_G26 <= 1
 0 <= bin_E2_R05 <= 1
 0 <= bin_G16 <= 1
 0 <= bin_G25 <= 1
 0 <= bin_G37 <= 1
 0 <= bin_R04 <= 1
 0 <= bin_E1_R04 <= 1
 0 <= bin_G04 <= 1
 0 <= bin_G24 <= 1
 0 <= bin_R10 <= 1
 0 <= bin_E1_R10 <= 1
 0 <= bin_G10 <= 1
 0 <= bin_G31 <= 1
 0 <= bin_G32 <= 1
 0 <= bin_R06 <= 1
 0 <= bin_E1_R06 <= 1
 0 <= bin_G06 <= 1
 0 <= bin_G27 <= 1
 0 <= bin_G28 <= 1
 0 <= bin_R09 <= 1
 0 <= bin_E1_R09 <= 1
 0 <= bin_G09 <= 1
 0 <= bin_G30 <= 1
 0 <= bin_R08 <= 1
 0 <= bin_E1_R08 <= 1
 0 <= bin_G08 <= 1
 0 <= bin_G25 <= 1
 0 <= bin_G26 <= 1
 0 <= bin_R11 <= 1
 0 <= bin_E1_R11 <= 1
 0 <= bin_G11 <= 1
 0 <= bin_G22 <= 1
 0 <= bin_G21 <= 1
 0 <= bin_E2_R11 <= 1
 0 <= bin_G02 <= 1
 0 <= bin_G20 <= 1
 0 <= bin_G21 <= 1
 0 <= bin_R14 <= 1
 0 <= bin_E1_R14 <= 1
 0 <= bin_G14 <= 1
 0 <= bin_G34 <= 1
 0 <= bin_G35 <= 1
 0 <= bin_R18 <= 1
 0 <= bin_E1_R18 <= 1
 0 <= bin_G18 <= 1
 0 <= bin_G08_or_G07 <= 1
 0 <= bin_R15 <= 1
 0 <= bin_E1_R15 <= 1
 0 <= bin_G15 <= 1
 0 <= bin_G36 <= 1
 0 <= bin_R12 <= 1
 0 <= bin_E1_R12 <= 1
 0 <= bin_G12 <= 1
 0 <= bin_G33 <= 1
 0 <= bin_R17 <= 1
 0 <= bin_E1_R17 <= 1
 0 <= bin_G17 <= 1
 0 <= bin_G38 <= 1
 0 <= bin_R13 <= 1
 0 <= bin_E1_R13 <= 1
 0 <= bin_G03_or_G08 <= 1
 0 <= bin_G23 <= 1
 0 <= bin_G03 <= 1
 0 <= bin_G08 <= 1
 0 <= bin_G08 <= 1
 0 <= bin_G07 <= 1
Binaries
 bin_R01  bin_E1_R01  bin_G01  bin_G19  bin_R16  bin_E1_R16  bin_G16  bin_G25 
 bin_G37  bin_R03  bin_E1_R03  bin_G23  bin_E2_R03  bin_G03  bin_R02 
 bin_E1_R02  bin_G02  bin_G20  bin_G21  bin_R05  bin_E1_R05  bin_G05  bin_G25 
 bin_G26  bin_E2_R05  bin_G16  bin_G25  bin_G37  bin_R04  bin_E1_R04  bin_G04 
 bin_G24  bin_R10  bin_E1_R10  bin_G10  bin_G31  bin_G32  bin_R06  bin_E1_R06 
 bin_G06  bin_G27  bin_G28  bin_R09  bin_E1_R09  bin_G09  bin_G30  bin_R08 
 bin_E1_R08  bin_G08  bin_G25  bin_G26  bin_R11  bin_E1_R11  bin_G11  bin_G22 
 bin_G21  bin_E2_R11  bin_G02  bin_G20  bin_G21  bin_R14  bin_E1_R14  bin_G14 
 bin_G34  bin_G35  bin_R18  bin_E1_R18  bin_G18  bin_G08_or_G07  bin_R15 
 bin_E1_R15  bin_G15  bin_G36  bin_R12  bin_E1_R12  bin_G12  bin_G33  bin_R17 
 bin_E1_R17  bin_G17  bin_G38  bin_R13  bin_E1_R13  bin_G03_or_G08  bin_G23 
 bin_G03  bin_G08  bin_G08  bin_G07 
End
