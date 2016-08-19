\ENCODING=ISO-8859-1
\Problem name: 1464012252.mip

Minimize
 rhs: - 10 l1 + 1000 l5 + 1000 l6 + 1000 l11 + 1000 l12 + 1000 l15 + 1000 u1
      + 100 u2 + 1000 u3 + 1000 u4 + 1000 u5 + 1000 u6 + 1000 u7 + 1000 u8
      + 1000 u9 + 1000 u10 + 1000 u11 + 1000 u12 + 1000 u13 + 1000 u14
      + 1000 u15 + 1000 u16 + 1000 u17 + 1000 u18
Subject To
 dR17: - pM5 - l1 + u1  = 1
 dR01: pADP - pATP + pM1 + u2 >= 0
 dR03: - pADP + pATP - pM2 + pM3 + u3 >= 0
 dR13: - pADP + pATP + pM11 - pM9 + u4 >= 0
 dR02: - pM1 + pM2 - l5 + u5  = 0
 dR10: - pM1 + pM9 - l6 + u6  = 0
 dR12: pM10 - pM9 + u7 >= 0
 dR14: - pM11 + 0.7 pM12 + 0.3 pM13 + u8 >= 0
 dR16: - pM12 + pM5 + u9 >= 0
 dR15: - pM13 + u10 >= 0
 dR11: pM2 - pM9 - l11 + u11  = 0
 dR04: - pM3 + pM4 - l12 + u12  = 0
 dR06: - pM3 + pM6 + u13 >= 0
 dR18: - pM3 + 2 pM8 + u14 >= 0
 dR05: - pM4 + pM5 - l15 + u15  = 0
 dR08: 0.5 pM5 - pM7 + 0.5 pM8 + u16 >= 0
 dR07: - pM6 + pM7 + u17 >= 0
 dR09: - pM8 + u18 >= 0
Bounds
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
End
