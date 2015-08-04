@model:3.1.1=PAEON_V2__OvF_ne_0


@compartments
  default = 1.0                           "body (volume?)"


@species
  default:[LH_pit]    = 259696.0          "LH in the pituitary"                                     "[IU]"
  default:[LH_blood]  =      2.71         "LH in the blood"                                         "[IU/l]"
  default:[R_LH]      =      8.398        "Receptors LH"                                            "[nmol/l]"
  default:[LH_R]      =      0.266        "LH-receptor-complex"                                     "[nmol/l]"
  default:[R_LH_des]  =      0.71         "Receptors LH internalised"                               "[nmol/l]"
  default:[FSH_pit]   =  59734.0          "FSH in the pituitary"                                    "[IU]"
  default:[FSH_blood] =      5.28         "FSH in the blood"                                        "[IU/l]"
  default:[R_FSH]     =      5.903        "Receptors FSH"                                           "[nmol/l]"
  default:[FSH_R]     =      0.796        "FSH-receptor-complex"                                    "[nmol/l]"
  default:[R_FSH_des] =      1.8          "Receptors FSH internalised"                              "[nmol/l]"
  default:[S_foll]    =      0.156        "Follicular sensitivity to LH"                            "[-]"
  default:[AF1]       =      2.593        "Antral follicle develop. stage 1"                        "[AF1]"
  default:[AF2]       =     23.14         "Antral follicle develop. stage 2"                        "[AF2]"
  default:[AF3]       =      0.492        "Antral follicle develop. stage 3"                        "[AF3]"
  default:[AF4]       =      0.0000161    "Antral follicle develop. stage 4"                        "[AF4]"
  default:[PrF]       =      0.2449       "Pre-ovulatory follicle stage"                            "[PrF]"
  default:[OvF]       =      1.0E-16      "Ovulatory follicle stage"                                "[OvF]"
  default:[Sc1]       =      1.0E-08      "Ovulatory scar 1"                                        "[Sc1]"
  default:[Sc2]       =      2.0E-06      "Ovulatory scar 2"                                        "[Sc2]"
  default:[Lut1]      =      2.0E-05      "Corpus luteum develop. stage 1"                          "[Lut1]"
  default:[Lut2]      =      0.0003       "Corpus luteum develop. stage 2"                          "[Lut2]"
  default:[Lut3]      =      0.003        "Corpus luteum develop. stage 3"                          "[Lut3]"
  default:[Lut4]      =      0.011        "Corpus luteum develop. stage 4"                          "[Lut4]"
  default:[E2]        =     41.36         "Estradiol"                                               "[pg/ml]"
  default:[P4]        =      1.978        "Progesterone"                                            "[ng/ml]"
  default:[InhA]      =      0.93         "Inhibin A"                                               "[IU/ml]"
  default:[InhB]      =     60.53         "Inhibin B"                                               "[pg/ml]"
  default:[InhA_tau]  =     71.14         "Inhibin A delay component"                               "[IU/ml]"
###  default:[freq]      =     11.15         "Frequency of GnRH pulses"                                "[1/d]"
###  default:[mass]      =      0.0012       "Mass of GnRH pulses"                                     "[nmol]"
  default:[GnRH_G]    =      0.0152       "GnRH"                                                    "[nmol/l]"
  default:[R_GnRH_a]  =      0.009        "active Receptors GnRH"                                   "[nmol/l]"
  default:[R_GnRH_i]  =      0.00096      "inactive Receptors GnRH"                                 "[nmol/l]"
  default:[GnRH_R_a]  =      0.000065     "active GnRH-receptor-complex"                            "[nmol/l]"
  default:[GnRH_R_i]  =      0.000059     "inactive GnRH-receptor-complex"                          "[nmol/l]"


@parameters
#
# 001 (LH_pit)
            p_001_001 =   7309.92         "basal LH synthesis rate constant"                        "[IU/d]"                   "    1000.0        15000.0       " 
            p_001_002 =   7309.92         "E2 promoted LH synthesis rate constant"                  "[IU/d]"                   "    1000.0        15000.0       "
            p_001_003 =    192.2          "threshold of E2 on LH synthesis"                         "[pg/ml]"                  "       1.0E-01      400.0       "
            p_001_004 =     10.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_001_005 =      2.371        "threshold of P4 on LH synthesis"                         "[ng/ml]"                  "       1.0E-01       16.0       "
            p_001_006 =      1.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_001_007 =      0.00476      "basal LH release rate constant"                          "[1/d]"                    "       1.0E-06        0.1       "
            p_001_008 =      0.1904       "GnRH_R_a promoted LH release rate constant"              "[1/d]"                    "       0.01           1.0       "
            p_001_009 =      0.0003       "threshold of GnRH_R_a on LH release"                     "[nmol/l]"                 "       1.0E-06        0.0004    " 
            p_001_010 =      5.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 002 (LH_blood)
            p_002_001 =      5.0          "blood volume"                                            "[l]"                      "       4.0           10.0       "
            p_002_002 =      2.143        "binding rate of LH to its receptor"                      "[l/(IU d)]"               "       0.1           20.0       "
            p_002_003 =     74.851        "clearance rate of LH from the blood"                     "[1/d]"                    "       5.0          300.0       "
#
# 003 (R_LH)
            p_003_001 =     68.949        "formation rate of free LH receptors"                     "[1/d]"                    "       5.0          500.0       "
#
# 004 (LH_R)
            p_004_001 =    183.36         "desensitisation rate of LH receptor complex"             "[1/d]"                    "      20.0          500.0       "
#
# 005 (R_LH_des)
#
# 006 (FSH_pit)
            p_006_001 =      2.213E+04    "basal FSH synthesis rate constant"                       "[IU/d]"                   "       1.0E+03        1.0E+05   "
            p_006_002 =     95.81         "threshold of InhA_tau in FSH synthesis"                  "[IU/ml]"                  "       1.0E-01      150.0       "
            p_006_003 =      5.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_006_004 =     70.0          "threshold of InhB in FSH synthesis"                      "[pg/ml]"                  "       1.0E-01      120.0       "
            p_006_005 =      3.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_006_006 =     10.0          "threshold of GnRH frequency in FSH synthesis"            "[1/d]"                    "       1.0E-01       20.0       "
            p_006_007 =      3.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_006_008 =      0.057        "basal FSH release rate constant"                         "[1/d]"                    "       0.01           1.0       "
            p_006_009 =      0.272        "stimulation of FSH release by GnRH rec."                 "[1/d]"                    "       0.05           1.0       "
            p_006_010 =      0.0003       "threshold of GnRH on FSH release"                        "[nmol/l]"                 "       1.0E-03        0.1       "
            p_006_011 =      3.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 007 (FSH_blood)
            p_007_001 =      3.529        "binding rate of FSH to its receptor"                     "[l/(IU d)]"               "       0.5           10.0       "
            p_007_002 =    114.25         "clearance rate of FSH from blood"                        "[1/d]"                    "      50.0          250.0       "  
#
# 008 (R_FSH)
            p_008_001 =     61.029        "formation rate of free FSH receptors"                    "[1/d]"                    "       5.0          500.0       "
#
# 009 (FSH_R)
            p_009_001 =    138.3          "desensitisation rate of FSH receptos complex"            "[1/d]"                    "      20.0         1000.0       "
#
# 010 (R_FSH_des)
#
# 011 (S_foll)
            p_011_001 =      0.219        "synth. rate const. of FSH_blood to S_foll"               "[[S_foll]/d]"             "       0.02           2.0       "
            p_011_002 =      3.0          "threshold of FSH_blood to stimulate S_foll"              "[nmol/l]"                 "       1.0E-01       10.0       "
            p_011_003 =      5.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_011_004 =      1.343        "S_foll clearance rate constant"                          "[1/d]"                    "       0.1           15.0       "
            p_011_005 =      1.235        "threshold of P4 on clearance of S_foll"                  "[ng/ml]"                  "       1.0E-01       16.0       "
            p_011_006 =      5.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_011_007 =      0.2          "additional clearance rate"                               "[]"                       "                                "
#
# 012 (AF1)
            p_012_001 =      3.663        "growth rate of AF1"                                      "[[AF1]/d]"                "       1.0           10.0       "
            p_012_002 =      0.608        "threshold of FSH_R for stimul. of AF1"                   "[nmol/l]"                 "       1.0E-01        1.2       "
            p_012_003 =      3.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_012_004 =      1.221        "transition rate constant from AF1 to AF2"                "[l/(d nmol)]"             "       0.1           10.0       "
#
# 013 (AF2)
            p_013_001 =      4.882        "transition rate constant from AF2 to AF3"                "[1/d]"                    "       0.5           20.0       "
            p_013_002 =      2.726        "scaling of LH receptor complex"                          "[nmol/l]"                 "       0.2           10.0       "
            p_013_003 =      3.689        "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 014 (AF3)
            p_014_001 =      0.122        "self growth rate of AF3"                                 "[l/(d nmol)]"             "       0.01           2.0       "
            p_014_002 =     10.0          "maximum size of AF3 and AF4"                             "[[AF3]]"                  "       1.0           50.0       "
            p_014_003 =    122.06         "transition rate constant from AF3 to AF4"                "[1/d]"                    "      10.0          500.0       "
            p_014_004 =      5.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 015 (AF4)
            p_015_001 =     12.206        "self growth rate of AF4"                                 "[1/d]"                    "       1.0           50.0       "
            p_015_002 =      2.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "    
            p_015_003 =    332.75         "transition rate constant from AF4 to PrF"                "[1/d]"                    "      30.0         1000.0       "
#
# 016 (PrF)
            p_016_001 =    122.06         "elimination rate of PrF"                                 "[1/d]"                    "      10.0          500.0       "
            p_016_002 =      6.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 017 (OvF)
            p_017_001 =      7.984        "growth rate of OvF"                                      "[1/d]"                    "       1.0           20.0       "
            p_017_002 =      3.0          "threshold of PrF for OvF formation"                      "[[PrF]]"                  "       1.0E-01        4.0       "
            p_017_003 =     10.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       " 
            p_017_004 =     12.206        "elimination rate constant of OvF"                        "[1/d]"                    "       1.0           50.0       "
#
# 018 (Sc1)
            p_018_001 =      1.208        "growth rate of Sc1 stimulated by OvF"                    "[1/d]"                    "       0.1           10.0       "
            p_018_002 =      0.02         "threshold of OvF to form Sc1"                            "[[OvF]]"                  "       1.0E-03        0.1       "
            p_018_003 =     10.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_018_004 =      1.221        "transition rate constant from Sc1 to Sc2"                "[1/d]"                    "       0.1           10.0       "
#
# 019 (Sc2)
            p_019_001 =      0.958        "transition rate from Sc2 to Lut1"                        "[1/d]"                    "       0.1           10.0       "
#
# 020 (Lut1)
            p_020_001 =      0.925        "transition rate constant from Lut1 to Lut2"              "[1/d]"                    "       0.1           10.0       "
            p_020_002 =     20.0          "rate of GnRH_R_a to the luteal development"              "[-]"                      "       2.0          100.0       "
            p_020_003 =      0.0008       "threshold of GnRH_R_a to the luteal development"         "[nmol/l]"                 "       0.0001         0.0004    "
            p_020_004 =      5.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 021 (Lut2)
            p_021_001 =      0.7567       "transition rate constant from Lut2 to Lut3"              "[1/d]"                    "       0.1           10.0       "
#
# 022 (Lut3)
            p_022_001 =      0.61         "transition rate constant from Lut3 to Lut4"              "[1/d]"                    "       0.1           10.0       "
#
# 023 (Lut4)
            p_023_001 =      0.543        "clearance rate constant of Lut4"                         "[1/d]"                    "       0.1           10.0       "
#
# 024 (E2)
            p_024_001 =     51.558        "basal E2 production"                                     "[pg/(ml d)]"              "       5.0          100.0       "
            p_024_002 =      2.0945       "production of E2 by AF2"                                 "[g/(mol [AF2])]"          "       0.2           10.0       "
            p_024_003 =      9.28         "production of E2 by LH and AF3"                          "[pg/(ml [AF3] [LH] d)]"   "       1.0          100.0       "
            p_024_004 =   3480.27         "production of E2 by AF4"                                 "[pg/(ml [AF4] d)]"        "     100.0        10000.0       "
            p_024_005 =      0.972        "production of E2 by LH and PrF"                          "[pg/(ml [AF3] [LH] d)]"   "       0.1           10.0       "
            p_024_006 =   1713.71         "production of E2 by Lut1"                                "[pg/(ml [Lut1] d)]"       "     100.0         5000.0       "
            p_024_007 =   8675.14         "production of E2 by Lut4"                                "[pg/(ml [Lut4] d)]"       "     100.0        10000.0       "
            p_024_008 =      5.235        "E2 clearance rate constant"                              "[1/d]"                    "       0.5           20.0       "
#
# 025 (P4)
            p_025_001 =      0.943        "basal P4 production"                                     "[ng/(ml d)]"              "       0.1           10.0       "
            p_025_002 =    761.64         "production of P4 by Lut4"                                "[ng/(ml [Lut4] d)]"       "      50.0         5000.0       "
            p_025_003 =      5.13         "P4 clearance rate constant"                              "[1/d]"                    "       1.0           20.0       "
#
# 026 (InhA)
            p_026_001 =      1.445        "basal InhA production"                                   "[IU/(ml d)]"              "       0.1           10.0       "
            p_026_002 =      2.285        "production of InhA by PrF"                               "[IU/(ml [PrF] d)]"        "       0.2           20.0       "
            p_026_003 =     60.0          "production of InhA by Sc1"                               "[pg/(ml [Sc1] d)]"        "       5.0          100.0       "
            p_026_004 =    180.0          "production of InhA by Lut1"                              "[pg/(ml [Lut1] d)]"       "      20.0          500.0       "
            p_026_005 =     28.211        "production of InhA by Lut2"                              "[pg/(ml [Lut2] d)]"       "      10.0          100.0       "
            p_026_006 =    194.07         "production of InhA by Lut3"                              "[pg/(ml [Lut3] d)]"       "      20.0          500.0       "
            p_026_007 =    114.25         "production of InhA by Lut4"                              "[pg/(ml [Lut4] d)]"       "      10.0          200.0       "
            p_026_008 =      4.287        "InhA clearance rate constant"                            "[1/d]"                    "       0.5           50.0       "
#
# 027 (InhB)
            p_027_001 =     89.943        "basal InhB production"                                   "[pg/(ml d)]"              "      10.0          500.0       "
            p_027_002 =    447.47         "production of InhB by AF2"                               "[pg/(ml [AF2] d)]"        "     100.0         1000.0       "
            p_027_003 = 132240.2          "production of InhB by AF3"                               "[pg/(ml [AF3] d)]"        "   10000.0       200000.0       "
            p_027_004 =    172.45         "InhB clearance rate constant"                            "[1/d]"                    "      20.0         1000.0       "
#
# 028 (InhA_tau)
            p_028_001 =      0.199        "clearance of InhA in delayed compartment"                "[1/d]"                    "       0.02          10.0       "
#
# 029 (freq)
            p_029_001 =     16.0          "mean GnRH_G frequency"                                   "[1/d]"                    "       0.1           10.0       "
            p_029_002 =      3.0          "threshold of P4 for inhibition of GnRH_G frequency"      "[ng/ml]"                  "       1.0E-01       16.0       "
            p_029_003 =      2.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_029_004 =      1.0          "factor for stimulation by E2"                            "[-]"                      "       0.1           10.0       "
            p_029_005 =    220.0          "threshold of E2 for stimulation of GnRH_G frequency"     "[pg/ml]"                  "       1.0E-01      400.0       "
            p_029_006 =     10.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 030 (mass)
            p_030_001 =      0.005593     "amount of GnRH released by one puls at E2 level height"  "[nmol]"                   "       0.0005         0.05      "
            p_030_002 =    220.0          "threshold of E2 for stimulation of GnRH mass"            "[pg/ml]"                  "       1.0E-01      400.0       "
            p_030_003 =      2.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
            p_030_004 =      9.6          "threshold of E2 for inhibition of GnRH mass"             "[pg/ml]"                  "       1.0E-01      400.0       "
            p_030_005 =      1.0          "Hill exponent"                                           "[-]"                      "       1.0           10.0       "
#
# 031 (GnRH_G)
            p_031_001 =    322.18         "binding rate of GnRH_G to its receptor"                  "[l/(d nmol)]"             "      20.0         1000.0       "
            p_031_002 =    644.35         "breakup rate of GnRH-receptor-complex"                   "[1/d]"                    "      20.0         1000.0       "
            p_031_003 =      0.447        "degradation rate of GnRH_G"                              "[1/d]"                    "       0.05           2.0       "
#
# 032 (R_GnRH_a)
            p_032_001 =      3.222        "rate of receptor inactivation"                           "[1/d]"                    "       0.2           20.0       "
            p_032_002 =     32.218        "rate of receptor activation"                             "[1/d]"                    "       2.0          200.0       "
#
# 033 (R_GnRH_i)
            p_033_001 =      0.00008949   "synthesis rate of inactive receptors"                    "[nmol/(l d)]"             "       0.000005       0.001     "
            p_033_002 =     32.218        "dissociation rate of inactive receptor-complex"          "[1/d]"                    "      20.0          500.0       "
            p_033_003 =      0.0895       "degradation rate of inactive receptors"                  "[1/d]"                    "       0.001          0.1       "
#
# 034 (GnRH_R_a)
            p_034_001 =     32.218        "rate of receptor-complex inactivation"                   "[1/d]"                    "       2.0          200.0       "
            p_034_002 =      3.222        "rate of receptor-complex activation"                     "[1/d]"                    "       0.2           20.0       "
#
# 035 (GnRH_R_i)
            p_035_001 =      0.00895      "degradation rate of inactive receptor-complex"           "[1/d]"                    "       0.001          0.1       "


@functions
#
  Hplus( y, thrs, n ) = pow( y/thrs, n ) / (1.0 + pow( y/thrs, n ))
#
  Hminus( y, thrs, n ) = 1.0 / (1.0 + pow( y/thrs, n ))


@rules
#
  @assign: freq = p_029_001 * Hminus( P4, p_029_002, p_029_003 ) * (1.0 + p_029_004 * Hplus( E2, p_029_005, p_029_006 ))
#
  @assign: mass = p_030_001 * (Hplus( E2, p_030_002, p_030_003 ) + Hminus( E2, p_030_004, p_030_005 ))


@reactions
#
#
#  Eq. 01 : LH in the pituitary (LH_pit)
#  ------ 
#
@r=re001
     () -> LH_pit
     (p_001_001 + p_001_002 * Hplus( E2, p_001_003, p_001_004 )) * Hminus( P4, p_001_005, p_001_006 )
#
@r=re002
     LH_pit -> ()
     (p_001_007 + p_001_008 * Hplus( GnRH_R_a, p_001_009, p_001_010 )) * LH_pit
#
#
#  Eq. 02 : LH in the blood (LH_blood)
#  ------
#
@r=re003
     () -> LH_blood
     (p_001_007 + p_001_008 * Hplus( GnRH_R_a, p_001_009, p_001_010 )) * LH_pit / p_002_001
#
@r=re004
     LH_blood + R_LH -> LH_R
     p_002_002 * LH_blood * R_LH
#
@r=re005
     LH_blood -> ()
     p_002_003 * LH_blood
#
#
#  Eq. 03 : Receptors LH (R_LH)
#  ------
#
@r=re006
     R_LH_des -> R_LH
     p_003_001 * R_LH_des
#
#
#  Eq. 04 : LH-receptor-complex (LH_R)
#  ------
#
@r=re007
     LH_R -> R_LH_des
     p_004_001 * LH_R
#
#
#  Eq. 05 : Receptors LH internalised (R_LH_des)
#  ------
#
#
#  Eq. 06 : FSH in the pituitary (FSH_pit)
#  ------
#
@r=re008
     () -> FSH_pit
     p_006_001 / (1.0 + pow( InhA_tau/p_006_002, p_006_003 ) + pow( InhB / p_006_004, p_006_005 )) * Hminus( freq, p_006_006, p_006_007 )
#
@r=re009
     FSH_pit -> ()
     (p_006_008 + p_006_009 * Hplus( GnRH_R_a, p_006_010, p_006_011 )) * FSH_pit
#
#
#  Eq. 07 : FSH in the blood (FSH_blood)
#  ------
#
@r=re010
     () -> FSH_blood
     (p_006_008 + p_006_009 * Hplus( GnRH_R_a, p_006_010, p_006_011 )) * FSH_pit / p_002_001
#
@r=re011
     FSH_blood + R_FSH -> FSH_R
     p_007_001 * R_FSH * FSH_blood
#
@r=re012
     FSH_blood -> ()
     p_007_002 * FSH_blood
#
#
#  Eq. 08 : Receptors FSH (R_FSH)
#  ------
#
@r=re013
     R_FSH_des -> R_FSH
     p_008_001 * R_FSH_des
#
#
#  Eq. 09 : FSH-receptor-complex (FSH_R)
#  ------
#
@r=re014
    FSH_R -> R_FSH_des
    p_009_001 * FSH_R
#
#
#  Eq. 10 : Receptors FSH internalised (R_FSH_des)
#  ------
#
#
#  Eq. 11 : Follicular sensitivity to LH (S_foll)
#  ------
#
@r=re015
     () -> S_foll
     p_011_001 * Hplus( FSH_blood, p_011_002, p_011_003 )
#
@r=re016
     S_foll -> ()
     p_011_004 * Hplus( P4, p_011_005, p_011_006 ) * S_foll
#
#
#  Eq. 12 : Antral follicle development stage 1 (AF1)
#  ------
#
@r=re017
     () -> AF1
     p_012_001 * Hplus( FSH_R, p_012_002, p_012_003 )
#
@r=re018
     AF1 -> AF2
     p_012_004 * FSH_R * AF1
#
#
#  Eq. 13 : Antral follicle development stage 2 (AF2)
#  ------
#
@r=re019
     AF2 -> AF3
     p_013_001 * pow( LH_R/p_013_002, p_013_003 ) * S_foll * AF2
#
#
#  Eq. 14 : Antral follicle development stage 3 (AF3)
#  ------
# 
@r=re020
     () -> AF3
     p_014_001 * FSH_R * AF3 * (1.0 - AF3/p_014_002)
@r=re021
     AF3 -> AF4
     p_014_003 * pow( LH_R/p_013_002, p_014_004 ) * S_foll * AF3
#
#
#  Eq. 15 : Antral follicle development stage 4 (AF4)
#  ------
#
@r=re022
     () -> AF4
     p_015_001 * pow( LH_R/p_013_002, p_015_002) * AF4 * (1.0 - AF4/p_014_002)
@r=re023
     AF4 -> PrF
     p_015_003 * (LH_R/p_013_002) * S_foll * AF4
#
#
#  Eq. 16 : Pre-ovulatory follicular stage (PrF)
#  ------
#
@r=re024
     PrF -> ()
     p_016_001 * pow( LH_R/p_013_002, p_016_002 ) * S_foll * PrF
#
#
#  Eq. 17 : Ovulatory folliclular stage (OvF)
#  ------
#
@r=re025
     () -> OvF
     p_017_001 * pow( LH_R/p_013_002, p_016_002 ) * S_foll * Hplus( PrF, p_017_002, p_017_003 )
#
@r=re026
     OvF -> ()
     p_017_004 * OvF
#
#
#  Eq. 18 : Ovulatory scar 1 (Sc1)
#  ------
#
@r=re027
     () -> Sc1
     p_018_001 * Hplus( OvF, p_018_002, p_018_003 )
#
@r=re028
     Sc1 -> Sc2
     p_018_004 * Sc1
#
#
#  Eq. 19 : Ovulatory scar 2 (Sc2)
#  ------
#
@r=re029
     Sc2 -> Lut1
     p_019_001 * Sc2
#
#
#  Eq. 20 : Corpus luteum development stage 1 (Lut1)
#  ------
#
@r=re030
     Lut1 -> Lut2
     p_020_001 * Lut1
@r=re031
     Lut1 -> ()
     p_020_001 * p_020_002 * Hplus( GnRH_R_a, p_020_003, p_020_004 ) * Lut1
#
#
#  Eq. 21 : Corpus luteum development stage 2 (Lut2)
#  ------
#
@r=re032
     Lut2 -> Lut3
     p_021_001 * Lut2
@r=re033
     Lut2 -> ()
     p_021_001 * p_020_002 * Hplus( GnRH_R_a, p_020_003, p_020_004 ) * Lut2
#
#
#  Eq. 22 : Corpus luteum development stage 3 (Lut3)
#  ------
#
@r=re034
     Lut3 -> Lut4
     p_022_001 * Lut3
@r=re035
     Lut3 -> ()
     p_022_001 * p_020_002 * Hplus( GnRH_R_a, p_020_003, p_020_004 ) * Lut3
#
#
#  Eq. 23 : Corpus luteum development stage 4 (Lut4)
#  ------
#
@r=re036
     Lut4 -> ()
     p_023_001 * (1.0 + p_020_002 * Hplus( GnRH_R_a, p_020_003, p_020_004 )) * Lut4
#
#
#  Eq. 24 : Estradiol blood level (E2)
#  ------
#
@r=re037
     () -> E2
     p_024_001 + p_024_002 * AF2 + p_024_003 * LH_blood * AF3 + p_024_004 * AF4 + p_024_005 * LH_blood * PrF + p_024_006 * Lut1 + p_024_007 * Lut4
#
@r=re038
     E2 -> ()
     p_024_008 * E2
#
#
#  Eq. 25 : Progesterone blood level (P4)
#  ------
#
@r=re039
     () -> P4
     p_025_001 + p_025_002 * Lut4
#
@r=re040
     P4 -> ()
     p_025_003 * P4
#
#
#  Eq. 26 : Inhibin A blood level (InhA)
#  ------
#
@r=re041
     () -> InhA
     p_026_001 + p_026_002 * PrF + p_026_003 * Sc1 + p_026_004 * Lut1 + p_026_005 * Lut2 + p_026_006 * Lut3 + p_026_007 * Lut4
#
@r=re042
     InhA -> InhA_tau
     p_026_008 * InhA
#
#
#  Eq. 27 : Inhibin B blood level (InhB)
#  ------
#
@r=re043
     () -> InhB
     p_027_001 + p_027_002 * AF2 + p_027_003 * Sc2
#
@r=re044
     InhB -> ()
     p_027_004 * InhB
#
#
#  Eq. 28 : Effective (i.e. delayed) inhibin A (InhA_tau)
#  ------
#
@r=re045
     InhA_tau -> ()
     p_028_001 * InhA_tau
#
#
#  Eq. 29 : GnRH frequency (freq)
#  ------
#  (cf. assignment rule)
#
#
#  Eq. 30 : GnRH mass (mass)
#  ------
#  (cf. assignment rule)
#
#
#  Eq. 31 : GnRH (GnRH_G)
#  ------
#
@r=re046
     () -> GnRH_G
     mass * freq
#
@r=re047
     GnRH_G + R_GnRH_a -> GnRH_R_a
     p_031_001 * R_GnRH_a * GnRH_G
#
@r=re048
     GnRH_R_a -> GnRH_G + R_GnRH_a
     p_031_002 * GnRH_R_a
#
@r=re049
      GnRH_G -> ()
      p_031_003 * GnRH_G
#
#
#  Eq. 32 : active Receptors GnRH (R_GnRH_a)
#  ------
#
@r=re050
     R_GnRH_a -> R_GnRH_i
     p_032_001 * R_GnRH_a
#
@r=re051
     R_GnRH_i -> R_GnRH_a
     p_032_002 * R_GnRH_i
#
#
#  Eq. 33 : inactive Receptors GnRH (R_GnRH_i)
#  ------
#
@r=re052
     () -> R_GnRH_i
     p_033_001
#
@r=re053
     GnRH_R_i -> R_GnRH_i
     p_033_002 * GnRH_R_i
#
@r=re054
     R_GnRH_i -> ()
     p_033_003 * R_GnRH_i
#
#
#  Eq. 34 : active GnRH-receptor-complex (GnRH_R_a)
#  ------
#
@r=re055
     GnRH_R_a -> GnRH_R_i
     p_034_001 * GnRH_R_a
#
@r=re056
     GnRH_R_i -> GnRH_R_a
     p_034_002 * GnRH_R_i
#
#
#  Eq. 35 : inactive GnRG-receptor-complex (GnRH_R_i)
#  ------
#
@r=re057
     GnRH_R_i -> ()
     p_035_001 * GnRH_R_i
#

