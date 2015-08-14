@model:3.1.1=potassium

@units
 substance  = mole                             "substance"
 volume     = litre                            "volume"
 area       = metre : e=2                      "area"
 length     = metre                            "length"
 time       = second                           "time"

@compartments
 default    =  1.0
 c1<default = 10.2                             "intracell"
 c2<default = 22.8                             "extracell"

@species
      c2:[K_ECF]       =    0.17513682777      "K_ECF"
      c1:[K_ICF]       =    0.897325991903     "K_ICF"
 default:K_urin        =    1.0e-9             "K_urin"
 default:K_tiss        = 1509.6             s  "K_tiss"
 default:K_sal         =    4.53818901679   s  "K_sal"
 default:met_act       =    6.24885128313   s  "met_act"
# default:K_milk        =    0.0             s  "K_milk"
      c2:[K_ECF_mmol]  =    4.6                "K_ECF_mmol"
      c1:[K_ICF_mmol]  =   22.9                "K_ICF_mmol"
 default:K_git         =   43.1365559451    s  "K_git"
      c2:[Gluc_b]      =    0.546566462071     "Gluc_b"
      c2:[ins_b]       =   21.8532954203       "ins_b"
 default:Gluc_stor     = 3647.74674177s        "Gluc_stor"
 default:Gluc_prod     =   34.0337183015s      "Gluc_prod"
 default:[s24]         =    1.0e-9             "Gluc_use"
 default:sa4_degraded  =    1.0e-9          s  "sa4_degraded"
 default:src_metact    =    1.0e-9          s  "src_metact"
 default:snk_metact    =    1.0e-9          s  "snk_metact"
      c2:src_insnew    =    1.0e-9          s  "src_insnew"
      c2:snk_insnew    =    1.0e-9          s  "snk_insnew"
 default:src_Kgit      =    1.0e-9          s  "src_Kgit"
      c2:src_Glucblood =    1.0e-9          s  "src_Glucb"
 default:s23           =    1.0e-9          s  "src_Gep"

@parameters
                    p3 =   22.0                "p3"
                    p4 =    0.0001783          "p4"
                    p5 =    0.563              "p5"
                    p6 =    0.051              "p6"
                    p8 =    0.5994             "p8"
                    p9 =    0.4016             "p9"
                   p13 =    6.08               "p13"
                   p14 =    0.15639932         "p14"
                   p15 =   32.0                "p15"
                   p16 =   10.0                "p16"
                   p18 =    0.0223             "p18"
                   p19 =    0.1373             "p19"
                   p20 =    7.38               "p20"
                   p21 =    0.1085             "p21"
                   p22 =    0.25               "p22"
                   p23 =    0.58649745         "p23"
                   p24 =    0.1962811466       "p24"
                   p25 =    6.6085             "p25"
                   p26 =   49.0279             "p26"
                   p27 =    0.232634885        "p27"
                   p29 =    0.11729949         "p29"
                   p30 =    0.95               "p30"
                   p31 =    0.353              "p31"
                   p32 =    1.31               "p32"
                   p33 =    0.07819966         "p33"
##                    pH =    0.0             v  "pH"
                   p42 =    4.5                "src_metact"
                   p43 =   10.0                "T_src_metact"
                   p44 =    1.2                "snc_metact"
                   p45 =    3.0                "basic_metact"
                   p46 =    0.0854             "src_glucfeed"
                   p48 =    0.08               "scale_blood_gluc"
                   p49 =  800.0                "p49"
                   p41 =   20.0                "p41"
###                 pi =    3.14159265359
##                   dmi =    0.0             v  "dmi"
##                K_feed =    0.0             v  "K_feed"
##             Gluc_feed =    0.0             v  "Gluc_feed"
                   p52 =    1.4                "kmilk_content"
                   p53 =    0.01               "basal_KeKu"
                   p54 =    1.0                "dmi_percentage"
                   p55 =    1.0e-9             "Milk_in_l/h"
                   p56 =    0.01137            "perctK_in_dmi"
                   p34 =    0.0390983          "con"
                    p1 =   40.4                "p1"
                   p11 =  860.0                "p11"
                   p12 =    0.0105             "p12"
                   p17 =   16.0                "p17"
                   p28 =   25.0                "p28"
                   p37 =    0.0                "p37"
                   p39 =   30.0                "p39"
                   p50 =    0.4                "p50"
                    p7 =    0.925              "p7"
                   p51 =    0.1                "p51"
                   p57 = 3800.0                "p57"
                   p59 =    1.0                "p59"
                   p60 =    2.0                "p60"
                   p61 =    5.0                "p61"
                   p40 =    0.01               "p40"
                   p47 =   72.0                "p47"
                   p58 =    0.01               "p58"
                   p63 =  800.0                "p63"
                   p64 =    0.2                "p64"
                   p65 =    0.3                "p65"
                   p66 =    1.4                "p66"
                    p2 = 1600.0                "p2"
                   p10 = 1400.0                "p10"
                   p35 =  100.0                "p35"
                   p62 =    1.0                "p62"
                   p38 = 1509.6                "p38"
                   p36 =    1.0                "p36"

@rules
 pH = 7.5 - met_act / 40
 dmi = p54 * 487.5 * (1 - sin(pi * t / 12))
 K_feed = p56 * dmi
 Gluc_feed = p46 * dmi
 K_ECF_mmol = K_ECF / p34
 K_ICF_mmol = K_ICF / p34

@reactions
@r=re20 "DmMa"
 src_metact -> met_act
 p45 + p42 * pow(dmi, 2) / (pow(dmi, 2) + pow(p43, 2))
@r=re37 "GbGs"
 Gluc_b -> Gluc_stor : ins_b, Gluc_prod
 pow(Gluc_feed, 10) / (pow(Gluc_feed, 10) + pow(p1, 10)) * p12 * Gluc_prod * pow(p57, 10) / (pow(Gluc_stor, 10) + pow(p57, 10)) * ins_b * Gluc_b * pow(p36, 10) / (pow(p36, 10) + pow(p55, 10))
@r=re41 "GpGb"
 Gluc_prod -> Gluc_b
 (p37 + p39 * Gluc_prod) * (1 / (1 + pow(Gluc_b / p50, 10)))
@r=re39 "GpGs"
 Gluc_prod -> Gluc_stor
 p51 * Gluc_prod * (pow(p57, 10) / (pow(Gluc_stor, 10) + pow(p57, 10))) * pow(p62, 10) / (pow(p62, 10) + pow(p55, 10))
@r=re38 "GsGb"
 Gluc_stor -> Gluc_b
 1 / (1 + pow(Gluc_feed / p1, 10)) * p17 * (p7 - Gluc_b) * pow(Gluc_stor, 10) / (pow(Gluc_stor, 10) + pow(p35, 10)) + 1 / (1 + pow(p64 / p55, 10)) * pow(p65, 10) / (pow(p65, 10) + pow(Gluc_b, 10)) * p63 * pow(Gluc_stor, 10) / (pow(Gluc_stor, 10) + pow(p35, 10))
@r=re43 "GsGp"
 Gluc_stor -> Gluc_prod
 p60 * 1 / (1 + pow(Gluc_feed / p61, 5)) * pow(Gluc_stor, 10) / (pow(Gluc_stor, 10) + pow(p35, 10))
@r=re4 "KeKi"
 K_ECF -> K_ICF : ins_b
 (p8 + p9 * pow(ins_b, 8) / (pow(ins_b, 8) + pow(p3, 8))) * pow(K_ECF, 2) / (pow(K_ECF, 2) + pow(p33, 2)) * (1 + p21 * pow(pH, 10) / (pow(pH, 10) + pow(p20, 10)))
@r=re5 "KeKs"
 K_ECF -> K_sal
 p18 * K_ECF * dmi
@r=re7 "KeKt"
 K_ECF -> K_tiss : K_git
 pow(K_git, 10) / (pow(K_git, 10) + pow(p15, 10)) * p25 * K_ECF * pow(p2 - p10, 10) / (pow(K_tiss - p10, 10) + pow(p2 - p10, 10))
@r=re9 "KeKu"
 K_ECF -> K_urin : K_git
 (1 + p13 * pow(K_ECF, 5) / (pow(K_ECF, 5) + pow(p24, 5))) * p6 * K_git * (1 + p16 * pow(K_ECF, 10) / (pow(K_ECF, 10) + pow(p22, 10))) + p53 * K_ECF
#@r=re35 "KeMi"
# K_ECF -> K_milk
# p52 * p55 * pow(K_ECF, 10) / (pow(K_ECF, 10) + pow(p40, 10))
@r=re6 "KeR"
 K_ECF -> sa4_degraded
 p4 * K_ECF
@r=re34 "KfKg"
 src_Kgit -> K_git
 p30 * K_feed
@r=re10 "KgKe"
 K_git -> K_ECF
 p31 * K_git
@r=re3 "KiKe"
 K_ICF -> K_ECF
 (1 + 1 / (1 + pow(K_ECF / p29, 10))) * p5 * (1 + 1 / (1 + pow(pH / p20, 10))) * pow(K_ICF, 2) / (pow(K_ICF, 2) + pow(p23, 2)) * (1 + p19 * pow(K_ICF, 2) / (pow(K_ICF, 2) + pow(p14, 2)))
@r=re11 "KsKg"
 K_sal -> K_git
 p32 * K_sal
@r=re8 "KtKe"
 K_tiss -> K_ECF : K_git
 1 / (1 + pow(K_git / p15, 10)) * p26 * (p27 - K_ECF) * pow(K_tiss - p10, 10) / (pow(K_tiss - p10, 10) + pow(p38 * 0.99 - p10, 10))
@r=re21 "MaSnk"
 met_act -> snk_metact
 p44 * met_act
@r=re40 "SnkGb"
 Gluc_b -> s24 : Gluc_prod
 pow(Gluc_b, 10) / (pow(Gluc_b, 10) + pow(p58, 10)) * (p28 * Gluc_b + p59 * Gluc_prod * exp(-p66 * p55) + p55 * p47)
@r=re32 "SnkIn"
 ins_b -> snk_insnew
 p41 * ins_b
@r=re33 "SrcGb"
 src_Glucblood -> Gluc_b
 p48 * Gluc_feed
@r=re36 "SrcGp"
 s23 -> Gluc_prod
 (1 - p48) * Gluc_feed
@r=re31 "SrcIn"
 src_insnew -> ins_b : Gluc_b
 p49 * Gluc_b
