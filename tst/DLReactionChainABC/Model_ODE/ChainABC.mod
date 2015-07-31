@model:3.1.1=ChainABC


@compartments
  default = 1.0                   "standard"


@species
  default:[A]    =    1.0         "Initial Concentration A" 
  default:[B]    =    0.0         "Initial Concentration B" 
  default:[C]    =    0.0         "Initial Concentration C"


@parameters
#
# Reaction Rates
            k1   =    2.0         "Reaction Rate A -> B"
            k_1  =    0.003       "Reaction Rate B -> A"
            k2   =    1.0         "Reaction Rate B -> C"
            k_2   =   0.002       "Reaction Rate C -> B"


@reactions
#
#
#  Eq. 01 :  A <--> B
#  ------ 
#
@r=re001
     A -> B
     k1 * A
#
@r=r002
     B -> A
     k_1 * B
#
#
#  Eq. 02 : B <--> C
#  ------
#
@r=re003
     B -> C
     k2 * B
#
@r=re004
     C -> B
     k_2 * C
#
