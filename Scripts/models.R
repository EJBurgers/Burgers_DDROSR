baseDDR <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * (1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2))
  })
}

baseOSR <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2 - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_E0 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_E1 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg+ct_M15*p21) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_E2 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/((Km_K1+ct_M23*P53P)**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_E3 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/((Km_K1+ct_M322*P53P)**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg+ct_M321*p21) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D0 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D1 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + sc_S*S/(Km_S + S) - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * (1 + (degradP5_byM2 * MDM2)) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * MDM2)) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D2 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * ((1/(1 + KiP5_sc*SRXN1)) + (degradP5_byM2 * MDM2)) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * MDM2)) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D3 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + sc_S*S/(Km_S + S) - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * ((1/(1 + KiP5_sc*SRXN1)) + (degradP5_byM2 * MDM2)) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * MDM2)) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D4 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + sc_S*S/(Km_S + S) - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * (1 + (degradP5_byM2 * ((MDM2**hillP5_byM2)/((KmP5_byM2**hillP5_byM2) + (MDM2**hillP5_byM2))))) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * ((MDM2**hillP5p_byM2)/((KmP5p_byM2**hillP5p_byM2) + (MDM2**hillP5p_byM2))))) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D5 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * ((1/(1 + KiP5_sc*SRXN1)) + (degradP5_byM2 * ((MDM2**hillP5_byM2)/((KmP5_byM2**hillP5_byM2) + (MDM2**hillP5_byM2))))) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * ((MDM2**hillP5p_byM2)/((KmP5p_byM2**hillP5p_byM2) + (MDM2**hillP5p_byM2))))) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_D6 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + sc_S*S/(Km_S + S) - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * ((1/(1 + KiP5_sc*SRXN1)) + (degradP5_byM2 * ((MDM2**hillP5_byM2)/((KmP5_byM2**hillP5_byM2) + (MDM2**hillP5_byM2))))) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * ((MDM2**hillP5p_byM2)/((KmP5p_byM2**hillP5p_byM2) + (MDM2**hillP5p_byM2))))) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_CE <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg+ct_M15*p21) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_CD <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + sc_S*S/(Km_S + S) - degradD_byP5 * DD * P53P
    dP53 = buildP5 + dephos * P53P - phos * P53 * DD - degradP5 * (1 + (degradP5_byM2 * ((MDM2**hillP5_byM2)/((KmP5_byM2**hillP5_byM2) + (MDM2**hillP5_byM2))))) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * ((MDM2**hillP5p_byM2)/((KmP5p_byM2**hillP5p_byM2) + (MDM2**hillP5p_byM2))))) * P53P
    dMDM2 = buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2)) - degradM2 * MDM2
    dp21 = buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*(1)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*(1)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = buildN2  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg+ct_M15*p21) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    #
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_CEKD <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + S - degradD_byP5 * DD * P53P
    dP53 = KD_P53*(buildP5) + dephos * P53P - phos * P53 * DD - degradP5*(1 + degradP5_byM2 * MDM2) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p*(1 + degradP5p_byM2 * MDM2) * P53P
    dMDM2 = KD_M2*(buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2))) - degradM2 * MDM2
    dp21 = KD_P21*(buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n)) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*S*(0)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = KD_N2*(buildN2)  - (VmaxN2_deg*KEAP1*(NRF2))/(((KmN2_deg+ct_M15*p21) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

Model_CDKD <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1 * S
    dDD = buildD + sc_S*S/(Km_S + S) - degradD_byP5 * DD * P53P
    dP53 = KD_P53*(buildP5) + dephos * P53P - phos * P53 * DD - degradP5 * (1 + (degradP5_byM2 * ((MDM2**hillP5_byM2)/((KmP5_byM2**hillP5_byM2) + (MDM2**hillP5_byM2))))) * P53
    dP53P = phos * P53 * DD - dephos * P53P - degradP5p * (1 + (degradP5p_byM2 * ((MDM2**hillP5p_byM2)/((KmP5p_byM2**hillP5p_byM2) + (MDM2**hillP5p_byM2))))) * P53P
    dMDM2 = KD_M2*(buildM2 + (VmaxM2_byP5p * P53P**n) / (KmM2_byP5p**n + P53P**n) + ((VmaxM2_byN2 * (NRF2**hillM2_byN2))/(KmM2_byN2**hillM2_byN2 + NRF2**hillM2_byN2))) - degradM2 * MDM2
    dp21 = KD_P21*(buildP2 + (VmaxP2 * P53P**n) / (KmP2**n + P53P**n)) - degradP2 * p21
    dBTG2 = buildB2 + (VmaxB2 * P53P**n) / (KmB2**n + P53P**n) - degradB2 * BTG2
    dKEAP1 = buildK1 + ((Vmax_K1*(NRF2**hillK1))/(Km_K1**hillK1 + NRF2**hillK1)) - (k_K1_modification*(S)*(1)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* KEAP1
    dmodifiedK1 = (k_K1_modification*(S)*(1)*KEAP1 - k_K1_unmodification*modifiedK1) - degradK1* modifiedK1
    dNRF2 = KD_N2*(buildN2)  - ((VmaxN2_deg*KEAP1*(NRF2))/((KmN2_deg+ct_M15*p21) + (NRF2))) - degradN2 * NRF2
    dSRXN1 = buildS1 + ((Vmax_S1 * (NRF2**hillS1))/(Km_S1**hillS1 + NRF2**hillS1)) - degradS1 * SRXN1
    #
    list(c(dS, dDD, dP53, dP53P, dMDM2, dp21, dBTG2, dKEAP1, dmodifiedK1, dNRF2, dSRXN1))
  })
}

# Model_KD
