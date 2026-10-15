"""Measured optical constants at 550 nm, read out of nk_JH_total.mat
(index 151 of the 400-800 nm grid)."""
AG = 0.04382 + 3.818978372838121j              # material.l_Ag_McPeak
AL = 0.958336812040867 + 6.68678039868499j     # material.l_Al_JO
ITO = 1.86362074031846 + 0.00322847509359885j  # material.l_ITO  (Koenig 2014)

# uniaxial ETL / transport materials: (n_o, n_e)
ETL = {
    'B3PyMPM': (1.820533808087541, 1.608639620557389),   # material.l_B3_o_JO / l_B3_e_JO
    'B4PyMPM': (1.82448446278145,  1.56003496787215),    # material.l_B4_o  / l_B4_e
    'TPBi':    (1.73909863830168,  1.71380169903695),    # material.TPBi_o  / TPBi_e
    'TCTA':    (1.803789230108173, 1.71407015059591),    # material.l_TCTA_o_JO / l_TCTA_e_JO
}
