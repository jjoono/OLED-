"""Pull the optical constants and emission spectra used in this study out of the group's
nk library (nk_JH_total.mat, the library used for the round-trip-loss paper) into a small
npz on a 400-700 nm, 1 nm grid.  The library itself is not redistributed; pass its path."""
import sys, numpy as np, scipy.io as sio

LAM = np.arange(400, 701)
MATS = {'Ag': 'l_Ag_McPeak', 'Al': 'l_Al_JO', 'ITO': 'l_ITO', 'IZO': 'l_IZO', 'ZnS': 'l_ZnS', 'LiF': 'l_LiF',
        'TiO2': 'l_TiO2_SJ_RTP', 'SiO2': 'l_SiO2_SJ_RTP', 'SiNx': 'Si3N4', 'Glass': 'l_Glass',
        'B3_o': 'l_B3_o_JO', 'B3_e': 'l_B3_e_JO', 'EML_o': 'l_TCTA_B3_o_JO', 'EML_e': 'l_TCTA_B3_e_JO',
        'TCTA_o': 'l_TCTA_o_JO', 'TCTA_e': 'l_TCTA_e_JO', 'TAPC_o': 'l_TAPC_o_JO', 'TAPC_e': 'l_TAPC_e_JO'}
SPECS = {'Ir(ppy)2acac': 'l_I_Irppy2acac', 'Ir(dmppy-ph)2tmd': 'l_I_Irdmppyph2tmd',
         'v-DABNA (HF)': 'l_mCPCN_SpiroAC_20wt_v_DABNA_1_0wt', 'blue QD': 'I_blue_QD'}

if __name__ == '__main__':
    m = sio.loadmat(sys.argv[1], squeeze_me=True, struct_as_record=False)
    mat, spec = m['material'], m['spectrum']
    out = {'lam': LAM}
    for k, v in MATS.items():
        out['n_' + k] = np.asarray(getattr(mat, v)).astype(complex)[:301]
    for k, v in SPECS.items():
        s = np.clip(np.asarray(getattr(spec, v), float)[:301], 0, None)
        out['s_' + k] = s / s.max()
    np.savez('data/inputs.npz', **out)
    for k in SPECS:
        s = out['s_' + k]; pk = LAM[s.argmax()]; hm = LAM[s >= 0.5]
        print(f'{k:18s} peak {pk} nm  FWHM ~{hm.max() - hm.min()} nm')
