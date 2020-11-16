import os
import numpy as np
import pandas as pd


def get_vs(df, rb_vs_const: bool, rb_vs_step: bool, par_vs: str, par_vs_gk: str, par_vs_ns: str,
           par_slp_in: str, par_l_gk: str, par_l_ns: str, par_vs_trans: str):
    vf = pd.DataFrame()
    vf['Winkel'] = df['Winkel']
    vf = vf.reset_index(drop=True)

    vf['vs'] = int(par_vs)

    if rb_vs_const:
        vf['vs'] = int(par_vs)
    elif rb_vs_step:
        v_s_gk = int(par_vs_gk)
        v_s_nh = int(par_vs_ns)
        slp_in = int(par_slp_in)
        l_gk = slp_in + int(par_l_gk)
        l_nh = int(par_l_ns)
        l_trans = int(par_vs_trans)

        # bestimme die Positionen, an denen nur die Kontur aufhört
        pos_nh_a = l_gk + l_trans
        pos_nh_e = pos_nh_a + l_nh

        # Bestimme die Positionen, an denen der Übergang stattfindet
        pos_trans_1 = l_gk
        pos_trans_2 = pos_nh_e + l_trans

        # initialisiere vs-Spalte
        vf['vs'] = v_s_gk

        # ändere Werte, die im Bereich der Nockenkontur liegen
        idx = vf[vf['Winkel'] <= pos_nh_a].idxmax().values[0]
        idx_2 = vf[vf['Winkel'] > pos_nh_e].idxmin().values[0]

        vf['vs'].loc[idx:idx_2] = v_s_nh

        # Wie viele Punkte liegen im Übergang?
        idx_3 = vf[vf['Winkel'] <= pos_trans_1].idxmax().values[0]
        idx_4 = vf[vf['Winkel'] > pos_trans_2].idxmin().values[0]

        anz_trans_1 = len(vf['vs'].loc[idx_3:idx])
        anz_trans_2 = len(vf['vs'].loc[idx_2:idx_4])

        # Ersetze den Übergang zwischen Kontur und GK mit einer linearen Interpolation
        vf['vs'].loc[idx_3:idx] = np.linspace(v_s_gk, v_s_nh, anz_trans_1)
        vf['vs'].loc[idx_2:idx_4] = np.linspace(v_s_nh, v_s_gk, anz_trans_2)
    return vf['vs'].values


def get_slope_pos(slp_out):
    """
    :param slp_out:
    :return: pos
    """
    # Aufaddieren der Slopepositionen, um Startpositionen zu berechnen
    # zuerst umkehren der Matrix slp_out, da cumsum positiv aufaddiert
    if len(slp_out) > 1:
        pos = np.flipud(slp_out[:, 0])
        # Bildung von cumsum, danach umkehren der Matrix, um ursprüngliche
        # Reihenfolge nicht durcheinander zu bringen
        pos = np.flipud(np.cumsum(pos))
    else:
        pos = [slp_out[0][0]]
    return pos


def save_csv(df, versuch, directory):
    file_out = versuch + '_out.csv'
    fname = os.path.join(directory + os.sep, file_out)
    df_out = df[['Temperatur', 'Winkel', 'IB', 'IB-korr']].copy()
    df_out.to_csv(path_or_buf=fname, sep=';', decimal=',',
                  encoding='utf-8')


