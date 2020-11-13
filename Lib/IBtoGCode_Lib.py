import math
import mmap
import os
import numpy as np
import pandas as pd
from scipy.signal import medfilt, savgol_filter


def test_if_string_helper(path):
    with open(path, 'rb', 0) as file, \
            mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
        # hier muss lediglich der Rückgabewert der Funktion abgefragt werden.
        if s.find(b',') != -1:
            return True
        else:
            return False


def trim_csv_helper(path, versuch, directory):
    with open(path, 'r') as file:
        filedata = file.read().replace(',', '.')
    # CSV wird abgespeichert
    path = os.path.normpath(os.path.join(directory, versuch + '_ed.csv'))
    with open(path, 'w') as file:
        file.write(filedata)
    pass


def init_data_helper(path, ang):
    first_line = pd.read_csv(path, sep=";", nrows=1, encoding='latin-1')
    if 'Aufnahmezeit' in first_line.columns:
        spalten = ['Aufnahmezeit', 'Temperatur', 'Sollwert', 'P-Ausgabe']
        print('Arbeite mit Aufnahmezeit')
        data_ini = pd.read_csv(path, sep=";", usecols=spalten, encoding='latin-1')
        start = data_ini.loc[data_ini['Aufnahmezeit'] == '00:00:00.000'].index[0]
        data_ini = data_ini.loc[start:]
    else:
        spalten = ['Zeit', 'Temperatur', 'Sollwert', 'P-Ausgabe']
        data_ini = pd.read_csv(path, sep=";", usecols=spalten, encoding='latin-1')

    df = pd.DataFrame(data_ini)
    # Zeit-Spalte ins Datumsformat konvertieren
    if 'Aufnahmezeit' in first_line.columns:
        df['Zeit'] = pd.to_datetime(df['Aufnahmezeit'])
    else:
        df['Zeit'] = pd.to_datetime(df['Zeit'])

    # Löschen des Startbereichs, falls vorhanden
    if any(df['Sollwert'].isin(['-----'])):
        df = df.loc[df['Sollwert'] != '-----']
    # Sollwertspalte wird in numerische Werte zurückgewandelt
    df.iloc[:, 2] = pd.to_numeric(df['Sollwert'])
    # self.df['Sollwert'] = self.df['Sollwert'].convert_objects(convert_numeric=True)
    # Sollwertspalte wird in numerische Werte zurückgewandelt
    # self.df['dT']=self.df['Sollwert'].diff()

    # Alle Zeilen, vor und nach der Regelung
    # (wenn Solltemp. < Grenztemp.) werden abgeschnitten
    df = df[df['Sollwert'] >= 550]

    # Setzen des Nullpunktes der Zeitmessung
    t1 = df['Zeit'].iloc[0]
    df['Zeit2'] = df['Zeit'] - t1

    # Erstellung der Winkelmatrix
    # (Gesamtlänge der Winkel wird auf komplettes Prozessfenster verteilt)
    matrix = np.linspace(0, ang, num=len(df['Temperatur']))
    # Erzeugte Matrix wird zu Dataframe hinzugefügt
    # (für leichteren Zugriff auf Daten)
    df['Winkel'] = matrix

    return df


def out_to_np_helper(laengen, werte):
    slp_out = []
    for lngn, wrt in zip(laengen, werte):
        try:
            slp_out.append([int(lngn), int(wrt) * 0.01])
        except ValueError:
            continue
    if len(slp_out) == 0:
        print('Hier fehlen Werte')
    return slp_out


def proc_path_helper(path):
    file = os.path.basename(path)
    versuch = os.path.splitext(file)[0]
    return file, versuch


def create_cnc_helper(ang, d_ang, df):
    # Anzahl der Zwischenschritte, aus der Strahlstromkurve besteht
    j = math.floor(len(df['Temperatur']) / (ang / d_ang))
    # jeder j. Schritt wird aus der geglätteten Strahlstromkurve übernommen
    cnc = df.iloc[::j, 5:]
    cnc = cnc.append(df.iloc[-1, 5:])
    # Der berechnete Strahlstrom wird auf 1 Stelle hinter dem Komma gerundet
    cnc['IB-korr2'] = round(cnc['IB-korr'], 1)
    cnc['SL'] = cnc['SL'].astype(int)
    return cnc


def filter_data_helper(met_med_checked: float, n_med_text: str, met_sav_checked: float, n_sav_text: str, p_sav_text: str, df):
    if met_med_checked:
        # Medianfilter über 99 Werte
        n_med = int(n_med_text)
        df['IB-korr'] = medfilt(df['IB'], n_med)
    elif met_sav_checked:
        n_sav = int(n_sav_text)
        poly_sav = int(p_sav_text)
        # Savgol-Filter über 401 Werte mit Polynom 5. Grades
        df['IB-korr'] = savgol_filter(df['IB'], n_sav, poly_sav)
    else:
        print('Keine/falsche Filtermethode angegeben')
    return df


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


def calc_sq_helper(df, i_0, foc_of, foc_ruhe, vs, slp_in, slp_out, ang, v_s_slope, di):
    sq = pd.DataFrame()

    '''
    Winkelangaben werden aus df übernommen (vorsicht: auch Indexe)
    Ansprechen über .loc (spricht Index an) nicht sinnvoll
    Besser: indizierung via .iloc (kann jedoch nicht mit Spaltennamen arbeiten)
    '''
    sq['Winkel'] = df['Winkel']

    # Die Spalten werden mit den Grundwerten gefüllt
    sq['SQ'] = np.repeat(i_0, len(sq))
    sq['SL'] = np.repeat(foc_of, len(sq))
    # Hier wird die Geschwindigkeitsspalte gefüllt
    sq['vs'] = vs
    ''' --- Slope IN ---'''
    # Der Slope-In wird für Winkel bestimmt, die <= dem Grenzwinkel sind
    slope_in = np.array(np.where(sq['Winkel'] <= slp_in))
    zaehl1 = slope_in.shape[1]  # Index Letztes Element Slope-In
    # SQ
    sq.iloc[:zaehl1, 1] = np.linspace(0, i_0, num=zaehl1)  # Slope IN SQ
    # SL
    sq.iloc[:zaehl1, 2] = np.linspace(foc_ruhe, foc_of, num=zaehl1)

    ''' --- SLOPE OUT --- '''
    pos_slope = get_slope_pos(slp_out)
    # Berechnung der Elemente von slp_out
    # (falls Anpassung der Slopepositionen notwendig)
    n = len(slp_out)
    # Initialisierung von Hilfvariablen
    slope_out = np.zeros(n)  # Position des Slopebeginns
    zaehl2 = np.zeros(n)  # Indexvariable
    for i in range(n):
        thresh = ang - int(pos_slope[i])
        # da np.where Tupel zurückgibt, muss die Anzahl der
        # betroffenen Elemente separat bestimmt werden
        slope_out[i] = np.array(np.where(sq['Winkel'] >= thresh)).shape[1]
        zaehl2[i] = len(sq) - slope_out[i]  # Index 1. Element Slope Out

    # Berechnung des Strahlstroms, auf den gesloped werden soll
    if n > 1:
        sq_out = slp_out[:, 1] * i_0
    else:
        sq_out = [slp_out[0][1] * i_0]

    sq_out.insert(0, i_0)

    # Berechnung der Zwischenschritte zwischen
    anz_out = np.abs(np.diff(slope_out))
    anz_out = np.append(anz_out, slope_out[-1])
    anz_out = anz_out.astype(int)
    # Zähler in integer konvertieren, um Typenverträglichkeit zu gewährleisten
    zaehl2 = zaehl2.astype(int)
    # SQ

    for i in range(len(zaehl2)):
        if (len(zaehl2) == 1) or (i == len(zaehl2)):
            print('ich mache das')
            sq.iloc[zaehl2[i]:, 1] = np.linspace(sq_out[i], sq_out[i + 1], anz_out[i])
        else:
            print('nein, das mach ich')
            sq.iloc[zaehl2[i]:zaehl2[i + 1], 1] = np.linspace(sq_out[i], sq_out[i + 1], anz_out[i])

    # SL
    sq.iloc[zaehl2[-1]:, 2] = np.linspace(foc_of, foc_ruhe, anz_out[-1])
    # FS
    sq.iloc[zaehl2[-1]:, 3] = np.repeat(v_s_slope, anz_out[-1])

    df['SQ'], df['SL'], df['vs'] = sq['SQ'], sq['SL'], sq['vs']

    df['IB'] = df['SQ'] + (df['P-Ausgabe'] * 0.01 * di)
    return df


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
