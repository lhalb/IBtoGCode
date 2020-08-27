# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 14:13:17 2018

@author: halbauer
"""

from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow, QFileDialog
import sys   # We need sys so that we can pass argv to QApplication
import numpy as np
import pandas as pd
import os.path
import mmap
import math
import matplotlib.pyplot as plt
from scipy.signal import medfilt, savgol_filter

import axesMPL

class MyApp(QMainWindow, axesMPL.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)   # Laden der UI-Datei
        self.setup_triggers()   # lade Events für Buttons
        self.statusBar().showMessage('Programm erfolgreich geladen.')

    def setup_triggers(self):
        self.but_trim_csv.clicked.connect(self.trim_csv)
        self.but_init_data.clicked.connect(self.init_data)
        self.but_calc_ib.clicked.connect(self.calc_sq)
        self.but_filter_data.clicked.connect(self.filter_data)
        self.but_create_cnc.clicked.connect(self.create_cnc)
        self.but_QuickPlot.clicked.connect(self.quick_plot)
        self.but_save_data.clicked.connect(self.save_data)
        self.but_replot_data.clicked.connect(self.replot)
        self.actionOpen_File.triggered.connect(self.file_open)
        self.but_open_file.clicked.connect(self.file_open)
        self.filt_met_sav.toggled.connect(self.enable_butt)
        self.filt_met_med.toggled.connect(self.disable_butt)
        self.but_out.clicked.connect(self.file_save)


    # Funktionen, die in ein extra Modul gehören
    
    def test_if_string(self):
        with open(self.path, 'rb', 0) as file, \
        mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
            # hier muss lediglich der Rückgabewert der Funktion abgefragt werden. 
            if s.find(b',') != -1:
                self.statusBar().showMessage('Datei enthält ungültige Zeichen!')
                # self.status_line.setStyleSheet("color: rgb(181, 18, 62);")
                self.but_trim_csv.setEnabled(True)
            else:
                self.statusBar().showMessage('Datei erfolgreich geladen')
                # self.status_line.setStyleSheet("color: rgb(0, 0, 0);")
                self.but_init_data.setEnabled(True)
                self.but_QuickPlot.setEnabled(True)


    def trim_csv(self):
        with open(self.path, 'r') as file:
            filedata = file.read()

        filedata = filedata.replace(',', '.')
        # CSV wird abgespeichert
        self.file = self.versuch + '_ed.csv'
        self.path = os.path.join(self.dir, self.file)
        self.path = os.path.normpath(self.path)
        with open(self.path, 'w') as file:
            file.write(filedata)
            self.proc_path()
        self.test_if_string()

    def init_data(self):
        first_line = pd.read_csv(self.path, sep=";", nrows=1, encoding='latin-1')
        if 'Aufnahmezeit' in first_line.columns:
            spalten = ['Aufnahmezeit', 'Temperatur', 'Sollwert', 'P-Ausgabe']
            print('Arbeite mit Aufnahmezeit')
            
            self.data_ini = pd.read_csv(self.path, sep=";", usecols=spalten, encoding='latin-1')
            start = self.data_ini.loc[self.data_ini['Aufnahmezeit'] == '00:00:00.000'].index[0]
            self.data_ini = self.data_ini.loc[start:]

        else:
            spalten = ['Zeit', 'Temperatur', 'Sollwert', 'P-Ausgabe']
            self.data_ini = pd.read_csv(self.path, sep=";", usecols=spalten, encoding='latin-1')
            
        self.df = pd.DataFrame(self.data_ini)
        # Zeit-Spalte ins Datumsformat konvertieren
        if 'Aufnahmezeit' in first_line.columns:
            self.df['Zeit'] = pd.to_datetime(self.df['Aufnahmezeit'])
        else:
            self.df['Zeit'] = pd.to_datetime(self.df['Zeit'])


        # Löschen des Startbereichs, falls vorhanden
        if any(self.df['Sollwert'].isin(['-----'])):
            self.df = self.df.loc[self.df['Sollwert'] != '-----']
        # Sollwertspalte wird in numerische Werte zurückgewandelt
        self.df.iloc[:, 2] = pd.to_numeric(self.df['Sollwert'])
        # self.df['Sollwert'] = self.df['Sollwert'].convert_objects(convert_numeric=True)
        # Sollwertspalte wird in numerische Werte zurückgewandelt
        # self.df['dT']=self.df['Sollwert'].diff()

        # Alle Zeilen, vor und nach der Regelung
        # (wenn Solltemp. < Grenztemp.) werden abgeschnitten
        self.df = self.df[self.df['Sollwert'] >= 550]

        # Setzen des Nullpunktes der Zeitmessung
        t1 = self.df['Zeit'].iloc[0]

        self.df['Zeit2'] = self.df['Zeit'] - t1

        # Erstellung der Winkelmatrix
        # (Gesamtlänge der Winkel wird auf komplettes Prozessfenster verteilt)
        ang = int(self.par_ang.text())
        matrix = np.linspace(0, ang, num=len(self.df['Temperatur']))
        # Erzeugte Matrix wird zu Dataframe hinzugefügt
        # (für leichteren Zugriff auf Daten)
        self.df['Winkel'] = matrix

        self.but_calc_ib.setEnabled(True)
        self.statusBar().showMessage('Datenmatrix erfolgreich erzeugt. Weiter mit CALC IB!', 2000)

    def out_to_np(self):
        laengen = [self.par_slp_out_0_1.text(), self.par_slp_out_0_2.text(), self.par_slp_out_0_3.text()]
        werte = [self.par_slp_out_1_1.text(), self.par_slp_out_1_2.text(), self.par_slp_out_1_3.text()]

        self.slp_out = []

        for l, w in zip(laengen, werte):
            try:
                self.slp_out.append([int(l), int(w)*0.01])
            except ValueError:
                continue

        if len(self.slp_out) == 0:
            print('Hier fehlen Werte')

        return self.slp_out

    def get_vs(self):
        vf = pd.DataFrame()
        vf['Winkel'] = self.df['Winkel']
        vf = vf.reset_index(drop=True)

        vf['vs'] = int(self.par_vs.text())

        if self.rb_vs_const.isChecked():
            vf['vs'] = int(self.par_vs.text())
        elif self.rb_vs_step.isChecked():
            v_s_gk = int(self.par_vs_gk.text())
            v_s_nh = int(self.par_vs_ns.text())
            slp_in = int(self.par_slp_in.text())
            l_gk = slp_in + int(self.par_l_gk.text())
            l_nh = int(self.par_l_ns.text())
            l_trans = int(self.par_vs_trans.text())

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
        
    def calc_sq(self):
        sq = pd.DataFrame()

        '''
        Winkelangaben werden aus df übernommen (vorsicht: auch Indexe)
        Ansprechen über .loc (spricht Index an) nicht sinnvoll
        Besser: indizierung via .iloc (kann jedoch nicht mit Spaltennamen arbeiten)
        '''
        sq['Winkel'] = self.df['Winkel']

        # Übernehmen der Einträge aus der GUI
        i_0 = int(self.par_i0.text())
        foc_of = int(self.par_il.text())
        

        slp_in = int(self.par_slp_in.text())
        foc_ruhe = foc_of + int(self.par_slp_out_foc.text())
        v_s_slope = int(self.par_slp_out_vs.text())

        ang = int(self.par_ang.text())

        slp_out = self.out_to_np()

        # Die Spalten werden mit den Grundwerten gefüllt
        sq['SQ'] = np.repeat(i_0, len(sq))
        sq['SL'] = np.repeat(foc_of, len(sq))
        # Hier wird die Geschwindigkeitsspalte gefüllt
        sq['vs'] = self.get_vs()
        ''' --- Slope IN ---'''
        # Der Slope-In wird für Winkel bestimmt, die <= dem Grenzwinkel sind
        slope_in = np.array(np.where(sq['Winkel'] <= slp_in))
        zaehl1 = slope_in.shape[1]    # Index Letztes Element Slope-In
        # SQ
        sq.iloc[:zaehl1, 1] = np.linspace(0, i_0, num=zaehl1)  # Slope IN SQ
        # SL
        sq.iloc[:zaehl1, 2] = np.linspace(foc_ruhe, foc_of, num=zaehl1)

        ''' --- SLOPE OUT --- '''
        pos_slope = self.get_slope_pos()
        # Berechnung der Elemente von slp_out
        # (falls Anpassung der Slopepositionen notwendig)
        n = len(slp_out)
        # Initialisierung von Hilfvariablen
        slope_out = np.zeros(n)   # Position des Slopebeginns
        zaehl2 = np.zeros(n)    # Indexvariable
        for i in range(n):
            thresh = ang-int(pos_slope[i])
            # da np.where Tupel zurückgibt, muss die Anzahl der
            # betroffenen Elemente separat bestimmt werden
            slope_out[i] = np.array(np.where(sq['Winkel'] >= thresh)).shape[1]
            zaehl2[i] = len(sq) - slope_out[i]      # Index 1. Element Slope Out

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
                sq.iloc[zaehl2[i]:, 1] = np.linspace(sq_out[i], sq_out[i+1], anz_out[i])
            else:
                print('nein, das mach ich')
                sq.iloc[zaehl2[i]:zaehl2[i+1], 1] = np.linspace(sq_out[i], sq_out[i+1], anz_out[i])
          
        # SL
        sq.iloc[zaehl2[-1]:, 2] = np.linspace(foc_of, foc_ruhe, anz_out[-1])
        # FS
        sq.iloc[zaehl2[-1]:, 3] = np.repeat(v_s_slope, anz_out[-1])

        self.df['SQ'], self.df['SL'], self.df['vs'] = \
        sq['SQ'], sq['SL'], sq['vs']

        di = int(self.par_di.text())
        self.df['IB'] = self.df['SQ'] + (self.df['P-Ausgabe']*0.01*di)

        self.but_filter_data.setEnabled(True)
        self.statusBar().showMessage('Strahlstromkurven berechnet. Weiter mit FILTER DATA!')

    def get_slope_pos(self):
        
        # Aufaddieren der Slopepositionen, um Startpositionen zu berechnen
        # zuerst umkehren der Matrix slp_out, da cumsum positiv aufaddiert
        if len(self.slp_out) > 1:
            pos = np.flipud(self.slp_out[:, 0])
            # Bildung von cumsum, danach umkehren der Matrix, um ursprüngliche
            # Reihenfolge nicht durcheinander zu bringen    
            pos = np.flipud(np.cumsum(pos))
        
        else:
            pos = [self.slp_out[0][0]]

        return pos

    def filter_data(self):
        # --- Filterung der Daten ---
        if self.filt_met_med.isChecked():
            # Medianfilter über 99 Werte
            n_med = int(self.par_n_med.text())
            self.df['IB-korr'] = medfilt(self.df['IB'], n_med)
        elif self.filt_met_sav.isChecked():
            n_sav = int(self.par_n_sav.text())
            poly_sav = int(self.par_p_sav.text())
            # Savgol-Filter über 401 Werte mit Polynom 5. Grades
            self.df['IB-korr'] = savgol_filter(self.df['IB'], n_sav, poly_sav)
        else:
            print('Keine/falsche Filtermethode angegeben')

        self.but_create_cnc.setEnabled(True)
        self.statusBar().showMessage('Daten erfolgreich gefiltert. Weiter mit CALC CNC oder PLOT DATA!')

    def create_cnc(self):
        ang = float(self.par_ang.text())
        d_ang = float(self.par_d_ang.text())
        # Anzahl der Zwischenschritte, aus der Strahlstromkurve besteht
        j = math.floor(len(self.df['Temperatur'])/(ang / d_ang))
        # jeder j. Schritt wird aus der geglätteten Strahlstromkurve übernommen
        self.cnc = self.df.iloc[::j, 5:]
        self.cnc = self.cnc.append(self.df.iloc[-1, 5:])
        # Der berechnete Strahlstrom wird auf 1 Stelle hinter dem Komma gerundet
        self.cnc['IB-korr2'] = round(self.cnc['IB-korr'], 1)
        self.cnc['SL'] = self.cnc['SL'].astype(int)
        self.but_save_data.setEnabled(True)


    # Hilfsfunktionen, die nur in der GUI benötigt werden

    def file_open(self):
        self.path = QFileDialog.getOpenFileName(self, 'Open File')[0]
        if not self.path:
            return
        
        self.dir = os.path.dirname(self.path)
        self.proc_path()
        self.test_if_string()

    def file_save(self):
        fname = f'{self.txt_versuch.text()}.MPF'
        out_file = QFileDialog.getSaveFileName(self, 'Speicherort wählen', fname, 'CNC-Files (*.MPF);;All Files (*)')[0]
        if not out_file:
            return
        self.txt_outfolder.setText(out_file)

    def proc_path(self):
        self.txt_path.setText(self.path)
        self.file = os.path.basename(self.path)
        self.txt_file.setText(self.file)
        self.versuch = os.path.splitext(self.file)[0]
        self.txt_versuch.setText(self.versuch)

    def quick_plot(self):
        self.init_data()
        self.calc_sq()
        self.filter_data()
        self.create_cnc()
        self.plot_data()
        self.but_replot_data.setEnabled(True)

    def save_data(self):
        cnc = self.cb_save_cnc.isChecked()
        csv = self.cb_save_csv.isChecked()
        if csv and cnc:
            self.save_cnc()
            self.save_csv()
            self.statusBar().showMessage('CSV und CNC erfolgreich gesichert.', 2000)
        elif csv:
            self.save_csv()
            self.statusBar().showMessage('CSV erfolgreich gesichert.', 2000)
        elif cnc:
            self.save_cnc()
            self.statusBar().showMessage('CNC-Code erfolgreich exportiert.', 2000)
        else:
            self.statusBar().showMessage('Keine Daten gesichert.', 2000)

    def save_csv(self):
        file_out = self.versuch + '_out.csv'
        fname = os.path.join(self.dir + os.sep, file_out)
        df_out = self.df[['Temperatur', 'Winkel', 'IB', 'IB-korr']].copy()
        df_out.to_csv(path_or_buf=(fname), sep=';', decimal=',',
                  encoding='utf-8')

    def save_cnc(self):
        pos_slope = int(self.get_slope_pos()[0])
        file_out = self.versuch + '_CNC.MPF'
        if self.txt_outfolder.text() == '':
            fname = os.path.join(self.dir, file_out)
        else:  
            fname = os.path.join(self.txt_outfolder.text())
        
        

        with open(fname, 'w') as f:     # direkte Erzeugung einer .MPF Datei
            # f.write("IB_CURVE:\nG0 SL SLs)\n")    # Schreiben des Headers
            f.write("G1 G91 G64 SL _SLo)\n")

            cnc = self.cnc

            for i in range(len(cnc)):
                '''
                zeilenweises Auslesen des Arrays und Erzeugung des CNC-Syntaxes
                - 0: Spalte Winkel
                - 2: Spalte Linsenstrom
                - 3: Spalte Vorschubgeschwindigkeit
                - 6: Spalte Strahlstrom
                '''

                ang = int(self.par_ang.text())
                d_ang = float(self.par_d_ang.text())
                slp_in = int(self.par_slp_in.text())

                _A = f'A={d_ang}'
                _SQ = f' SQ ({cnc["IB-korr2"].iloc[i]} + _SQ_off))'
                _SL = f' SL {cnc["SL"].iloc[i]})'
                _Fs = f' Fms {cnc["vs"].iloc[i]}'
                if i == 0:
                    f.write(_A + _SQ + _SL + _Fs + "\n")
                elif cnc['Winkel'].iloc[i] <= slp_in:
                    f.write(_A + _SQ + _SL + "\n")
                elif cnc['Winkel'].iloc[i] >= (ang - pos_slope):
                    f.write(_A + _SQ + _SL + _Fs + "\n")
                else:
                    f.write(_A + _SQ + _Fs + "\n")

    def enable_butt(self):
        self.par_n_med.setEnabled(False)
        self.par_n_sav.setEnabled(True)
        self.par_p_sav.setEnabled(True)

    def disable_butt(self):
        self.par_n_med.setEnabled(True)
        self.par_n_sav.setEnabled(False)
        self.par_p_sav.setEnabled(False)

    def save_var_state(self):
        self.var_save = [
            self.par_i0.text(),
            self.par_di.text(),
            self.par_t_soll.text(),
            self.par_d_ang.text(),
            self.par_slp_in.text(),
            self.par_slp_out_0_1.text(),
            self.par_slp_out_0_2.text(),
            self.par_slp_out_0_3.text(),
            self.par_slp_out_1_1.text(),
            self.par_slp_out_1_2.text(),
            self.par_slp_out_1_3.text(),
            self.par_il.text(),
            self.par_slp_out_foc.text(),
            self.par_vs.text(),
            self.par_slp_out_vs.text(),
            self.par_ang.text(),
            self.filt_met_med.isChecked(),
            self.filt_met_sav.isChecked(),
            self.par_n_med.text(),
            self.par_n_sav.text(),
            self.par_p_sav.text()        
            ]

    def get_act_state(self):
        self.var_act = [
            self.par_i0.text(),
            self.par_di.text(),
            self.par_t_soll.text(),
            self.par_d_ang.text(),
            self.par_slp_in.text(),
            self.par_slp_out_0_1.text(),
            self.par_slp_out_0_2.text(),
            self.par_slp_out_0_3.text(),
            self.par_slp_out_1_1.text(),
            self.par_slp_out_1_2.text(),
            self.par_slp_out_1_3.text(),
            self.par_il.text(),
            self.par_slp_out_foc.text(),
            self.par_vs.text(),
            self.par_slp_out_vs.text(),
            self.par_ang.text(),
            self.filt_met_med.isChecked(),
            self.filt_met_sav.isChecked(),
            self.par_n_med.text(),
            self.par_n_sav.text(),
            self.par_p_sav.text()        
            ]

    def replot(self):
        self.get_act_state()
        a = self.var_act
        b = self.var_save

        if not np.array_equal(a[:14], b[:14]):
            self.calc_sq()
            self.filter_data()
            self.create_cnc()
            self.replot_all()
            print('Der Strahlstrom muss neu berechnet werden!')
        elif a[15] != b[15]:
            self.init_data()
            self.calc_sq()
            self.filter_data()
            self.create_cnc()
            self.replot_all()
            print('Alles muss neu berechnet werden!')
        elif not np.array_equal(a[16:], b[16:]):
            self.filter_data()
            self.create_cnc()
            self.replot_all()
            print('Die Ergebnisse müssen neu gefiltert werden!')

    def replot_all(self):
        w = self.graphicsView.canvas
        x_lim = self.ax1.get_xlim()
        y1_lim = self.ax1.get_ylim()
        y2_lim = self.ax2.get_ylim()
        w.figure.clf()

        x = self.df['Winkel'].values
        y1 = self.df['Temperatur'].values
        y2_1 = self.df['IB'].values
        y2_2 = self.df['IB-korr'].values
        x2 = self.cnc['Winkel'].values
        y2_3 = self.cnc['IB-korr2'].values

        ang = int(self.par_ang.text())
        t_soll = int(self.par_t_soll.text())

        self.ax1 = w.figure.add_subplot(111)
        self.ax2 = self.ax1.twinx()

        self.ax1.plot(x, y1, 'r-')
        self.ax1.plot((0, ang), (t_soll, t_soll), 'b:')

        # self.ax1.set_title('Strahlstromkurve')
        self.ax1.set_xlabel('Winkel [°]')
        self.ax1.axis([x_lim[0], x_lim[1], y1_lim[0], y1_lim[1]])

        self.ax1.set_xlabel('Winkel [°]')
        # Make the y-axis label, ticks and tick labels match the line color.
        self.ax1.set_ylabel('Temperatur [°C]', color='r')
        self.ax1.tick_params('y', colors='r')

        text = "T_Soll={:.0f} °C".format(t_soll)
        self.ax1.annotate(text, xy=(0.1, t_soll),
                          xytext=(15, (t_soll+(t_soll*0.03))), color='b')

        self.ax2.plot(x, y2_1, 'g-', alpha=0.3, label='Rohdaten')
        self.ax2.plot(x, y2_2, 'g-', label='Savgol-Filter')
        self.ax2.plot(x2, y2_3, 'k--', label='CNC')
        self.ax2.set_ylabel('Strahlstrom [mA]', color='g')
        self.ax2.set_ylim([y2_lim[0], y2_lim[1]])
        self.ax2.tick_params('y', colors='g')
        self.ax2.legend(loc=(.05, .05))

        self.save_var_state()

        w.draw_idle()

    def plot_data(self):
        w = self.graphicsView.canvas

        w.figure.clf()

        x = self.df['Winkel'].values
        y1 = self.df['Temperatur'].values
        y2_1 = self.df['IB'].values
        y2_2 = self.df['IB-korr'].values
        x2 = self.cnc['Winkel'].values
        y2_3 = self.cnc['IB-korr2'].values

        ang = int(self.par_ang.text())
        t_soll = int(self.par_t_soll.text())

        self.ax1 = w.figure.add_subplot(111)
        self.ax2 = self.ax1.twinx()

        self.ax1.plot(x, y1, 'r-')
        self.ax1.plot((0, ang), (t_soll, t_soll), 'b:')

        self.ax1.set_xlabel('Winkel [°]')
        self.ax1.axis([0, ang, 0, t_soll*1.1])

        self.ax1.set_xlabel('Winkel [°]')
        self.ax1.set_ylabel('Temperatur [°C]', color='r')
        self.ax1.tick_params('y', colors='r')

        text = "T_Soll={:.0f} °C".format(t_soll)
        self.ax1.annotate(text, xy=(0.1, t_soll),
                          xytext=(15, (t_soll+(t_soll*0.03))), color='b')

        self.ax2.plot(x, y2_1, 'g-', alpha=0.3, label='Rohdaten')
        self.ax2.plot(x, y2_2, 'g-', label='Savgol-Filter')
        self.ax2.plot(x2, y2_3, 'k--', label='CNC')
        self.ax2.set_ylabel('Strahlstrom [mA]', color='g')
        self.ax2.set_ylim([0, 90])
        self.ax2.tick_params('y', colors='g')
        self.ax2.legend(loc=(.05, .05))

        self.save_var_state()

        w.draw_idle()

    def print_data(self):
        # Hilfsfunktion, um den aktuellen Wert des Dataframes und der CNC-Ausgabe auszugeben
        print(d)



def main():
    if not QApplication.instance():
        app = QApplication(sys.argv)
    else:
        app = QApplication.instance()
    form = MyApp()             # We set the form to be our ExampleApp (bsp1)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app

if __name__ == '__main__':
    main()