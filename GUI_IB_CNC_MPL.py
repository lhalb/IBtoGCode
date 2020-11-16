# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 14:13:17 2018

@author: halbauer
"""
import os

from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
import sys  # We need sys so that we can pass argv to QApplication
import numpy as np
import axesMPL
from Lib import IBtoGCode_Lib
from Lib import IBtoGCode_Helper
from Lib import dataObject as dOj


class MyApp(QMainWindow, axesMPL.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)  # Laden der UI-Datei
        self.setup_triggers()  # lade Events für Buttons
        self.statusBar().showMessage('Programm erfolgreich geladen.')
        self.df = None
        self.var_save = None
        self.data = dOj.DataFile()

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
        if IBtoGCode_Helper.test_if_string_helper(self.data):
            self.statusBar().showMessage('Datei enthält ungültige Zeichen!')
            # self.status_line.setStyleSheet("color: rgb(181, 18, 62);")
            self.but_trim_csv.setEnabled(True)
        else:
            self.statusBar().showMessage('Datei erfolgreich geladen')
            # self.status_line.setStyleSheet("color: rgb(0, 0, 0);")
            self.but_init_data.setEnabled(True)
            self.but_QuickPlot.setEnabled(True)

    def trim_csv(self):
        IBtoGCode_Helper.trim_csv_helper(self.data)
        self.proc_path()
        self.test_if_string()

    def init_data(self):
        df = IBtoGCode_Helper.init_data_helper(self.data.pfad, int(self.par_ang.text()))
        self.data.set_data_frame(df)

        self.but_calc_ib.setEnabled(True)
        self.statusBar().showMessage('Datenmatrix erfolgreich erzeugt. Weiter mit CALC IB!', 2000)
        return df

    def out_to_np(self):
        laengen = [self.par_slp_out_0_1.text(), self.par_slp_out_0_2.text(), self.par_slp_out_0_3.text()]
        werte = [self.par_slp_out_1_1.text(), self.par_slp_out_1_2.text(), self.par_slp_out_1_3.text()]
        return IBtoGCode_Helper.out_to_np_helper(laengen, werte)

    def calc_sq(self):
        df = getattr(self.data, "df")

        # Übernehmen der Einträge aus der GUI
        i_0 = int(self.par_i0.text())
        foc_of = int(self.par_il.text())
        slp_in = int(self.par_slp_in.text())
        foc_ruhe = foc_of + int(self.par_slp_out_foc.text())
        v_s_slope = int(self.par_slp_out_vs.text())
        ang = int(self.par_ang.text())
        di = int(self.par_di.text())
        slp_out = self.out_to_np()
        self.data.set_slp_out(slp_out)

        # Daten aus UI für die Funktion get_vs #################################
        par_vs = self.par_vs.text()
        par_vs_gk = self.par_vs_gk.text()
        par_vs_ns = self.par_vs_ns.text()
        par_slp_in = self.par_slp_in.text()
        par_l_gk = self.par_l_gk.text()
        par_l_ns = self.par_l_ns.text()
        par_vs_trans = self.par_vs_trans.text()
        vs = IBtoGCode_Lib.get_vs(df, self.rb_vs_const.isChecked(), self.rb_vs_step.isChecked(),
                              par_vs, par_vs_gk, par_vs_ns, par_slp_in, par_l_gk,
                              par_l_ns, par_vs_trans)

        df = IBtoGCode_Helper.calc_sq_helper(df, i_0, foc_of, foc_ruhe, vs, slp_in, slp_out, ang, v_s_slope, di)

        self.data.set_data_frame(df)
        self.but_filter_data.setEnabled(True)
        self.statusBar().showMessage('Strahlstromkurven berechnet. Weiter mit FILTER DATA!')

    def filter_data(self):
        df = getattr(self.data, "df")
        # --- Filterung der Daten ---
        df = IBtoGCode_Helper.filter_data_helper(self.filt_met_med.isChecked(), self.par_n_med.text(),
                                          self.filt_met_sav.isChecked(), self.par_n_sav.text(), self.par_p_sav.text(),
                                          df)
        self.but_create_cnc.setEnabled(True)
        self.statusBar().showMessage('Daten erfolgreich gefiltert. Weiter mit CALC CNC oder PLOT DATA!')
        return df

    def create_cnc(self):
        ang = float(self.par_ang.text())
        d_ang = float(self.par_d_ang.text())
        cnc = IBtoGCode_Helper.create_cnc_helper(ang, d_ang, getattr(self.data, "df"))
        self.data.set_cnc(cnc)
        self.but_save_data.setEnabled(True)
        return cnc

    # Hilfsfunktionen, die nur in der GUI benötigt werden

    def file_open(self):
        path = QFileDialog.getOpenFileName(self, 'Open File')[0]
        self.data.set_pfad(path)
        if not path:
            return
        self.data.set_directory(os.path.dirname(path))
        self.proc_path()
        self.test_if_string()

    def file_save(self):
        fname = f'{self.txt_versuch.text()}.MPF'
        out_file = QFileDialog.getSaveFileName(self, 'Speicherort wählen', fname, 'CNC-Files (*.MPF);;All Files (*)')[0]
        if not out_file:
            return
        self.txt_outfolder.setText(out_file)

    def proc_path(self):
        file, versuch = IBtoGCode_Helper.proc_path_helper(self.data)
        self.txt_path.setText(self.data.pfad)
        self.txt_file.setText(file)
        self.txt_versuch.setText(versuch)

    def quick_plot(self):
        self.init_data()
        self.calc_sq()
        df = self.filter_data()
        cnc = self.create_cnc()
        self.plot_data(df, cnc)
        self.but_replot_data.setEnabled(True)

    def save_data(self):
        cnc = self.cb_save_cnc.isChecked()
        csv = self.cb_save_csv.isChecked()
        if csv and cnc:
            self.save_cnc(self.data.versuch, self.data.directory, self.data.cnc)
            IBtoGCode_Lib.save_csv(getattr(self.data, "df"), getattr(self.data, "versuch"), getattr(self.data, "directory"))
            self.statusBar().showMessage('CSV und CNC erfolgreich gesichert.', 2000)
        elif csv:
            IBtoGCode_Lib.save_csv(getattr(self.data, "df"), getattr(self.data, "versuch"), getattr(self.data, "directory"))
            self.statusBar().showMessage('CSV erfolgreich gesichert.', 2000)
        elif cnc:
            self.save_cnc(getattr(self.data, "versuch"), getattr(self.data, "directory"), getattr(self.data, "cnc"))
            self.statusBar().showMessage('CNC-Code erfolgreich exportiert.', 2000)
        else:
            self.statusBar().showMessage('Keine Daten gesichert.', 2000)

    def save_cnc(self, versuch, direction, cnc):
        pos_slope = int(IBtoGCode_Lib.get_slope_pos(getattr(self.data, "slp_out"))[0])
        txt_out_folder = self.txt_outfolder.text()
        ang = int(self.par_ang.text())
        d_ang = float(self.par_d_ang.text())
        slp_in = int(self.par_slp_in.text())

        IBtoGCode_Helper.save_cnc_helper(direction, versuch, txt_out_folder, d_ang, ang, slp_in, pos_slope, cnc)

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
        var_act = [
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
        return var_act

    def replot(self):
        a = self.get_act_state()
        b = self.var_save

        if not np.array_equal(a[:14], b[:14]):
            self.calc_sq()
            df = self.filter_data()
            cnc = self.create_cnc()
            self.replot_all(df, cnc)
            print('Der Strahlstrom muss neu berechnet werden!')
        elif a[15] != b[15]:
            self.init_data()
            self.calc_sq()
            df = self.filter_data()
            cnc = self.create_cnc()
            self.replot_all(df, cnc)
            print('Alles muss neu berechnet werden!')
        elif not np.array_equal(a[16:], b[16:]):
            df = self.filter_data()
            cnc = self.create_cnc()
            self.replot_all(df, cnc)
            print('Die Ergebnisse müssen neu gefiltert werden!')

    def replot_all(self, df, cnc):
        w = self.graphicsView.canvas
        x_lim = self.ax1.get_xlim()
        y1_lim = self.ax1.get_ylim()
        y2_lim = self.ax2.get_ylim()
        w.figure.clf()

        x = df['Winkel'].values
        y1 = df['Temperatur'].values
        y2_1 = df['IB'].values
        y2_2 = df['IB-korr'].values
        x2 = cnc['Winkel'].values
        y2_3 = cnc['IB-korr2'].values

        ang = int(self.par_ang.text())
        t_soll = int(self.par_t_soll.text())

        # todo ist die Definition von ax1 und ax2 an self in diesem Teil nötig?
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
                          xytext=(15, (t_soll + (t_soll * 0.03))), color='b')

        self.ax2.plot(x, y2_1, 'g-', alpha=0.3, label='Rohdaten')
        self.ax2.plot(x, y2_2, 'g-', label='Savgol-Filter')
        self.ax2.plot(x2, y2_3, 'k--', label='CNC')
        self.ax2.set_ylabel('Strahlstrom [mA]', color='g')
        self.ax2.set_ylim([y2_lim[0], y2_lim[1]])
        self.ax2.tick_params('y', colors='g')
        self.ax2.legend(loc=(.05, .05))

        self.save_var_state()

        w.draw_idle()

    def plot_data(self, df, cnc):
        w = self.graphicsView.canvas

        w.figure.clf()

        x = df['Winkel'].values
        y1 = df['Temperatur'].values
        y2_1 = df['IB'].values
        y2_2 = df['IB-korr'].values
        x2 = cnc['Winkel'].values
        y2_3 = cnc['IB-korr2'].values

        ang = int(self.par_ang.text())
        t_soll = int(self.par_t_soll.text())

        # todo ist die Definition von ax1 und ax2 an self in diesem Teil nötig?
        self.ax1 = w.figure.add_subplot(111)
        self.ax2 = self.ax1.twinx()

        self.ax1.plot(x, y1, 'r-')
        self.ax1.plot((0, ang), (t_soll, t_soll), 'b:')

        self.ax1.set_xlabel('Winkel [°]')
        self.ax1.axis([0, ang, 0, t_soll * 1.1])

        self.ax1.set_xlabel('Winkel [°]')
        self.ax1.set_ylabel('Temperatur [°C]', color='r')
        self.ax1.tick_params('y', colors='r')

        text = "T_Soll={:.0f} °C".format(t_soll)
        self.ax1.annotate(text, xy=(0.1, t_soll),
                          xytext=(15, (t_soll + (t_soll * 0.03))), color='b')

        self.ax2.plot(x, y2_1, 'g-', alpha=0.3, label='Rohdaten')
        self.ax2.plot(x, y2_2, 'g-', label='Savgol-Filter')
        self.ax2.plot(x2, y2_3, 'k--', label='CNC')
        self.ax2.set_ylabel('Strahlstrom [mA]', color='g')
        self.ax2.set_ylim([0, 90])
        self.ax2.tick_params('y', colors='g')
        self.ax2.legend(loc=(.05, .05))

        self.save_var_state()

        w.draw_idle()

    @staticmethod  # todo diese Funktion hat keinen nutzen
    def print_data():
        # Hilfsfunktion, um den aktuellen Wert des Dataframes und der CNC-Ausgabe auszugeben
        print()


def main():
    if not QApplication.instance():
        app = QApplication(sys.argv)
    else:
        app = QApplication.instance()
    form = MyApp()  # We set the form to be our ExampleApp (bsp1)
    form.show()  # Show the form
    app.exec_()  # and execute the app


if __name__ == '__main__':
    main()
