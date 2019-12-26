# -*- coding: utf-8 -*-

import os
import sys
import re
import csv

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_AppGenEpi(object):
    ### *******************************************************
    def __init__(self):
        self.out = os.path.expanduser("~")
        self.cmd_exec = "GenEpi"
        self.cmd_args = ["-g", "example", "-p", "example", "-o", "./"]
        self.process = QtCore.QProcess()
        self.process.readyRead.connect(self.dataReady)
        self.process.started.connect(lambda: self.bt_run.setEnabled(False))
        self.process.finished.connect(self.presentResult)
    ### *******************************************************

    def setupUi(self, AppGenEpi):
        AppGenEpi.setObjectName("AppGenEpi")
        AppGenEpi.resize(1146, 606)
        self.centralwidget = QtWidgets.QWidget(AppGenEpi)
        self.centralwidget.setObjectName("centralwidget")
        self.gb_io = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_io.setGeometry(QtCore.QRect(2, 0, 490, 108))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_io.setFont(font)
        self.gb_io.setObjectName("gb_io")
        self.tx_geno = QtWidgets.QPlainTextEdit(self.gb_io)
        self.tx_geno.setGeometry(QtCore.QRect(165, 25, 240, 24))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.tx_geno.setFont(font)
        self.tx_geno.setObjectName("tx_geno")
        self.lb_geno = QtWidgets.QLabel(self.gb_io)
        self.lb_geno.setGeometry(QtCore.QRect(10, 25, 140, 24))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_geno.setFont(font)
        self.lb_geno.setObjectName("lb_geno")
        self.bt_geno = QtWidgets.QPushButton(self.gb_io)
        self.bt_geno.setGeometry(QtCore.QRect(402, 21, 90, 36))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.bt_geno.setFont(font)
        self.bt_geno.setAutoDefault(False)
        self.bt_geno.setObjectName("bt_geno")
        self.tx_pheno = QtWidgets.QPlainTextEdit(self.gb_io)
        self.tx_pheno.setGeometry(QtCore.QRect(165, 52, 240, 24))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.tx_pheno.setFont(font)
        self.tx_pheno.setObjectName("tx_pheno")
        self.lb_pheno = QtWidgets.QLabel(self.gb_io)
        self.lb_pheno.setGeometry(QtCore.QRect(10, 52, 140, 24))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_pheno.setFont(font)
        self.lb_pheno.setObjectName("lb_pheno")
        self.bt_pheno = QtWidgets.QPushButton(self.gb_io)
        self.bt_pheno.setGeometry(QtCore.QRect(402, 48, 90, 36))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.bt_pheno.setFont(font)
        self.bt_pheno.setObjectName("bt_pheno")
        self.bt_out = QtWidgets.QPushButton(self.gb_io)
        self.bt_out.setGeometry(QtCore.QRect(402, 75, 90, 36))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.bt_out.setFont(font)
        self.bt_out.setObjectName("bt_out")
        self.lb_out = QtWidgets.QLabel(self.gb_io)
        self.lb_out.setGeometry(QtCore.QRect(10, 79, 140, 24))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_out.setFont(font)
        self.lb_out.setObjectName("lb_out")
        self.tx_out = QtWidgets.QPlainTextEdit(self.gb_io)
        self.tx_out.setGeometry(QtCore.QRect(165, 79, 240, 24))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.tx_out.setFont(font)
        self.tx_out.setPlainText("")
        self.tx_out.setObjectName("tx_out")
        self.gb_model = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_model.setGeometry(QtCore.QRect(2, 108, 490, 110))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_model.setFont(font)
        self.gb_model.setObjectName("gb_model")
        self.lb_model = QtWidgets.QLabel(self.gb_model)
        self.lb_model.setGeometry(QtCore.QRect(10, 25, 80, 28))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_model.setFont(font)
        self.lb_model.setObjectName("lb_model")
        self.lb_fold = QtWidgets.QLabel(self.gb_model)
        self.lb_fold.setGeometry(QtCore.QRect(10, 52, 80, 28))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_fold.setFont(font)
        self.lb_fold.setObjectName("lb_fold")
        self.lb_thread = QtWidgets.QLabel(self.gb_model)
        self.lb_thread.setGeometry(QtCore.QRect(10, 79, 80, 28))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_thread.setFont(font)
        self.lb_thread.setObjectName("lb_thread")
        self.cb_model = QtWidgets.QComboBox(self.gb_model)
        self.cb_model.setGeometry(QtCore.QRect(162, 26, 327, 26))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cb_model.setFont(font)
        self.cb_model.setObjectName("cb_model")
        self.cb_model.addItem("")
        self.cb_model.addItem("")
        self.cb_fold = QtWidgets.QComboBox(self.gb_model)
        self.cb_fold.setGeometry(QtCore.QRect(162, 53, 327, 26))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cb_fold.setFont(font)
        self.cb_fold.setObjectName("cb_fold")
        self.cb_fold.addItem("")
        self.cb_fold.addItem("")
        self.cb_fold.addItem("")
        self.cb_thread = QtWidgets.QComboBox(self.gb_model)
        self.cb_thread.setGeometry(QtCore.QRect(162, 80, 327, 26))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cb_thread.setFont(font)
        self.cb_thread.setObjectName("cb_thread")
        self.cb_thread.addItem("")
        self.cb_thread.addItem("")
        self.cb_thread.addItem("")
        self.cb_thread.addItem("")
        self.cb_thread.addItem("")
        self.gb_option = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_option.setGeometry(QtCore.QRect(2, 218, 490, 110))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_option.setFont(font)
        self.gb_option.setObjectName("gb_option")
        self.cb_db = QtWidgets.QComboBox(self.gb_option)
        self.cb_db.setGeometry(QtCore.QRect(162, 53, 327, 26))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cb_db.setFont(font)
        self.cb_db.setObjectName("cb_db")
        self.cb_db.addItem("")
        self.cb_db.addItem("")
        self.cb_db.addItem("")
        self.cb_db.addItem("")
        self.ck_region = QtWidgets.QCheckBox(self.gb_option)
        self.ck_region.setGeometry(QtCore.QRect(6, 28, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.ck_region.setFont(font)
        self.ck_region.setObjectName("ck_region")
        self.bt_region = QtWidgets.QPushButton(self.gb_option)
        self.bt_region.setGeometry(QtCore.QRect(402, 21, 90, 36))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.bt_region.setFont(font)
        self.bt_region.setObjectName("bt_region")
        self.tx_region = QtWidgets.QPlainTextEdit(self.gb_option)
        self.tx_region.setGeometry(QtCore.QRect(165, 25, 240, 24))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.tx_region.setFont(font)
        self.tx_region.setPlainText("")
        self.tx_region.setObjectName("tx_region")
        self.ck_db = QtWidgets.QCheckBox(self.gb_option)
        self.ck_db.setGeometry(QtCore.QRect(6, 55, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.ck_db.setFont(font)
        self.ck_db.setObjectName("ck_db")
        self.ck_ld = QtWidgets.QCheckBox(self.gb_option)
        self.ck_ld.setGeometry(QtCore.QRect(6, 82, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.ck_ld.setFont(font)
        self.ck_ld.setObjectName("ck_ld")
        self.lb_d = QtWidgets.QLabel(self.gb_option)
        self.lb_d.setGeometry(QtCore.QRect(165, 78, 80, 28))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_d.setFont(font)
        self.lb_d.setObjectName("lb_d")
        self.sb_r = QtWidgets.QDoubleSpinBox(self.gb_option)
        self.sb_r.setGeometry(QtCore.QRect(405, 80, 82, 24))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.sb_r.setFont(font)
        self.sb_r.setMinimum(0.8)
        self.sb_r.setMaximum(1.0)
        self.sb_r.setSingleStep(0.01)
        self.sb_r.setProperty("value", 0.9)
        self.sb_r.setObjectName("sb_r")
        self.sb_d = QtWidgets.QDoubleSpinBox(self.gb_option)
        self.sb_d.setGeometry(QtCore.QRect(240, 80, 82, 24))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.sb_d.setFont(font)
        self.sb_d.setMinimum(0.8)
        self.sb_d.setMaximum(1.0)
        self.sb_d.setSingleStep(0.01)
        self.sb_d.setProperty("value", 0.9)
        self.sb_d.setObjectName("sb_d")
        self.lb_r = QtWidgets.QLabel(self.gb_option)
        self.lb_r.setGeometry(QtCore.QRect(330, 78, 80, 28))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.lb_r.setFont(font)
        self.lb_r.setObjectName("lb_r")
        self.gb_console = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_console.setGeometry(QtCore.QRect(2, 328, 490, 225))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_console.setFont(font)
        self.gb_console.setObjectName("gb_console")
        self.tx_console = QtWidgets.QTextEdit(self.gb_console)
        self.tx_console.setGeometry(QtCore.QRect(2, 26, 485, 195))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.tx_console.setFont(font)
        self.tx_console.setAutoFillBackground(False)
        self.tx_console.setStyleSheet("background-color: rgb(0, 0, 0); border-color: rgb(0, 0, 0); color: white")
        self.tx_console.setObjectName("tx_console")
        self.gb_hist = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_hist.setGeometry(QtCore.QRect(498, 0, 320, 328))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_hist.setFont(font)
        self.gb_hist.setObjectName("gb_hist")
        self.lb_hist = QtWidgets.QLabel(self.gb_hist)
        self.lb_hist.setGeometry(QtCore.QRect(5, 25, 310, 300))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.lb_hist.setFont(font)
        self.lb_hist.setStyleSheet("background-color: rgb(220, 220, 220);")
        self.lb_hist.setAlignment(QtCore.Qt.AlignCenter)
        self.lb_hist.setObjectName("lb_hist")
        self.gb_prs = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_prs.setGeometry(QtCore.QRect(824, 0, 320, 328))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_prs.setFont(font)
        self.gb_prs.setObjectName("gb_prs")
        self.lb_prs = QtWidgets.QLabel(self.gb_prs)
        self.lb_prs.setGeometry(QtCore.QRect(5, 25, 310, 300))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.lb_prs.setFont(font)
        self.lb_prs.setStyleSheet("background-color: rgb(220, 220, 220);")
        self.lb_prs.setAlignment(QtCore.Qt.AlignCenter)
        self.lb_prs.setObjectName("lb_prs")
        self.gb_result = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_result.setGeometry(QtCore.QRect(498, 328, 646, 225))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.gb_result.setFont(font)
        self.gb_result.setObjectName("gb_result")
        self.tv_result = QtWidgets.QTableView(self.gb_result)
        self.tv_result.setGeometry(QtCore.QRect(5, 26, 636, 195))
        self.tv_result.setStyleSheet("background-color: rgb(220, 220, 220); border: none;")
        self.tv_result.setObjectName("tv_result")
        self.bt_run = QtWidgets.QPushButton(self.centralwidget)
        self.bt_run.setGeometry(QtCore.QRect(-2, 552, 1146, 36))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(False)
        font.setWeight(50)
        self.bt_run.setFont(font)
        self.bt_run.setObjectName("bt_run")
        AppGenEpi.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(AppGenEpi)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1146, 22))
        self.menubar.setObjectName("menubar")
        AppGenEpi.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(AppGenEpi)
        self.statusbar.setObjectName("statusbar")
        AppGenEpi.setStatusBar(self.statusbar)

        self.retranslateUi(AppGenEpi)

        ### *******************define signals**********************
        self.tx_consoleInit()
        self.bt_geno.clicked.connect(self.btn_chooseGenoFile)
        self.bt_pheno.clicked.connect(self.btn_choosePhenoFile)
        self.bt_out.clicked.connect(self.btn_chooseOutDir)
        self.bt_region.clicked.connect(self.btn_chooseRegionFile)
        self.bt_run.clicked.connect(self.btn_callProgram)
        ### *******************************************************

        QtCore.QMetaObject.connectSlotsByName(AppGenEpi)

    def retranslateUi(self, AppGenEpi):
        _translate = QtCore.QCoreApplication.translate
        AppGenEpi.setWindowTitle(_translate("AppGenEpi", "AppGenEpi"))
        self.gb_io.setTitle(_translate("AppGenEpi", "Setting I/O"))
        self.tx_geno.setPlainText(_translate("AppGenEpi", "example"))
        self.lb_geno.setText(_translate("AppGenEpi", "Input genotype file:"))
        self.bt_geno.setText(_translate("AppGenEpi", "Browse"))
        self.tx_pheno.setPlainText(_translate("AppGenEpi", "example"))
        self.lb_pheno.setText(_translate("AppGenEpi", "Input phenotype file:"))
        self.bt_pheno.setText(_translate("AppGenEpi", "Browse"))
        self.bt_out.setText(_translate("AppGenEpi", "Browse"))
        self.lb_out.setText(_translate("AppGenEpi", "Output directory:"))
        self.gb_model.setTitle(_translate("AppGenEpi", "Setting Model Args"))
        self.lb_model.setText(_translate("AppGenEpi", "Model type:"))
        self.lb_fold.setText(_translate("AppGenEpi", "K-fold CV:"))
        self.lb_thread.setText(_translate("AppGenEpi", "Threads:"))
        self.cb_model.setItemText(0, _translate("AppGenEpi", "Classification"))
        #self.cb_model.setItemText(1, _translate("AppGenEpi", "Regression"))
        self.cb_fold.setItemText(0, _translate("AppGenEpi", "2"))
        self.cb_fold.setItemText(1, _translate("AppGenEpi", "4"))
        self.cb_fold.setItemText(2, _translate("AppGenEpi", "8"))
        self.cb_thread.setItemText(0, _translate("AppGenEpi", "2"))
        self.cb_thread.setItemText(1, _translate("AppGenEpi", "4"))
        self.cb_thread.setItemText(2, _translate("AppGenEpi", "8"))
        self.cb_thread.setItemText(3, _translate("AppGenEpi", "16"))
        self.cb_thread.setItemText(4, _translate("AppGenEpi", "32"))
        self.gb_option.setTitle(_translate("AppGenEpi", "Setting Optional Args"))
        self.cb_db.setCurrentText(_translate("AppGenEpi", "hg18"))
        self.cb_db.setItemText(0, _translate("AppGenEpi", "hg18"))
        self.cb_db.setItemText(1, _translate("AppGenEpi", "hg19"))
        self.cb_db.setItemText(2, _translate("AppGenEpi", "hg38"))
        self.cb_db.setItemText(3, _translate("AppGenEpi", "mm10"))
        self.ck_region.setText(_translate("AppGenEpi", "Self-defined regions:"))
        self.bt_region.setText(_translate("AppGenEpi", "Browse"))
        self.ck_db.setText(_translate("AppGenEpi", "Update database:"))
        self.ck_ld.setText(_translate("AppGenEpi", "LD prune:"))
        self.lb_d.setText(_translate("AppGenEpi", "d-prime:"))
        self.lb_r.setText(_translate("AppGenEpi", "r-square"))
        self.gb_console.setTitle(_translate("AppGenEpi", "Console"))
        self.gb_hist.setTitle(_translate("AppGenEpi", "Histogram of samples"))
        self.lb_hist.setText(_translate("AppGenEpi", "Not Ready"))
        self.gb_prs.setTitle(_translate("AppGenEpi", "PRS to prevalence plot"))
        self.lb_prs.setText(_translate("AppGenEpi", "Not Ready"))
        self.gb_result.setTitle(_translate("AppGenEpi", "Result table"))
        self.bt_run.setText(_translate("AppGenEpi", "Run GenEpi !!"))

    def tx_consoleInit(self):
        ### initial console
        str_help = "usage: GenEpi [-h] -g G -p P [-s S] [-o O] [-m {c,r}] [-k K] [-t T]\r\n"
        str_help += "\t[--updatedb] [-b {hg19,hg38}] [--compressld] [-d D] [-r R]\r\n"
        str_help += "\r\n"
        str_help += "\r\n"
        str_help += "\r\n"
        str_help += "optional arguments:\r\n"
        str_help += "  -h, --help\t\tshow this help message and exit\r\n"
        str_help += "  -g G\t\tfilename of the input .gen file\r\n"
        str_help += "  -p P\t\tfilename of the input phenotype\r\n"
        str_help += "  -s S\t\tself-defined genome regions\r\n"
        str_help += "  -o O\t\toutput file path\r\n"
        str_help += "  -m {c,r}\t\tchoose model type: c for classification; r for regression\r\n"
        str_help += "  -k K\t\tk of k-fold cross validation\r\n"
        str_help += "  -t T\t\tnumber of threads\r\n"
        str_help += "\r\n"
        str_help += "\r\n"
        str_help += "update UCSC database:\r\n"
        str_help += "  --updatedb\t\tenable this function\r\n"
        str_help += "  -b {hg19,hg38}\thuman genome build\r\n"
        str_help += "\r\n"
        str_help += "\r\n"
        str_help += "compress data by LD block:\r\n"
        str_help += "  --compressld\tenable this function\r\n"
        str_help += "  -d D\t\tthreshold for compression: D prime\r\n"
        str_help += "  -r R\t\tthreshold for compression: R square\r\n"
        str_help += "\r\n"
        str_help += "\r\n"
        self.tx_console.setPlainText(str_help)

    def btn_chooseGenoFile(self):
        ### open file dialog by specific extension
        fileName_choose, filetype = QtWidgets.QFileDialog.getOpenFileName(self.centralwidget, "Choose a file", self.out, "Genotype Files (*.gen);;All Files (*)")
        ### if no file choosed than break
        if fileName_choose == "":
            return
        ### set filename to plain text
        self.tx_geno.setPlainText(fileName_choose)

    def btn_choosePhenoFile(self):
        ### open file dialog by specific extension
        fileName_choose, filetype = QtWidgets.QFileDialog.getOpenFileName(self.centralwidget, "Choose a file", self.out, "Phenotype Files (*.csv);;All Files (*)")
        ### if no file choosed than break
        if fileName_choose == "":
            return
        ### set filename to plain text
        self.tx_pheno.setPlainText(fileName_choose)
    
    def btn_chooseOutDir(self):
        ### open directory dialog by specific extension
        dir_choose = QtWidgets.QFileDialog.getExistingDirectory(self.centralwidget, "Choose a directory", self.out)
        ### if no directory choosed than break
        if dir_choose == "":
            return
        ### set filename to plain text
        self.tx_out.setPlainText(dir_choose)
    
    def btn_chooseRegionFile(self):
        ### open file dialog by specific extension
        fileName_choose, filetype = QtWidgets.QFileDialog.getOpenFileName(self.centralwidget, "Choose a file", self.out, "Self-defined Regions Files (*.txt);;All Files (*)")
        ### if no file choosed than break
        if fileName_choose == "":
            return
        ### set filename to plain text
        self.tx_region.setPlainText(fileName_choose)

    def dataReady(self):
        ### move cursor to the end of textedit
        cursor = self.tx_console.textCursor()
        cursor.movePosition(cursor.End)
        ### read line from process than add line to current position
        for line in self.process.readAll().data().decode().split('\r\n'):
            cursor.insertText(line)
        ### move cursor to rhe end again
        self.tx_console.ensureCursorVisible()
        self.tx_console.moveCursor(cursor.End)
    
    def btn_callProgram(self):
        ### concatenate command
        list_command = []

        list_command.append("-g")
        list_command.append(self.tx_geno.toPlainText())
    
        list_command.append("-p")
        list_command.append(self.tx_pheno.toPlainText())

        list_command.append("-o")    
        if self.tx_out.toPlainText() == "":
            list_command.append(self.out)
        else:
            list_command.append(self.tx_out.toPlainText())

        list_command.append("-m")
        if self.cb_model.currentText() == "Classification":
            list_command.append("c")
        else:
            list_command.append("r")
        
        list_command.append("-k")
        list_command.append(self.cb_fold.currentText())
    
        list_command.append("-t")
        list_command.append(self.cb_thread.currentText())

        if self.ck_region.isChecked():
            list_command.append("-s")
            list_command.append(self.tx_region.toPlainText())
        
        if self.ck_db.isChecked():
            list_command.append("--updatedb")
            list_command.append("-b")
            list_command.append(self.cb_db.currentText())
        
        if self.ck_ld.isChecked():
            list_command.append("--compressld")
            list_command.append("-d")
            list_command.append(str(self.sb_d.value()))
            list_command.append("-r")
            list_command.append(str(self.sb_r.value()))

        if self.tx_out.toPlainText() != "":
            self.out = self.tx_out.toPlainText()
        
        ### move cursor to the end of console
        cursor = self.tx_console.textCursor()
        cursor.movePosition(cursor.End)
        ### print command to console
        cursor.insertText("\r\n")
        cursor.insertText("GenEpi command: GenEpi " + " ".join(list_command))
        cursor.insertText("\r\n")
        cursor.insertText("\r\n")
        cursor.insertText("GenEpi analysis started")
        cursor.insertText("\r\n")
        cursor.insertText("\r\n")
        ### move cursor to rhe end again
        self.tx_console.ensureCursorVisible()
        self.tx_console.moveCursor(cursor.End)

        ### start process
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.start("GenEpi", list_command)
    
    def presentResult(self):
        ### load histogram
        pm_hist = QtGui.QPixmap(os.path.join(self.out, "crossGeneResult", "GenEpi_PGS.png"))
        if pm_hist.isNull():
            QtWidgets.QMessageBox.information(self.centralwidget, "Image Viewer", "Cannot load %s." % os.path.join(self.out, "crossGeneResult", "GenEpi_PGS.png"))
            return
        
        ### plot histogram
        self.lb_hist.setText("")
        self.lb_hist.setPixmap(pm_hist)
        self.lb_hist.setScaledContents(True)
        self.lb_hist.setAlignment(QtCore.Qt.AlignCenter)
        self.lb_hist.show()

        ### load prevalence
        pm_prs = QtGui.QPixmap(os.path.join(self.out, "crossGeneResult", "GenEpi_Prevalence.png"))
        if pm_prs.isNull():
            QtWidgets.QMessageBox.information(self.centralwidget, "Image Viewer", "Cannot load %s." % os.path.join(self.out, "crossGeneResult", "GenEpi_Prevalence.png"))
            return
        
        ### plot prevalence
        self.lb_prs.setText("")
        self.lb_prs.setPixmap(pm_prs)
        self.lb_prs.setScaledContents(True)
        self.lb_prs.setAlignment(QtCore.Qt.AlignCenter)
        self.lb_prs.show()

        ### show table
        self.model = QtGui.QStandardItemModel(self.centralwidget)
        self.tv_result.setModel(self.model)
        with open(os.path.join(self.out, "crossGeneResult", "Result.csv"), "r") as file_inputFile:
            list_header = file_inputFile.readline().strip().split(",")
            self.model.setHorizontalHeaderLabels(list_header)
            for line in file_inputFile: 
                list_line = line.strip().split(",")
                list_line = [QtGui.QStandardItem(field) for field in list_line]
                self.model.appendRow(list_line)
        self.tv_result.setSortingEnabled(True)
        self.tv_result.setStyleSheet("border: none;")

        ### move cursor to the end of console
        cursor = self.tx_console.textCursor()
        cursor.movePosition(cursor.End)
        ### print finish message to console
        cursor.insertText("\r\n")
        cursor.insertText("Successful finished GenEpi analysis.")
        cursor.insertText("\r\n")
        cursor.insertText("\r\n")
        ### move cursor to rhe end again
        self.tx_console.ensureCursorVisible()
        self.tx_console.moveCursor(cursor.End)

        ### release the run button
        self.bt_run.setEnabled(True)

if __name__ == '__main__':  
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_AppGenEpi()

    ui.setupUi(MainWindow) 
    MainWindow.show()
    sys.exit(app.exec_()) 
    
