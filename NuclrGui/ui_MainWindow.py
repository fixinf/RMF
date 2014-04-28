# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow.ui'
#
# Created: Mon Apr 28 18:26:48 2014
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(618, 512)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.treeView = QtGui.QTreeView(self.centralwidget)
        self.treeView.setObjectName(_fromUtf8("treeView"))
        self.horizontalLayout.addWidget(self.treeView)
        self.groupBox = QtGui.QGroupBox(self.centralwidget)
        self.groupBox.setCheckable(False)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.btnSolve = QtGui.QPushButton(self.groupBox)
        self.btnSolve.setObjectName(_fromUtf8("btnSolve"))
        self.gridLayout.addWidget(self.btnSolve, 5, 0, 1, 1)
        self.linec = QtGui.QLineEdit(self.groupBox)
        self.linec.setObjectName(_fromUtf8("linec"))
        self.gridLayout.addWidget(self.linec, 4, 1, 1, 1)
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 4, 0, 1, 1)
        self.lineCs = QtGui.QLineEdit(self.groupBox)
        self.lineCs.setObjectName(_fromUtf8("lineCs"))
        self.gridLayout.addWidget(self.lineCs, 0, 1, 1, 1)
        self.lineCr = QtGui.QLineEdit(self.groupBox)
        self.lineCr.setObjectName(_fromUtf8("lineCr"))
        self.gridLayout.addWidget(self.lineCr, 2, 1, 1, 1)
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 3, 0, 1, 1)
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.groupBox_2 = QtGui.QGroupBox(self.groupBox)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.formLayout = QtGui.QFormLayout(self.groupBox_2)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.rbP = QtGui.QRadioButton(self.groupBox_2)
        self.rbP.setObjectName(_fromUtf8("rbP"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.rbP)
        self.rbFeq = QtGui.QRadioButton(self.groupBox_2)
        self.rbFeq.setObjectName(_fromUtf8("rbFeq"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.rbFeq)
        self.rbNp = QtGui.QRadioButton(self.groupBox_2)
        self.rbNp.setObjectName(_fromUtf8("rbNp"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.rbNp)
        self.rbPconstr = QtGui.QRadioButton(self.groupBox_2)
        self.rbPconstr.setObjectName(_fromUtf8("rbPconstr"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.rbPconstr)
        self.rbM = QtGui.QRadioButton(self.groupBox_2)
        self.rbM.setObjectName(_fromUtf8("rbM"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.rbM)
        self.rbEtao = QtGui.QRadioButton(self.groupBox_2)
        self.rbEtao.setObjectName(_fromUtf8("rbEtao"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.rbEtao)
        self.btnPlot = QtGui.QPushButton(self.groupBox_2)
        self.btnPlot.setObjectName(_fromUtf8("btnPlot"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.btnPlot)
        self.gridLayout.addWidget(self.groupBox_2, 6, 0, 1, 1)
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.lineb = QtGui.QLineEdit(self.groupBox)
        self.lineb.setObjectName(_fromUtf8("lineb"))
        self.gridLayout.addWidget(self.lineb, 3, 1, 1, 1)
        self.LineCo = QtGui.QLineEdit(self.groupBox)
        self.LineCo.setObjectName(_fromUtf8("LineCo"))
        self.gridLayout.addWidget(self.LineCo, 1, 1, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.horizontalLayout.addWidget(self.groupBox)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 618, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.groupBox.setTitle(_translate("MainWindow", "Parameters", None))
        self.btnSolve.setText(_translate("MainWindow", "Solve", None))
        self.label_5.setText(_translate("MainWindow", "c", None))
        self.label_4.setText(_translate("MainWindow", "b", None))
        self.label_3.setText(_translate("MainWindow", "C_r", None))
        self.groupBox_2.setTitle(_translate("MainWindow", "Plotting", None))
        self.rbP.setText(_translate("MainWindow", "Pressure", None))
        self.rbFeq.setText(_translate("MainWindow", "F_eq", None))
        self.rbNp.setText(_translate("MainWindow", "N_p", None))
        self.rbPconstr.setText(_translate("MainWindow", "P Constraint SNM", None))
        self.rbM.setText(_translate("MainWindow", "M(n)", None))
        self.rbEtao.setText(_translate("MainWindow", "eta_omega(f)", None))
        self.btnPlot.setText(_translate("MainWindow", "Plot", None))
        self.label.setText(_translate("MainWindow", "C_s", None))
        self.label_2.setText(_translate("MainWindow", "C_o", None))

