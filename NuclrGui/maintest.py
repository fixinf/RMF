from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication, QMessageBox, QInputDialog, QMainWindow
from sys import argv
import sys

def f(x):
    return x*x

print f(1), f(2), f(3)

for x in range(10):
    print f(x)
    
app = QApplication(argv)  

mainwindow = QMainWindow()
mainwindow.show()  
msgbox = QInputDialog()
text, ok = msgbox.getText(mainwindow, 'Enter text', 'Enter text')

if ok:
    print text
    
    
sys.exit(app.exec_())   