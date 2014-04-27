from sys import argv, exit
from PyQt4.QtGui import QApplication
import MainWindow

if __name__ == '__main__':
    mainapp = QApplication(argv)
    mainForm = MainWindow.MainWindow()
    mainForm.show()
    exit(mainapp.exec_())
    
