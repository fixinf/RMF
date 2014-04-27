from PyQt4 import uic
from os import path
import sys
from PyQt4.QtGui import QMainWindow, QWidget, QSizePolicy
from matplotlib.figure import Figure
with open('MainWindow.ui', 'r') as infile:
    with open('ui_MainWindow.py', 'w') as outfile:
        uic.compileUi(infile, outfile)
import ui_MainWindow

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from DataModel.TreeModel import Node, TreeModel
from DataModel.Database import RMFClasses, RMFParams

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = ui_MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
#       self.main_frame = QWidget()
#         
#         # Create the mpl Figure and FigCanvas objects. 
#         # 5x4 inches, 100 dots-per-inch
#         #
#         self.dpi = 100
#         self.fig = Figure((5.0, 4.0), dpi=self.dpi)
#         self.canvas = FigureCanvas(self.fig)
#         self.canvas.setParent(self.main_frame)
#         self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
#         self.canvas.updateGeometry()
#         
#         # Since we have only one plot, we can use add_axes 
#         # instead of add_subplot, but then the subplot
#         # configuration tool in the navigation toolbar wouldn't
#         # work.
#         #
#         self.axes = self.fig.add_subplot(111)
#         self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
#         self.ui.vlay.addWidget(self.main_frame)
        
        self.database = RMFParams()
        self.classbase = RMFClasses()
        class1 = RMFClasses.create(class__='TEST')
        record1 = RMFParams.create(name='Test 1', massfile='', tabfile='',
                                   constants='', type=class1, class__='TEST', module='')
        
        
        self.rootNode = Node("root")
        self.model = TreeModel(self.rootNode, self.classbase, self.database)
        self.ui.treeView.setModel(self.model)
    
        
        
       
       
        
        