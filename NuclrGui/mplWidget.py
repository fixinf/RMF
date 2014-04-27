
from PyQt4.QtGui import QWidget
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

class mplWidget(QWidget):
    '''
    Matplotlib widget for PyQt4 with some capabilities of data manipulation
    '''


    def __init__(self, parent):
        '''
        Constructor
        '''
        QWidget.__init__(parent)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(parent)
        
        
        
        