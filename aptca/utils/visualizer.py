from itertools import cycle
import sys
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from aptca.utils import ColorTable
from pyqtgraph.Qt import QtCore, QtGui
import OpenGL.GL as ogl 



def visualize_atom_cloud(pos,GLView,size=1):
    w = GLView
    data = pos.astype(np.float32)
    data[:, 2] = -data[:, 2] # reverse the z axis to show needle right.

    size = size

    sp1 = gl.GLScatterPlotItem(pos=pos, size=size, color=[0.0,0.0,0.0,0.5], pxMode=True)
    sp1.setGLOptions('translucent')
    w.setBackgroundColor('w')
    w.opts['distance'] = 2000
    w.opts['fov'] = 1
    w.addItem(sp1)
    # add axis
    # axis = Custom3DAxis(w,color=(0.2,0.2,0.2,0.6))
    # maxX = max(pos[:,0])
    # maxY = max(pos[:,1])
    # maxZ = max(pos[:,2])
    # minX = min(pos[:,0])
    # minY = min(pos[:,1])
    # minZ = min(pos[:,2])
    # xticks = np.round(np.linspace(minX,maxX,num=10)).tolist()
    # yticks = np.round(np.linspace(minY,maxY,num=10)).tolist()
    # zticks = np.round(np.linspace(minZ,maxZ,num=10)).tolist()
    # axis.setSize(x=maxX,y=maxY,z=maxZ)
    # axis.add_labels()
    # axis.add_tick_values(xticks=xticks,yticks=yticks,zticks=zticks)
    # w.addItem(axis)
    #w.show()


def visualize_clusters(coord, cluster_id, GLView, background):
    """
    Plot and show cluster point clouds.

    !!!!!!!! Coord is byte sensitive, it will plot scatter so wrong if use default '>f4' data type as in our program.
    Change to default float solves the problem.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    """
    w = GLView
    data = coord.astype(np.float32)
    data[:, 2] = -data[:, 2] # reverse the z axis to show needle right.

    if not background:
        selection = cluster_id != -1
    else:
        selection = np.ones(len(cluster_id), dtype=bool)
        selection.fill(True)

    size = 3

    color_table = ColorTable.gen_color_table(style='float')
    color_table = np.delete(color_table, [0], axis=0)
    cycol = cycle(color_table)
    cluster_color = dict()
    unique_cluster_id = np.unique(cluster_id)

    for unique_id in unique_cluster_id:
        if unique_id == -1:
            cluster_color[unique_id] = [0.0, 0.0, 0.0]
        else:
            cluster_color[unique_id] = next(cycol)

    color = np.empty((data.shape[0], 4))
    color[:, 3] = 1.0
    for idx, c_id in enumerate(cluster_id):
        color[idx, 0:3] = cluster_color[c_id]

    sp1 = gl.GLScatterPlotItem(pos=data[selection], size=size, color=color[selection], pxMode=True)
    sp1.setGLOptions('translucent')
    w.setBackgroundColor('w')
    w.opts['distance'] = 2000
    w.opts['fov'] = 1
    w.addItem(sp1)

def visualization(pos, label_lst, background=True):
    """
    Visualization of clusters and RD plot.
    """
    app = pg.Qt.QtGui.QApplication.instance()
    if app == None:
        app = pg.Qt.QtGui.QApplication([])

    w = pg.Qt.QtGui.QWidget()
    layout = pg.Qt.QtGui.QGridLayout()
    w.setLayout(layout)
    #plot = pg.PlotWidget()
    #hist = pg.PlotWidget()
    view = gl.GLViewWidget()

    #plot.sizeHint = pg.QtCore.QSize(800, 600)
    view.sizeHint = lambda: pg.QtCore.QSize(800, 600)
    #view.setSizePolicy(plot.sizePolicy())
    layout.addWidget(view, 0, 0, 1, 1)
    #layout.addWidget(plot, 0, 1)
    #layout.addWidget(hist, 1, 1)

    #leaves = []
    #leaves = TreeNode.retrieve_leaf_nodes(Clusterer.hierarchy_RootNode, leaves)
    #visualize_RD_hist(Clusterer.RD, hist)
    #visualize_RD(leaves, Clusterer.RD[Clusterer.ordered_lst], plot)
    visualize_clusters(pos, label_lst, view, background)

    w.show()
    sys.exit(app.exec_())



class CustomTextItem(gl.GLGraphicsItem.GLGraphicsItem):
    def __init__(self, X, Y, Z, text):
        gl.GLGraphicsItem.GLGraphicsItem.__init__(self)
        self.text = text
        self.X = X
        self.Y = Y
        self.Z = Z

    def setGLViewWidget(self, GLViewWidget):
        self.GLViewWidget = GLViewWidget

    def setText(self, text):
        self.text = text
        self.update()

    def setX(self, X):
        self.X = X
        self.update()

    def setY(self, Y):
        self.Y = Y
        self.update()

    def setZ(self, Z):
        self.Z = Z
        self.update()

    def paint(self):
        self.GLViewWidget.qglColor(QtCore.Qt.black)
        self.GLViewWidget.renderText(self.X, self.Y, self.Z, self.text)


class Custom3DAxis(gl.GLAxisItem):
    """Class defined to extend 'gl.GLAxisItem'."""
    def __init__(self, parent, color=(0,0,0,.6)):
        gl.GLAxisItem.__init__(self)
        self.parent = parent
        self.c = color

    def add_labels(self):
        """Adds axes labels."""
        x,y,z = self.size()
        #X label
        self.xLabel = CustomTextItem(X=x/2, Y=-y/20, Z=-z/20, text="X")
        self.xLabel.setGLViewWidget(self.parent)
        self.parent.addItem(self.xLabel)
        #Y label
        self.yLabel = CustomTextItem(X=-x/20, Y=y/2, Z=-z/20, text="Y")
        self.yLabel.setGLViewWidget(self.parent)
        self.parent.addItem(self.yLabel)
        #Z label
        self.zLabel = CustomTextItem(X=-x/20, Y=-y/20, Z=z/2, text="Z")
        self.zLabel.setGLViewWidget(self.parent)
        self.parent.addItem(self.zLabel)

    def add_tick_values(self, xticks=[], yticks=[], zticks=[]):
        """Adds ticks values."""
        x,y,z = self.size()
        xtpos = np.linspace(0, x, len(xticks))
        ytpos = np.linspace(0, y, len(yticks))
        ztpos = np.linspace(0, z, len(zticks))
        #X label
        for i, xt in enumerate(xticks):
            val = CustomTextItem(X=xtpos[i], Y=-y/20, Z=-z/20, text=str(xt))
            val.setGLViewWidget(self.parent)
            self.parent.addItem(val)
        #Y label
        for i, yt in enumerate(yticks):
            val = CustomTextItem(X=-x/20, Y=ytpos[i], Z=-z/20, text=str(yt))
            val.setGLViewWidget(self.parent)
            self.parent.addItem(val)
        #Z label
        for i, zt in enumerate(zticks):
            val = CustomTextItem(X=-x/20, Y=-y/20, Z=ztpos[i], text=str(zt))
            val.setGLViewWidget(self.parent)
            self.parent.addItem(val)

    def paint(self):
        self.setupGLState()
        if self.antialias:
            ogl.glEnable(ogl.GL_LINE_SMOOTH)
            ogl.glHint(ogl.GL_LINE_SMOOTH_HINT, ogl.GL_NICEST)
        ogl.glBegin(ogl.GL_LINES)

        x,y,z = self.size()
        #Draw Z
        ogl.glColor4f(self.c[0], self.c[1], self.c[2], self.c[3])
        ogl.glVertex3f(0, 0, 0)
        ogl.glVertex3f(0, 0, z)
        #Draw Y
        ogl.glColor4f(self.c[0], self.c[1], self.c[2], self.c[3])
        ogl.glVertex3f(0, 0, 0)
        ogl.glVertex3f(0, y, 0)
        #Draw X
        ogl.glColor4f(self.c[0], self.c[1], self.c[2], self.c[3])
        ogl.glVertex3f(0, 0, 0)
        ogl.glVertex3f(x, 0, 0)
        ogl.glEnd()