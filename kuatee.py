import wx, os
import datetime
from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigureCanvas,
    NavigationToolbar2WxAgg as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from astropy.table import Table, vstack


import wx.lib.mixins.inspection as WIT

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from tools.nightsky import bscTable, orionTable, set_altaz, set_xyz_proj, Projection, getPlanetTable,getEclipticTable
from tools.nightsky import localtime2utctime
# https://wiki.wxpython.org/Getting%20Started#Working_with_Windows
# https://matplotlib.org/stable/gallery/user_interfaces/embedding_in_wx2_sgskip.html#sphx-glr-gallery-user-interfaces-embedding-in-wx2-sgskip-py
class MyFrame(wx.Frame):
    """ We simply derive a new class of Frame. """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=(800,600))
        #self.control = wx.TextCtrl(self, style=wx.TE_MULTILINE)

        self.bscTable = bscTable
        self.az0 = 195
        self.alt0 = 40
        self.fov = 175
        #self.datetime = '2021-02-22 08:33:00.724'
        self.datetime = "{}".format(datetime.datetime.now())
        self.timezone = 8
        #self.datetime = localtime2utctime(timestr,8)
        #print(self.datetime)
        self.lon = '121.5528'
        self.lat = '24.99136'
        self.currentStar=""
        self.projection = Projection.Lambert
        
        self.figure = Figure()
        self.axes = self.figure.add_subplot()
        self.axes.margins(0)
        self.canvas = FigureCanvas(self, -1, self.figure)
        #self.Bind(wx.EVT_KEY_DOWN, self.on_key_press, self)

        self.cid = self.canvas.mpl_connect('pick_event', self.on_pick)
        self.cidkey_press = self.canvas.mpl_connect('key_press_event', self.on_key_press)

        self.canvasDraw(None)


        self.CreateStatusBar() # A Statusbar in the bottom of the window
        self.SetStatusText("This is status bar")

        # Setting up the menu.
        filemenu= wx.Menu()

        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        menuOpen  = filemenu.Append(wx.ID_ANY, "&Open"," open file for edit")
        menuAbout  = filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")
        
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        # Set events.
        self.Bind(wx.EVT_MENU, self.OnOpen, menuOpen)
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)

        self.sizerControls = wx.BoxSizer(wx.VERTICAL)
        self.starValue = wx.StaticText(self, label="天體:{}".format(self.currentStar))
        self.sizerControls.Add(self.starValue,0,wx.EXPAND)

        self.sizerDateTime = wx.BoxSizer(wx.HORIZONTAL)
        #self.datepicker = wx.adv.DatePickerCtrl(self)
        #self.timepicker = wx.adv.TimePickerCtrl(self)
        #self.sizerDateTime.Add(self.datepicker,0,wx.EXPAND)
        #self.sizerDateTime.Add(self.timepicker,0,wx.EXPAND)
        self.datetimeLabel = wx.StaticText(self, label="日期時間:")
        self.datetimeValue = wx.TextCtrl(self, value=self.datetime)
        self.datetimeButton = wx.Button(self, label="update sky")
        self.Bind(wx.EVT_BUTTON, self.onUpdateTime, self.datetimeButton)
        self.sizerDateTime.Add(self.datetimeLabel,0,wx.EXPAND)
        self.sizerDateTime.Add(self.datetimeValue,0,wx.EXPAND)
        self.sizerDateTime.Add(self.datetimeButton,0,wx.EXPAND)
        self.sizerControls.Add(self.sizerDateTime,0,wx.EXPAND)

        self.sizerLon = wx.BoxSizer(wx.HORIZONTAL)
        self.lonLabel = wx.StaticText(self, label="經度:")
        self.lonValue = wx.TextCtrl(self, value=self.lon)
        self.lonButton = wx.Button(self, label="update sky")
        self.Bind(wx.EVT_BUTTON, self.onUpdateLocation, self.lonButton)
        self.sizerLon.Add(self.lonLabel,0)
        self.sizerLon.Add(self.lonValue,0)
        self.sizerLon.Add(self.lonButton,0)
        self.sizerControls.Add(self.sizerLon,0,wx.EXPAND)

        self.sizerLat = wx.BoxSizer(wx.HORIZONTAL)
        self.latLabel = wx.StaticText(self, label="緯度:")
        self.latValue = wx.TextCtrl(self, value=self.lat)
        self.latButton = wx.Button(self, label="update sky")
        self.Bind(wx.EVT_BUTTON, self.onUpdateLocation, self.latButton)
        self.sizerLat.Add(self.latLabel,0)
        self.sizerLat.Add(self.latValue,0)
        self.sizerLat.Add(self.latButton,0)
        self.sizerControls.Add(self.sizerLat,0,wx.EXPAND)

        self.sizerAz = wx.BoxSizer(wx.HORIZONTAL)
        self.azLabel = wx.StaticText(self, label="方位角:")
        self.azValue = wx.StaticText(self, label="{}".format(self.az0), size=(20,10))
        self.azButton = wx.Button(self, label="Add 15")
        self.Bind(wx.EVT_BUTTON, self.onAddAz, self.azButton)
        self.azSubButton = wx.Button(self, label="Sub 15")
        self.Bind(wx.EVT_BUTTON, self.onSubAz, self.azSubButton)
        self.sizerAz.Add(self.azLabel,0)
        self.sizerAz.Add(self.azValue,0)
        self.sizerAz.Add(self.azButton,0)
        self.sizerAz.Add(self.azSubButton,0)
        self.sizerControls.Add(self.sizerAz,0,wx.EXPAND)
        
        self.sizerAlt = wx.BoxSizer(wx.HORIZONTAL)
        self.altLabel = wx.StaticText(self, label="仰角:")
        self.altValue = wx.StaticText(self, label="{}".format(self.alt0), size=(20,10))
        self.altAddButton = wx.Button(self, label="Add 10")
        self.Bind(wx.EVT_BUTTON, self.onAddAlt, self.altAddButton)
        self.altSubButton = wx.Button(self, label="Sub 10")
        self.Bind(wx.EVT_BUTTON, self.onSubAlt, self.altSubButton)
        self.sizerAlt.Add(self.altLabel,0)
        self.sizerAlt.Add(self.altValue,0)
        self.sizerAlt.Add(self.altAddButton,0)
        self.sizerAlt.Add(self.altSubButton,0)
        self.sizerControls.Add(self.sizerAlt,0,wx.EXPAND)
        
        self.sizerFov = wx.BoxSizer(wx.HORIZONTAL)
        self.fovLabel = wx.StaticText(self, label="視角:")
        self.fovValue = wx.StaticText(self, label="{}".format(self.fov), size=(20,10))
        self.fovAddButton = wx.Button(self, label="Add 5")
        self.Bind(wx.EVT_BUTTON, self.onAddFov, self.fovAddButton)
        self.fovSubButton = wx.Button(self, label="Sub 5")
        self.Bind(wx.EVT_BUTTON, self.onSubFov, self.fovSubButton)
        self.sizerFov.Add(self.fovLabel,0)
        self.sizerFov.Add(self.fovValue,0)
        self.sizerFov.Add(self.fovAddButton,0)
        self.sizerFov.Add(self.fovSubButton,0)
        self.sizerControls.Add(self.sizerFov,0,wx.EXPAND)
        
        #self.drawButton = wx.Button(self, label="Draw")
        #self.Bind(wx.EVT_BUTTON, self.canvasDraw, self.drawButton)
        #self.sizerControls.Add(self.drawButton,0,wx.EXPAND)

        self.buttons = []
        for i in range(0, 6):
            self.buttons.append(wx.Button(self, -1, "Button &"+str(i)))
            self.sizerControls.Add(self.buttons[i], 1, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.OnPaint, self.buttons[0])
        # Use some sizers to see layout options
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        #self.sizer.Add(self.control, 1, wx.EXPAND)
        self.sizer.Add(self.canvas, 1, wx.EXPAND)
        self.sizer.Add(self.sizerControls, 0, wx.EXPAND)

        #Layout sizers
        self.SetSizer(self.sizer)
        self.SetAutoLayout(True)
        self.sizer.Fit(self)

        self.Show(True)

    def OnAbout(self, event):
        # A message dialog box with an OK button. wx.OK is a standard ID in wxWidgets.
        dlg = wx.MessageDialog( self, "A small text editor", "About Sample Editor", wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy() # finally destroy it when finished.
    
    def OnExit(self,event):
        self.Close(True)  # Close the frame.

    def OnOpen(self,e):
        """ Open a file"""
        self.dirname = ''
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = open(os.path.join(self.dirname, self.filename), 'r')
            #self.control.SetValue(f.read())
            wx.MessageBox(f.read())
            f.close()
        dlg.Destroy()

    def onUpdateTime(self,e):
        #wx.MessageBox(self.datetimeValue.GetValue())
        #self.datetime= self.datetimeValue.GetString()
        self.datetime = self.datetimeValue.GetValue()
        self.canvasDraw(None)

    def onUpdateLocation(self,e):
        self.lon = self.lonValue.GetValue()
        self.lat = self.latValue.GetValue()
        self.canvasDraw(None)

        
    def onAddAz(self,e):
        self.az0= (self.az0+15)%360
        self.azValue.SetLabel('{}'.format(self.az0))
        self.canvasDraw(None)
    def onSubAz(self,e):
        self.az0= (self.az0+345)%360
        self.azValue.SetLabel('{}'.format(self.az0))
        self.canvasDraw(None)
    def onAddAlt(self,e):
        if self.alt0 < 80 :
            self.alt0 = self.alt0 + 10
            self.altValue.SetLabel('{}'.format(self.alt0))
            self.canvasDraw(None)

    def onSubAlt(self,e):
        if self.alt0 > -80 :
            self.alt0 = self.alt0 - 10
            self.altValue.SetLabel('{}'.format(self.alt0))
            self.canvasDraw(None)

    def onAddFov(self,e):
        if self.fov < 175 :
            self.fov = self.fov + 5
            self.fovValue.SetLabel('{}'.format(self.fov))
            self.canvasDraw(None)

    def onSubFov(self,e):
        if self.fov > 5 :
            self.fov = self.fov - 5
            self.fovValue.SetLabel('{}'.format(self.fov))
            self.canvasDraw(None)

    def on_pick(self, event):
        line = event.artist
        ind = event.ind
        #self.SetStatusText('{}'.format(self.table[ind[0]]['name']))
        mag=self.table['size'][ind[0]]
        index=0
        for i in range(1,len(ind)):
            if self.table['size'][ind[i]]>mag:
                index=i
                mag=self.table['size'][ind[i]]
        self.starValue.SetLabel('天體：{}, az:{:.4f}, alt={:.4f}'.format(self.table[ind[index]]['name'],self.table[ind[index]]['az'],self.table[ind[index]]['alt']))
    def on_key_press(self, event):
        #print(event.key, event.xdata, event.ydata)
        if event.key=="=":
            self.onSubFov(event)
        if event.key=="-":
            self.onAddFov(event)
        if event.key=="4":
            self.onSubAz(event)
        if event.key=="6":
            self.onAddAz(event)
        if event.key=="2":
            self.onSubAlt(event)
        if event.key=="8":
            self.onAddAlt(event)
        
    def getAzTable(self, az):
        alts=[-80+i*5 for i in range(33)]
        azs= [az for i in range(33)]
        table = Table()
        table['az']=azs
        table['alt']=alts
        return table

    def getAltTable(self, alt):
        alts=[alt for i in range(73)]
        azs= [i*5 for i in range(73)]
        table = Table()
        table['az']=azs
        table['alt']=alts
        return table

    def plot(self,xs,ys,zs,color):
        valid = ys >= np.cos(self.fov/2*np.pi/180)
        preValid = False
        line=[]
        for i in range(len(valid)):
            if valid[i]==True:
                if preValid == False:
                    if i==0:
                        line = [i]
                    else:
                        line=[i-1,i]
                    preValid = True
                else:
                    line.append(i)
            else:
                if preValid==True:
                    line.append(i)
                    self.axes.plot(xs[line], zs[line],color=color)   
                    preValid = False
                    line = []
        if len(line) > 0:
           self.axes.plot(xs[line], zs[line],color=color)             
        #self.axes.plot(xs,zs, color=color)

    def canvasDraw(self, e):
        utctime = localtime2utctime(self.datetime,self.timezone)
        eclipTable = getEclipticTable(utctime)
        location = EarthLocation(lon=float(self.lon)*u.deg, lat= float(self.lat)*u.deg)
        # getPlanetTable is a altaz table with colnames: name, alt, az,Vmag,size
        planetTable = getPlanetTable(utctime, location) 
        self.table = vstack([self.bscTable])
        set_altaz(self.table,utctime, location)
        self.table= self.table['name','alt','az','Vmag','size']
        self.table = vstack([self.table,planetTable])
        self.scale = set_xyz_proj(self.table, self.az0, self.alt0, self.projection, self.fov)
        set_altaz(eclipTable,utctime, EarthLocation(lon=float(self.lon)*u.deg, lat= float(self.lat)*u.deg))
        set_xyz_proj(eclipTable, self.az0, self.alt0, self.projection, self.fov)
        
        self.axes.clear()
        xl = self.scale*16/np.sqrt(337)
        zl = self.scale*9/np.sqrt(337)

        for alt in [-80,-60, -40, -20, 0, 20, 40, 60, 80]:
            table = self.getAltTable(alt)
            set_xyz_proj(table, self.az0, self.alt0, self.projection, self.fov)
            #self.axes.plot(table['x'],table['z'],color='r')
            #mask = table['y'] > np.cos(self.fov/2*np.pi/180)
            #table=table[mask]
            self.plot(table['x'],table['y'],table['z'],'r')

        for az in [i*15 for i in range(24)]:
            table = self.getAzTable(az)
            set_xyz_proj(table, self.az0, self.alt0, self.projection, self.fov)
            #self.axes.plot(table['x'],table['z'],color='r')
            self.plot(table['x'],table['y'],table['z'],'r')

        #self.axes.plot(eclipTable['x'],eclipTable['z'], color='b')
        self.plot(eclipTable['x'],eclipTable['y'],eclipTable['z'],'b')

        #self.axes.scatter(self.table['x'],self.table['z'], s=(self.table['size']*(72./self.figure.dpi))**2, c='w',marker='.', picker=5)
        self.axes.scatter(self.table['x'],self.table['z'], s=self.table['size']*self.table['size'], c='w',marker='.', picker=5)
        self.axes.set_xlim(-xl,xl)
        self.axes.set_ylim(-zl,zl)
        self.axes.set_aspect('equal')
        self.axes.set_facecolor((0, 0, 0))
        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.figure.tight_layout()

        self.canvas.draw()

    def OnPaint(self, event):
        self.axes.clear()
        x=np.array([1,2,3,4])
        self.axes.plot(x, x*x)
        # draw and flush_events
        self.canvas.draw()
        #self.canvas.flush_events()

class App(WIT.InspectableApp):
    def OnInit(self):
        """Create the main window and insert the custom frame."""
        self.Init()
        frame = MyFrame(None, 'Star map')
        frame.Show(True)
        return True

#alternatively you could use
#class App(WIT.InspectableApp):
#    def OnInit(self):
#        """Create the main window and insert the custom frame."""
#        frame = MyFrame(None, 'Star map')
#        frame.Show(True)
#        return True

if __name__ == "__main__":
    app = App()
    app.MainLoop()

