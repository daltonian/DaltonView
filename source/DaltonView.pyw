import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import re
import sys

from FileManager import FileManager
from DataManager import DataManager
from LogManager import LogManager
from PlotManager import PlotManager

class DaltonViewApp:
    def __init__(self,root,figure,ax):
        self.root = root
        self.fileManager = FileManager()
        self.dataManager = DataManager()
        
        self.figure = figure
        self.ax = ax
        self.initializeVars()
        self.setupGUI()
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def initializeVars(self):
        self.openTheoryFileCounts = 0
        self.root.updates = 0
        self.root.spec_change = 0
        self.comp_chem_file_name = tk.StringVar()
        self.experiment_file_name = tk.StringVar()
        self.data = []
        self.experimental_data = []
        self.old_hwhm = [200]
        self.checkvarGaussian = tk.IntVar()
        self.checkvarGaussian.set(1)
        self.checkvarLorentzian = tk.IntVar()
        self.checkvarPlotBars = tk.IntVar()
        self.checkvarPlotBars.set(1)
        self.checkvarFullTDDFT = tk.IntVar()
        self.checkvarFullTDDFT.set(0)
        self.checkvarHoldPlots = tk.IntVar()
        self.checkvarHoldPlots.set(0)
        self.xmax = tk.IntVar()
        self.xmax.set(2000)
        self.xmin = tk.IntVar()
        self.xmin.set(100)
        self.hwhm = tk.IntVar()
        self.hwhm.set(200)
        self.ymax = tk.StringVar()
        self.ymax.set('auto')
        self.scale_exp = tk.StringVar()
        self.scale_exp.set('1')
    
    def setupGUI(self):
        # Set up the main window
        self.root.title('DaltonView 2.0')
        self.root.configure(bg='#f0f0f0')
        self.root.resizable(width=False, height=False)
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        # find the center point
        center_x = int(screen_width/2 - 700 / 2)
        center_y = int(screen_height/2 - 750 / 2)
        # set window dimensions
        root.geometry(f'{700}x{750}+{center_x}+{center_y}')

        # Top frame for file input and processing
        topframe = ttk.Frame(self.root)
        topframe.grid(column=0, row=0)
        
        ttk.Label(topframe, text='Output file:').grid(column=0, row=0)
        self.TheoryFilePath = tk.Entry(topframe, textvariable=self.comp_chem_file_name, width=50)
        self.TheoryFilePath.grid(column=1, row=0)
        
        ttk.Label(topframe, text='Experimental Reference file: ').grid(column=0, row=1)
        self.ExperimentFilePath = tk.Entry(topframe, textvariable=self.experiment_file_name, width=50)
        self.ExperimentFilePath.grid(column=1, row=1)
        
        buttonSelectTheoryFile = ttk.Button(topframe, text='Select Output File', command=self.pressedSelectTheoryFile, width=20)
        buttonSelectTheoryFile.grid(column=2, row=0)
        
        buttonSelectExperimentFile = ttk.Button(topframe, text='Select Experiment File', command=self.pressedSelectExperimentFile, width=20)
        buttonSelectExperimentFile.grid(column=2, row=1)
        
        buttonMakePlot = ttk.Button(topframe, text='Make Plot', command=self.pressedMakePlot)
        buttonMakePlot.grid(column=3, row=0)

        topframe.pack(padx=10, pady=10)

        midframe = ttk.LabelFrame(self.root, text='Plot Options')

        checkboxGaussian = ttk.Checkbutton(midframe, text='Gaussian    ', command=self.changed_checkboxGaussian, variable=self.checkvarGaussian)
        checkboxGaussian.grid(column=0, row=0)
        
        checkboxLorentzian = ttk.Checkbutton(midframe, text='Lorentzian ', command=self.changed_checkboxLorentzian, variable=self.checkvarLorentzian)
        checkboxLorentzian.grid(column=0, row=1)
        
        checkboxPlotBars = ttk.Checkbutton(midframe, text='Plot bars    ', variable=self.checkvarPlotBars)
        checkboxPlotBars.grid(column=0, row=2)

        # checkboxFullTDDFT = ttk.Checkbutton(midframe, text='Use full TD ', variable=self.checkvarFullTDDFT)
        # checkboxFullTDDFT.grid(column=0, row=3)

        ttk.Label(midframe, text='xmin: ').grid(column=1,row=0)
        xminbox = tk.Entry(midframe,textvariable=self.xmin,width=8)
        #xminbox.bind('<Return>',self.do_nothing())
        xminbox.grid(column=2,row=0)
        ttk.Label(midframe, text='xmax: ').grid(column=1,row=1)
        xmaxbox = tk.Entry(midframe,textvariable=self.xmax,width=8)
        #xmaxbox.bind('<Return>',self.do_nothing())
        xmaxbox.grid(column=2,row=1)
        ttk.Label(midframe, text='ymax: ').grid(column=1,row=2)
        ymaxbox = tk.Entry(midframe,textvariable=self.ymax,width=8)
        #ymaxbox.bind('<Return>',self.do_nothing())
        ymaxbox.grid(column=2,row=2)
        ttk.Label(midframe, text='scale_exp:').grid(column=1,row=3)
        ymaxbox = tk.Entry(midframe,textvariable=self.scale_exp,width=8)
        #ymaxbox.bind('<Return>',self.do_nothing())
        ymaxbox.grid(column=2,row=3)

        ttk.Label(midframe, text='     ').grid(column=3,row=0)
        ttk.Label(midframe, text='     ').grid(column=3,row=1)
        ttk.Label(midframe, text='0.100 eV ').grid(column=4,row=1)
        ttk.Label(midframe, text='Peak Width (HWHM)').grid(column=5,row=0)

        # Slider for adjusting peak width
        self.slider = ttk.Scale(
            midframe,
            from_=100,
            to=400,
            orient='horizontal',
            command=self.changed_slider,
            variable=self.hwhm
        )
        self.slider.grid(row=1,column=5)

        self.value_label = ttk.Label(
            midframe,
            text=self.get_current_hwhm()
        )
        self.value_label.grid(row=2,column=5)

        ttk.Label(midframe, text='0.400 eV    ').grid(column=6,row=1)

        # Button for updating the figure
        buttonUpdatePlot = ttk.Button(
            midframe,
            command=self.pressedUpdatePlot,
            text="Update Plot"
        )
        buttonUpdatePlot.grid(row=1,column=7)

        # Option for keeping old traces on figure
        checkboxHoldPlots = ttk.Checkbutton(
            midframe,
            variable=self.checkvarHoldPlots,
            text='Hold all plots'
        )
        checkboxHoldPlots.grid(row=3,column=7)

        midframe.pack(padx=10, pady=10, ipadx=10, ipady=10)

        self.initiateFigures()

        self.canvas = FigureCanvasTkAgg(figure, self.root)
        self.canvas.get_tk_widget().pack(pady=20)

        # Frame for exporting functions
        bottom_frame = ttk.Frame(self.root)

        buttonSaveFigure = ttk.Button(bottom_frame, text='Save Figure Image', command=self.pressedSaveFigure)
        buttonSaveFigure.grid(column=0, row=0, padx=5)

        buttonExportData = ttk.Button(bottom_frame, command=self.pressedExportData, text='Export Spectrum')
        buttonExportData.grid(column=1, row=0)

        bottom_frame.pack()

        # Create Text widget for console log
        console_frame = ttk.LabelFrame(self.root,text='Console Log')
        self.console_log = scrolledtext.ScrolledText(console_frame, wrap=tk.WORD, height=2, state=tk.DISABLED)
        self.console_log.grid(column=0,row=0,pady=0)
        console_frame.pack()

        # Redirect stdout to the console log
        sys.stdout = LogManager(self.console_log, "stdout")


    def pressedSelectTheoryFile(self):
        filetypes = (
            ('Output files', ['*.out','*.log']),
            ('All files', '*.*')
        )

        theoryFile = filedialog.askopenfile(
            title='Computational Output File',
            initialdir='/',
            filetypes=filetypes)
        
        try:
            theoryFile.readlines()
            self.comp_chem_file_name.set(theoryFile.name)
            theoryFile.close()
        except AttributeError:
            'NULL'
        return

    def pressedSelectExperimentFile(self):
        filetypes = [('Comma Separated Variables','*.csv')]
        experimentFile = filedialog.askopenfile(
            title='Spectrum CSV File [wavelength,abs]',
            initialdir='/',
            filetypes=filetypes)
        try:
            experimentFile.readlines()
            self.experiment_file_name.set(experimentFile.name)
            experimentFile.close()
        except AttributeError:
            'NULL'
    
        return

    def pressedMakePlot(self):
        linesFromTheoryFile = self.fileManager.openTheoryFile(TheoryFilePath=self.TheoryFilePath.get())
        self.eV, self.nm , self.OS  = self.dataManager.retrieveDataFromLines(linesFromTheoryFile)
        self.spectrum, self.sticks = self.dataManager.makeSpectrumFromRetrievedData(eV=self.eV,nm=self.nm,OS=self.OS,hwhm=self.hwhm.get(),peakShapeBool=self.checkvarGaussian.get())
        PlotManager.plotTheoryData(self.ax,self.spectrum,self.sticks,self.checkvarPlotBars.get())
        self.canvas.draw()

        if self.ExperimentFilePath.get():
            self.plotExperimentData(self.ExperimentFilePath.get().replace('"',''),self.scale_exp.get())
        return

    def pressedUpdatePlot(self):
        self.root.updates += 1
        self.old_hwhm.append(float(self.hwhm.get()))
        if ((float(self.hwhm.get())-self.old_hwhm[self.root.updates-1]) != 0) or (self.root.spec_change==True):
            if self.checkvarHoldPlots.get() == 0: 
                print('clearing')
                self.ax.clear()
            self.spectrum, self.sticks = self.dataManager.makeSpectrumFromRetrievedData(eV=self.eV,nm=self.nm,OS=self.OS,hwhm=self.hwhm.get(),peakShapeBool=self.checkvarGaussian.get())
            PlotManager.plotTheoryData(self.ax,self.spectrum,self.sticks,self.checkvarPlotBars.get())
        self.spec_change = False
        if self.ExperimentFilePath.get():
            PlotManager.plotExperimentData(self.ExperimentFilePath.get().replace('"',''),self.scale_exp.get())
        if self.ymax.get() != 'auto':
            self.ax.set_ylim(ymin=0,ymax=float(self.ymax.get()))

        self.ax.set_xlim(xmin=self.xmin.get(),xmax=self.xmax.get())
        self.figure.subplots_adjust(top=0.92,right=0.95,bottom=0.16,left=0.15)
        self.ax.set_xlabel('wavelength / nm',size=14)
        self.ax.set_ylabel('normalized absorbance',size=14)
        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.legend(edgecolor='none')
        self.ax.tick_params(axis='both',labelsize=14)
        self.figure.set(figheight=3.8,figwidth=7,dpi=100)
        
        self.canvas.draw()
        return 

    def changed_checkboxGaussian(self):
        self.checkvarLorentzian.set(0)
        self.checkvarGaussian.set(1)
        self.root.spec_change = True
        return

    def changed_checkboxLorentzian(self):
        self.checkvarGaussian.set(0)
        self.checkvarLorentzian.set(1)
        self.root.spec_change = True
        return

    def changed_slider(self, event):
        self.hwhm.set(round(self.slider.get(),3))
        self.value_label.configure(text=self.get_current_hwhm())
        return

    def get_current_hwhm(self):
        return str(self.hwhm.get()/1000)+' eV'

    def initiateFigures(self):
        self.figure.subplots_adjust(top=0.92,right=0.95,bottom=0.16,left=0.15)
        self.ax.set_xlabel('wavelength / nm',size=14)
        self.ax.set_ylabel('normalized absorbance',size=14)
        self.ax.set_xlim(100,2000)
        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.tick_params(axis='both',labelsize=14)
        self.figure.set(figheight=3.8,figwidth=7,dpi=100)
        return

    def pressedExportData(self):
        self.popupExportReminder()
        return

    def confirmedExportData(self, win):
        txtpath = self.comp_chem_file_name.get().replace('"','')
        basefilename = txtpath.split(".")
        spectrum,sticks = self.spectrum,self.sticks
        np.savetxt(basefilename[0]+'_spectrum.csv',spectrum.T,delimiter=',',header="nm,abs")
        np.savetxt(basefilename[0]+'_sticks.csv',sticks.T,delimiter=',',header='nm,f')
        win.destroy()
        return
    
    def popupExportReminder(self):
        win = tk.Toplevel()
        win.wm_title('Export Reminder')
        l0 = ttk.Label(win,text='')
        l0.grid(row=0,column=0,pady=5)
        l1 = ttk.Label(win,text='The exported data will reflect')
        l1.grid(row=1,column=0,padx=10,columnspan=2)
        l2 = ttk.Label(win,text='the last \'updated\' data!')
        l2.grid(row=2,column=0,columnspan=2)
        b1 = ttk.Button(win,text='Okay',command=lambda win= win:self.confirmedExportData(win))
        b1.grid(row=3,column=0,pady=10)
        b2 = ttk.Button(win,text='Cancel',command=lambda win = win:self.canceled(win))
        b2.grid(row=3,column=1)
        root_x = self.root.winfo_rootx()
        root_y = self.root.winfo_rooty()
        win.geometry(f'+{root_x+250}+{root_y+220}')
        return win

    def pressedSaveFigure(self):
        txtpath = self.comp_chem_file_name.get().replace('"','')
        plt.savefig(txtpath.replace('.out','.png'),dpi=500)
        return

    def on_closing(self):
        if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            plt.close()
            self.root.destroy()
        return


if __name__ == "__main__":
    root = tk.Tk()
    figure, ax = plt.subplots()
    app = DaltonViewApp(root,figure,ax)
    root.mainloop()

