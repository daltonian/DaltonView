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


class FileManager:
    def __init__(self):
        self.openTheoryFileCounts = 0
    
    def openTheoryFile(self,TheoryFilePath):
        self.openTheoryFileCounter()
        try: 
            lines = []
            txtpath = TheoryFilePath.replace('"','')
            with open(txtpath) as f:
                lines = f.readlines()

        except FileNotFoundError:
            print('File not found error... retry selecting the computational output file')

        return lines

    def openTheoryFileCounter(self):
        self.openTheoryFileCounts += 1
        return
    
class DataManager:
    def __init__(self):
        self.counts = 0
    
    def retrieveDataFromLines(self,lines):
        """ Retrieves excitation data from quantum chemical output file 
            Returns list with eV, nm, and Oscillator Strength values """
        #...
        # Find lines which contain eV, oscialltor strength, and FC factor

        if 'Q-Chem' in lines[0]:
            start,end = [x for x,y in zip(range(len(lines)),lines) if 'TDDFT Excitation Energies' in y or 'SETman timing' in y]
            eV = []
            OS = []
            nm = []
            for line in lines[start:end]:
                if 'Excited state' in line:
                    eV.append(float(re.findall('\d*\.?\d+',line)[1]))
                    nm.append(1239.8/eV[-1])
                elif 'Strength' in line:
                    OS.append(float(re.findall('\d*\.?\d+',line)[0]))
                           
        
        elif 'Gaussian' in lines[0]:
            # Find lines which contain eV, oscialltor strength, and FC factor
            start = [x for x in range(len(lines)) if 'Excitation energies and oscillator strengths:' in lines[x]][0]
            end = [x for x in range(len(lines)) if ' Population analysis using the SCF density.' in lines[x]][0]
            states = [x for x in range(start,end) if 'Excited State  ' in lines[x]]        # Find which lines have eV values

            eV = []
            nm = []
            OS = []
            for state in states:
                eV.append(float(re.findall('\d*\.?\d+',lines[state])[1]))
                nm.append(float(re.findall('\d*\.?\d+',lines[state])[2]))
                OS.append(float(re.findall('\d*\.?\d+',lines[state])[3]))
        return eV, nm, OS

    def makeSpectrumFromRetrievedData(self,eV,nm,OS,hwhm,peakShapeBool):
        # Define wavelength (or eV) region over which to create spectrum
        # # Position of wavelengths will effect the normalization of the
        # # spectrum, ie should ideally match experimental ref range
        min_eV = 1500 # Minimum energy in meV (~830 nm)
        max_eV = 12000 # Maximum wavelength in meV (~100 nm)
        num_eV = (max_eV-min_eV)*2+1
        space = np.linspace(min_eV,max_eV,num_eV)/1000

        peak_width = hwhm/1000

        spec_lines = []
        for i in range(len(nm)):
            spec_lines.append(self.peakShape(float(eV[i]),float(OS[i]),peak_width,space,peakShapeBool))
        peaks = np.array(spec_lines).T.tolist()

        # Take each energy's Lorentzian and sum
        spectrum = []
        maxAbs = 0.001
        for element in peaks:
            if sum(element) > maxAbs:
                maxAbs = sum(element)
        maxAbs = 1 # Normalizing turned off
        # Needs two loops to find the total maximum first

        for element in peaks:
            spectrum.append(sum(element)/maxAbs)
        
        return np.array([1239.8/space,spectrum]), np.array([nm,OS])

    def peakShape(self,eV,OS,width,space,peakShapeBool):
        """ Takes eV, osc. strength, hwhm (in nm, converts to eV), 
        and wavelenth range for plotting, and makes indididual peak shape """
        spectrum = []
        if peakShapeBool == 1:
            sigma = width/(2*np.log(2))**0.5
            for x in space:
                spectrum.append(OS*np.exp(-0.5*((x-eV)**2)/(sigma**2))/(sigma*np.sqrt(2*np.pi))) #GAUSSIAN
        elif peakShapeBool == 0:
            for x in space:
                spectrum.append(OS*(width)/(((x-eV)**2)+(width)**2)/np.pi)  #LORENTZIAN
        return spectrum

class PlotManager:
    def __init__(self):
        self.counts = 0

    @staticmethod
    def plotTheoryData(ax,spectrum,sticks,barsBool):
        #... 
        ax.plot(spectrum[0],spectrum[1],linewidth=2)
        if barsBool == 1:
            for x in range(len(sticks[0])):
                print(sticks[0,x],sticks[1,x])
            ax.bar(sticks[0],sticks[1],width=2)
        return
    
    @staticmethod
    def plotExperimentData(path_to_experimental_csv,scale_factor):
        experimental_data = np.loadtxt(path_to_experimental_csv,delimiter=',',skiprows=1)
        ax.plot(experimental_data[:,0],experimental_data[:,1]*float(scale_factor),linewidth=2,color='k')
        return
    
class TextRedirector:
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state=tk.NORMAL)
        self.widget.insert(tk.END, str, (self.tag,))
        self.widget.configure(state=tk.DISABLED)
        self.widget.see(tk.END)  # Autoscroll to the end
    
    def flush(self):
        pass


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
        self.checkvarHoldPlots = tk.IntVar()
        self.checkvarHoldPlots.set(0)
        self.xmax = tk.IntVar()
        self.xmax.set(700)
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

        # ... (continue with the rest of your frame construction)

        topframe.pack(padx=10, pady=10)

        # Middle frame for plot options
        midframe = ttk.LabelFrame(self.root, text='Plot Options')

        checkboxGaussian = ttk.Checkbutton(midframe, text='Gaussian   ', command=self.changed_checkboxGaussian, variable=self.checkvarGaussian)
        checkboxGaussian.grid(column=0, row=0)
        
        checkboxLorentzian = ttk.Checkbutton(midframe, text='Lorentzian', command=self.changed_checkboxLorentzian, variable=self.checkvarLorentzian)
        checkboxLorentzian.grid(column=0, row=1)
        
        checkboxPlotBars = ttk.Checkbutton(midframe, text='Plot bars   ', variable=self.checkvarPlotBars)
        checkboxPlotBars.grid(column=0, row=2)

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
        sys.stdout = TextRedirector(self.console_log, "stdout")


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
        self.ax.set_xlim(100,700)
        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.tick_params(axis='both',labelsize=14)
        self.figure.set(figheight=3.8,figwidth=7,dpi=100)

    def pressedExportData(self):
        self.popupExportReminder()
        return

    def confirmedExportData(self, win):
        #print('Exporting....')
        data = []
        txtpath = self.comp_chem_file_name.get().replace('"','')
        with open(txtpath) as f:
            data = f.readlines()
        spectrum,sticks = self.retrieveDataFromLines(data)
        np.savetxt(txtpath.replace('.out','_spectrum.csv'),spectrum.T,delimiter=',',header="nm,abs")
        np.savetxt(txtpath.replace('.out','_sticks.csv'),sticks.T,delimiter=',',header='nm,f')
        win.destroy()
        return
    
    def popupExportReminder(self):
        win = tk.Toplevel()
        win.wm_title('Export Reminder')
        l0 = ttk.Label(win,text='')
        l0.grid(row=0,column=0,pady=5)
        l1 = ttk.Label(win,text='The exported data will reflect')
        l1.grid(row=1,column=0,padx=10,columnspan=2)
        l2 = ttk.Label(win,text='the current export options!')
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

