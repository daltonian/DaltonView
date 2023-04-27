import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import filedialog as fd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk 
)
import numpy as np
import re

# DaltonView
# Version 1.2
# Release date: 27 April 2023

#@# Create app and set theme #@#
root = tk.Tk()
root.title('DaltonView 1.2')

style = ttk.Style(root)
style.theme_use('xpnative')
root.configure(bg='#f0f0f0')
root.resizable(width=False,height=False)

#@# Defining vars for GUI features #@#
#
root.opens = 0
root.updates = 0
root.spec_change = 0
pathname = tk.StringVar()
data = []
path = []
old_hwhm = [20]


c1var = tk.IntVar()
c1var.set(1)
c2var = tk.IntVar()
c3var = tk.IntVar()
c3var.set(1)
c4var = tk.IntVar()
c4var.set(0)
xmax = tk.IntVar()
xmax.set(700)
xmin = tk.IntVar()
xmin.set(100)
hwhm = tk.IntVar()
hwhm.set(20)
ymax = tk.StringVar()
ymax.set('auto')
#
#@# End definitions #@#

#@# Begin spectral functions #@#
# Functions used for converting txt spec
# data to numerical spectral data
def remove(string):
    return string.replace(" ", "")

def SpecShape(nm,OS,width,space):
    f = []
    for x in space:
        if c1var.get() == 1:
            width = width/(2*np.log(2))**0.5
            f.append(OS*np.exp(-0.5*((x-nm)**2)/(width**2))/(width*np.sqrt(2*np.pi))) #GAUSSIAN
        elif c1var.get() == 0:
            f.append(OS*(width)/(((x-nm)**2)+(width)**2)/np.pi)  #LORENTZIAN
    return f

def spec_gen(outfile):
    #...
    # Find lines which contain eV, oscialltor strength, and FC factor

    if 'Q-Chem' in outfile[0]:
        start = [x for x in range(len(outfile)) if 'TDDFT Excitation Energies' in outfile[x]][0]
        end = [x for x in range(len(outfile)) if 'SETman timing' in outfile[x]][0]
        indEV = [x for x in range(start,end) if 'Excited state' in outfile[x]]        # Find which lines have eV values
        indOS = [x for x in range(start,end) if 'Strength' in outfile[x]]    # Find which lines have oscillator Strength
        # Create list of eV values from TDDFT/TDA result
        eV = []
        for i in range(len(indEV)):
            eV.append(float(re.findall('\d*\.?\d+',outfile[int(indEV[i])])[1]))
        # Create list of nm values
        nm = []
        for i in eV:
            if i < 50:
                nm.append((1.2398e3)/i)            # Convert eV to nm
        # Generate list of oscillator strength values
        OS = []
        for i in range(len(nm)):
            OS.append(float(re.findall('\d*\.?\d+',outfile[int(indOS[i])])[0]))
    
    elif 'Gaussian' in outfile[0]:
        # Find lines which contain eV, oscialltor strength, and FC factor
        start = [x for x in range(len(outfile)) if 'Excitation energies and oscillator strengths:' in outfile[x]][0]
        end = [x for x in range(len(outfile)) if ' Population analysis using the SCF density.' in outfile[x]][0]
        states = [x for x in range(start,end) if 'Excited State  ' in outfile[x]]        # Find which lines have eV values

        eV = []
        nm = []
        OS = []
        for state in states:
            eV.append(float(re.findall('\d*\.?\d+',outfile[state])[1]))
            nm.append(float(re.findall('\d*\.?\d+',outfile[state])[2]))
            OS.append(float(re.findall('\d*\.?\d+',outfile[state])[3]))


    # Define wavelength (or eV) region over which to create spectrum
    # # Position of wavelengths will effect the normalization of the
    # # spectrum, ie should ideally match experimental ref range
    min_wl = 1 # Minimum wavelength in nm
    max_wl = 700 # Maximum wavelength in nm
    num_wl = (max_wl-min_wl)*2+1
    space = np.linspace(min_wl,max_wl,num_wl)

    gamma = int(hwhm.get())

    spec_lines = []
    for i in range(len(nm)):
        spec_lines.append(SpecShape(float(nm[i]),float(OS[i]),gamma,space))
    peaks = np.array(spec_lines).T.tolist()

    # Take each energy's Lorentzian and sum
    spectrum = []
    maxAbs = 0.001
    for element in peaks:
        if sum(element) > maxAbs:
            maxAbs = sum(element)
    maxAbs = 1
    # Needs two loops to find the total maximum first
    for element in peaks:
        spectrum.append(sum(element)/maxAbs)
    return np.array([space,spectrum]), np.array([nm,OS])

#@# End spectral functions #@#

#@# Begin functions #@#

def select_file():
    filetypes = (
        ('Output files', ['*.out','*.log']),
        ('All files', '*.*')
    )

    f = fd.askopenfile(
        title='Open a file',
        initialdir='/',
        filetypes=filetypes)
    
    try:
        data=f.readlines()
        pathname.set(f.name)
    except AttributeError:
        'NULL'
    
    return

def opened():
    root.opens += 1
    return

def do_nothing():
    return 'break'

def plottr(spectrum,sticks):
    #... 
    line = plt.plot(spectrum[0],spectrum[1]*70,linewidth=2)
    if c3var.get() == 1: plt.bar(sticks[0],sticks[1],width=2,color=ax.lines[-1]._color)
    ax.set_xlabel('wavelength / nm')
    ax.set_ylabel('absorbance (arb. units)')
    ax.set_xlim(xmin=xmin.get(),xmax=xmax.get())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    return

def c1_changed():
    c2var.set(0)
    c1var.set(1)
    root.spec_change = True
    return

def c2_changed():
    c1var.set(0)
    c2var.set(1)
    root.spec_change = True
    return

def open_file():
    try:
        opened()
        data = []
        txtpath = path.get().replace('"','')
        with open(txtpath) as f:
            data = f.readlines()
            f.close()
        
        spectrum, sticks = np.zeros((1,1)), np.zeros((1,1))
        spectrum, sticks  = spec_gen(data)
        
        plottr(spectrum, sticks)
        canvas.draw()
    except FileNotFoundError:
        pass
    return

def slider_changed(event):
    hwhm.set(round(slider.get()))
    value_label.configure(text=get_current_hwhm())
    return

def get_current_hwhm():
    return str(hwhm.get())+' nm'
    
def update_plot():
    root.updates += 1
    old_hwhm.append(int(hwhm.get()))
    ax.set_xlim(xmin=xmin.get(),xmax=xmax.get())
    try:
        ax.set_ylim(ymin=0,ymax=float(ymax.get()))
    except ValueError:
        pass
    if ((hwhm.get()-old_hwhm[root.updates-1]) != 0) or (root.spec_change==True):
        if c4var.get() != 1: ax.clear()
        open_file()
    spec_change = False
    return canvas.draw()

def confirm_export():
    popupBonus()
    return

def export(win):
    #print('Exporting....')
    data = []
    txtpath = pathname.get().replace('"','')
    with open(txtpath) as f:
        data = f.readlines()
        f.close()
    spectrum,sticks = spec_gen(data)
    np.savetxt(txtpath.replace('.out','_spectrum.csv'),spectrum.T,delimiter=',',header="nm,abs")
    np.savetxt(txtpath.replace('.out','_sticks.csv'),sticks.T,delimiter=',',header='nm,f')
    win.destroy()
    return


#@# Begin setting window size #@#
#

# get the screen dimension
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
# find the center point
center_x = int(screen_width/2 - 700 / 2)
center_y = int(screen_height/2 - 650 / 2)
# set window dimensions
root.geometry(f'{700}x{650}+{center_x}+{center_y}')

#
#@# End window size #@#


## Assembling the GUI ##
#

# Top frame for file input and processing
topframe = ttk.Frame(root)
topframe.grid(column=0,row=0)
ttk.Label(topframe, text='Output file path:').grid(column=0,row=0)
path = tk.Entry(topframe,textvariable=pathname,width=50)
path.grid(column=1,row=0)
f = ttk.Button(topframe,text='Select Output File',command=select_file).grid(column=2,row=0)
file_button = ttk.Button(topframe, text='Plot File',command=open_file).grid(column=3,row=0)
topframe.pack(padx=10,pady=10)

# middle frame for plot options
midframe = ttk.LabelFrame(root,text='Plot Options')
c1 = ttk.Checkbutton(midframe,
                text='Gaussian         ',
                command=c1_changed,
                variable=c1var
                )
c1.grid(column=0,row=0)
c2 = ttk.Checkbutton(midframe,
                text='Lorentzian      ',
                command=c2_changed,
                variable=c2var
                )
c2.grid(column=0,row=1)
c3 = ttk.Checkbutton(midframe,
                text='Plot bars         ',
                variable=c3var
                )
c3.grid(column=0,row=2)

ttk.Label(midframe, text='xmin: ').grid(column=1,row=0)
xminbox = tk.Entry(midframe,textvariable=xmin,width=8)
xminbox.bind('<Return>',do_nothing())
xminbox.grid(column=2,row=0)
ttk.Label(midframe, text='xmax: ').grid(column=1,row=1)
xmaxbox = tk.Entry(midframe,textvariable=xmax,width=8)
xmaxbox.bind('<Return>',do_nothing())
xmaxbox.grid(column=2,row=1)
ttk.Label(midframe, text='ymax: ').grid(column=1,row=3)
ymaxbox = tk.Entry(midframe,textvariable=ymax,width=8)
ymaxbox.bind('<Return>',do_nothing())
ymaxbox.grid(column=2,row=3)

ttk.Label(midframe, text='     ').grid(column=3,row=0)
ttk.Label(midframe, text='     ').grid(column=3,row=1)
ttk.Label(midframe, text='1 nm').grid(column=4,row=1)
ttk.Label(midframe, text='Peak Width (HWHM)').grid(column=5,row=0)

# Slider for adjusting peak width
slider = ttk.Scale(
    midframe,
    from_=1,
    to=80,
    orient='horizontal',
    command=slider_changed,
    variable=hwhm
)
slider.grid(row=1,column=5)

value_label = ttk.Label(
    midframe,
    text=get_current_hwhm()
)
value_label.grid(row=2,column=5)

ttk.Label(midframe, text='80 nm    ').grid(column=6,row=1)

# Button for updating the figure
update_button = ttk.Button(
    midframe,
    command=update_plot,
    text="Update Plot"
)
update_button.grid(row=1,column=7)
c4 = ttk.Checkbutton(
    midframe,
    variable=c4var,
    text='Hold all plots'
)
c4.grid(row=3,column=7)
midframe.pack(padx=10,pady=10,ipadx=10,ipady=10)


# plot section - free floating
figure,ax = plt.subplots()
global line
figure.subplots_adjust(top=0.92,right=0.95,bottom=0.16,left=0.15)
ax.set_xlabel('wavelength / nm',size=14)
ax.set_ylabel('absorbance (arb. units)',size=14)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
figure.set(figheight=3.8,figwidth=7,dpi=100)
canvas = FigureCanvasTkAgg(figure, root)
canvas.get_tk_widget().pack(pady=20)


def canceled(win):
    export_exception = True
    return win.destroy()

def popupBonus():
    win = tk.Toplevel()
    win.wm_title('Export Reminder')
    l0 = ttk.Label(win,text='')
    l0.grid(row=0,column=0,pady=5)
    l1 = ttk.Label(win,text='The exported data will reflect')
    l1.grid(row=1,column=0,padx=10,columnspan=2)
    l2 = ttk.Label(win,text='the current export options!')
    l2.grid(row=2,column=0,columnspan=2)
    b1 = ttk.Button(win,text='Okay',command=lambda win= win:export(win))
    b1.grid(row=3,column=0,pady=10)
    b2 = ttk.Button(win,text='Cancel',command=lambda win = win:canceled(win))
    b2.grid(row=3,column=1)
    root_x = root.winfo_rootx()
    root_y = root.winfo_rooty()
    win.geometry(f'+{root_x+250}+{root_y+220}')
    return win

def save_fig():
    txtpath = pathname.get().replace('"','')
    plt.savefig(txtpath.replace('.out','.png'),dpi=500)
    return


# Frame for exporting functions
bottom_frame = ttk.Frame(root)

# Save figure button
save_fig_button = ttk.Button(
    bottom_frame,
    text='Save Figure Image',
    command=save_fig
)
save_fig_button.grid(column=0,row=0,padx=5)

# Export button
global export_exception
export_exception = False
export_button = ttk.Button(
    bottom_frame,
    command=confirm_export,
    text='Export Spectrum'
)
export_button.grid(column=1,row=0)

bottom_frame.pack()


#@# End GUI #@#

#@# Fin #@#

# Allow for people to change their minds
def on_closing():
    if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
        plt.close()
        root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()







