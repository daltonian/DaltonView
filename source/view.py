"""View: Tkinter GUI - handles only display and input capture"""
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import AutoMinorLocator

class SpectrumView(tk.Tk):
    """Main application view"""
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        
        # store references for updating
        self._canvas = None
        self._console_log = None
        self._widgets = {}  # store widget references
        
        self.title('DaltonView 2.0')
        self.configure(bg='#f0f0f0')
        self.resizable(width=False, height=False)
        
        # center window
        self._center_window(700, 750)
        
        self._setup_gui()
        self.protocol("WM_DELETE_WINDOW", self._on_closing)
    
    def _center_window(self, width, height):
        """Center window on screen"""
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        center_x = int(screen_width / 2 - width / 2)
        center_y = int(screen_height / 2 - height / 2)
        self.geometry(f'{width}x{height}+{center_x}+{center_y}')
    
    def _setup_gui(self):
        """Build UI components"""
        # view is composed of 3 Tkinter frames, hereby called Top, Middle, and Bottom.

        # top frame - file selection
        self._setup_file_frame()
        
        # middle frame - plot options
        self._setup_options_frame()
        
        # canvas for matplotlib
        self._setup_canvas_frame()
        
        # bottom frame - export buttons
        self._setup_export_frame()
        
        # console log
        self._setup_console_frame()
    
    def _setup_file_frame(self):
        """File selection UI"""
        topframe = ttk.Frame(self)
        topframe.grid(column=0, row=0)
        
        # add text box for quantum file; allowing entry if user wants to paste it in manually
        ttk.Label(topframe, text='Output file:').grid(column=0, row=0)
        self._widgets['theory_path_var'] = tk.StringVar()
        theory_entry = tk.Entry(topframe, textvariable=self._widgets['theory_path_var'], width=50)
        theory_entry.grid(column=1, row=0)
        self._widgets['theory_entry'] = theory_entry
        
        # same for expermental file
        ttk.Label(topframe, text='Experimental Reference file: ').grid(column=0, row=1)
        self._widgets['exp_path_var'] = tk.StringVar()
        exp_entry = tk.Entry(topframe, textvariable=self._widgets['exp_path_var'], width=50)
        exp_entry.grid(column=1, row=1)
        self._widgets['exp_entry'] = exp_entry
        
        ttk.Button(topframe, text='Select Output File', 
                  command=self._on_select_theory_file, width=20).grid(column=2, row=0)
        ttk.Button(topframe, text='Select Experiment File', 
                  command=self._on_select_experiment_file, width=20).grid(column=2, row=1)
        ttk.Button(topframe, text='Make Plot', 
                  command=self._on_make_plot).grid(column=3, row=0)
        
        topframe.pack(padx=10, pady=10)
    
    def _setup_options_frame(self):
        """Plot options UI"""
        midframe = ttk.LabelFrame(self, text='Plot Options')
        
        # checkboxes for various options
        self._widgets['gaussian_var'] = tk.IntVar(value=1)
        self._widgets['lorentzian_var'] = tk.IntVar(value=0)
        self._widgets['plot_bars_var'] = tk.IntVar(value=1)
        self._widgets['hold_plots_var'] = tk.IntVar(value=0)
        
        ttk.Checkbutton(midframe, text='Gaussian', 
                       variable=self._widgets['gaussian_var'],
                       command=self._on_gaussian_changed).grid(column=0, row=0)
        ttk.Checkbutton(midframe, text='Lorentzian',
                       variable=self._widgets['lorentzian_var'],
                       command=self._on_lorentzian_changed).grid(column=0, row=1)
        ttk.Checkbutton(midframe, text='Plot bars',
                       variable=self._widgets['plot_bars_var']).grid(column=0, row=2)
        
        # range inputs for user specificed min and max fro axes
        ttk.Label(midframe, text='xmin:').grid(column=1, row=0)
        self._widgets['xmin_var'] = tk.IntVar(value=100)
        ttk.Entry(midframe, textvariable=self._widgets['xmin_var'], width=8).grid(column=2, row=0)
        
        ttk.Label(midframe, text='xmax:').grid(column=1, row=1)
        self._widgets['xmax_var'] = tk.IntVar(value=2000)
        ttk.Entry(midframe, textvariable=self._widgets['xmax_var'], width=8).grid(column=2, row=1)
        
        ttk.Label(midframe, text='ymax:').grid(column=1, row=2)
        self._widgets['ymax_var'] = tk.StringVar(value='auto')
        ttk.Entry(midframe, textvariable=self._widgets['ymax_var'], width=8).grid(column=2, row=2)
        
        ttk.Label(midframe, text='scale_exp:').grid(column=1, row=3)
        self._widgets['scale_exp_var'] = tk.StringVar(value='1')
        ttk.Entry(midframe, textvariable=self._widgets['scale_exp_var'], width=8).grid(column=2, row=3)
        
        # slider for HWHM (half width at half maximum)
        ttk.Label(midframe, text='     ').grid(column=3, row=0)
        ttk.Label(midframe, text='0.100 eV').grid(column=4, row=1)
        ttk.Label(midframe, text='Peak Width (HWHM)').grid(column=5, row=0)
        
        self._widgets['hwhm_var'] = tk.IntVar(value=200)
        slider = ttk.Scale(midframe, from_=100, to=400, orient='horizontal',
                          variable=self._widgets['hwhm_var'],
                          command=self._on_hwhm_changed)
        slider.grid(row=1, column=5)
        
        self._widgets['hwhm_label'] = ttk.Label(midframe, text=self._get_hwhm_text())
        self._widgets['hwhm_label'].grid(row=2, column=5)
        
        ttk.Label(midframe, text='0.400 eV').grid(column=6, row=1)
        
        ttk.Button(midframe, text='Update Plot',
                  command=self._on_update_plot).grid(row=1, column=7)
        ttk.Checkbutton(midframe, text='Hold all plots',
                       variable=self._widgets['hold_plots_var']).grid(row=3, column=7)
        
        midframe.pack(padx=10, pady=10, ipadx=10, ipady=10)
    
    def _setup_canvas_frame(self):
        """Matplotlib canvas"""
        # create figure
        figure, ax = plt.subplots()
        self.figure = figure
        self.ax = ax
        self._initialize_figure()
        
        # create canvas
        self._canvas = FigureCanvasTkAgg(figure, self)
        self._canvas.get_tk_widget().pack(pady=20)
    
    def _setup_export_frame(self):
        """Export buttons"""
        bottom_frame = ttk.Frame(self)
        ttk.Button(bottom_frame, text='Save Figure Image',
                  command=self._on_save_figure).grid(column=0, row=0, padx=5)
        ttk.Button(bottom_frame, text='Export Spectrum',
                  command=self._on_export_data).grid(column=1, row=0)
        bottom_frame.pack()
    
    def _setup_console_frame(self):
        """Console log output"""
        console_frame = ttk.LabelFrame(self, text='Console Log')
        self._console_log = scrolledtext.ScrolledText(console_frame, wrap=tk.WORD, 
                                                      height=2, state=tk.DISABLED)
        self._console_log.grid(column=0, row=0, pady=0)
        console_frame.pack()
        
        # make program sdout go into the gui
        import sys
        sys.stdout = ConsoleRedirector(self._console_log)
    
    def _initialize_figure(self):
        """Set up matplotlib figure"""
        self.figure.subplots_adjust(top=0.92, right=0.95, bottom=0.16, left=0.15)
        self.ax.set_xlabel('wavelength / nm', size=14)
        self.ax.set_ylabel('normalized absorbance', size=14)
        self.ax.set_xlim(100, 2000)
        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.tick_params(axis='both', labelsize=14)
        self.figure.set(figheight=3.8, figwidth=7, dpi=100)
    
    # event handlers - all delegate to controller

    # runs when user closes the open file dialog
    def _on_select_theory_file(self):
        filetypes = (('Output files', ['*.out', '*.log']), ('All files', '*.*'))
        file = filedialog.askopenfile(title='Computational Output File', 
                                      initialdir='/', filetypes=filetypes)
        if file:
            self._widgets['theory_path_var'].set(file.name)
            self.controller.on_theory_file_selected(file.name)
            file.close()
    
    # same but when user closes open experimental file dialog
    def _on_select_experiment_file(self):
        filetypes = [('Comma Separated Variables', '*.csv')]
        file = filedialog.askopenfile(title='Spectrum CSV File [wavelength,abs]',
                                      initialdir='/', filetypes=filetypes)
        if file:
            self._widgets['exp_path_var'].set(file.name)
            self.controller.on_experiment_file_selected(file.name)
            file.close()
    
    # user presses to load plot manually
    def _on_make_plot(self):
        self.controller.on_make_plot_clicked()
    
    # auto update plot if peak shape type is changed
    # this could potentially be combined with method below?
    def _on_gaussian_changed(self):
        if self._widgets['gaussian_var'].get():
            self._widgets['lorentzian_var'].set(0)
        self.controller.on_peak_shape_changed('gaussian', True)
    
    def _on_lorentzian_changed(self):
        if self._widgets['lorentzian_var'].get():
            self._widgets['gaussian_var'].set(0)
        self.controller.on_peak_shape_changed('lorentzian', True)
    
    # update plot on hwhm changed
    def _on_hwhm_changed(self, value):
        self._widgets['hwhm_label'].config(text=self._get_hwhm_text())
    
    # user presses to change plot manualy
    def _on_update_plot(self):
        self.controller.on_update_plot_clicked()
    
    # saves figure
    def _on_save_figure(self):
        self.controller.on_save_figure_clicked()
    
    # export to csv (our bread-and-butter :D)
    def _on_export_data(self):
        if messagebox.askyesno("Export", "Export spectrum data?"):
            self.controller.on_export_data_clicked()
    
    # closes DaltonView
    def _on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.destroy()
    
    def _get_hwhm_text(self):
        return f"{self._widgets['hwhm_var'].get() / 1000:.3f} eV"
    
    # Update methods called by controller
    def update_canvas(self):
        """Redraw matplotlib canvas"""
        self._canvas.draw()
    
    def clear_plot(self):
        """Clear the plot"""
        self.ax.clear()
        self._initialize_figure()
    
    def display_message(self, message):
        """Show message in console"""
        print(message)
    
    def show_error(self, error_msg):
        """Show error dialog"""
        messagebox.showerror("Error", error_msg)
    
    def get_all_settings(self):
        """Return all current UI settings"""
        return {
            'xmin': self._widgets['xmin_var'].get(),
            'xmax': self._widgets['xmax_var'].get(),
            'ymax': self._widgets['ymax_var'].get(),
            'hwhm': self._widgets['hwhm_var'].get(),
            'gaussian_enabled': self._widgets['gaussian_var'].get() == 1,
            'lorentzian_enabled': self._widgets['lorentzian_var'].get() == 1,
            'plot_bars': self._widgets['plot_bars_var'].get() == 1,
            'hold_plots': self._widgets['hold_plots_var'].get() == 1,
            'scale_exp': float(self._widgets['scale_exp_var'].get()),
        }


class ConsoleRedirector:
    """Redirect print output to Tkinter widget"""
    def __init__(self, widget):
        self.widget = widget
    
    def write(self, message):
        self.widget.config(state=tk.NORMAL)
        self.widget.insert(tk.END, message)
        self.widget.config(state=tk.DISABLED)
        self.widget.see(tk.END)
    
    def flush(self):
        pass