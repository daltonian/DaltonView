"""Controller: Orchestrates Model and View interaction"""

class SpectrumController:
    """Handles communication between Model and View"""
    
    def __init__(self, model, view):
        self.model = model
        self.view = view
        
        # Register as observer
        self.model.attach_observer(self)
        
        # Reference to matplotlib axis
        self.ax = view.ax
    
    # File operations
    def on_theory_file_selected(self, path):
        """Handle theory file selection"""
        if self.model.set_theory_file(path):
            self.view.display_message(f"Loaded theory file: {path}")
        else:
            self.view.show_error("Failed to load theory file")
    
    def on_experiment_file_selected(self, path):
        """Handle experiment file selection"""
        self.model.set_experiment_file(path)
        self.view.display_message(f"Set experiment file: {path}")
    
    # Plot operations
    def on_make_plot_clicked(self):
        """Generate initial plot"""
        settings = self.view.get_all_settings()
        self._apply_settings(settings)
        
        if self.model.regenerate_spectrum():
            self._plot_spectrum()
            self._load_experiment_if_available(settings)
        else:
            self.view.show_error("Failed to generate spectrum")
    
    def on_update_plot_clicked(self):
        """Update existing plot with new settings"""
        settings = self.view.get_all_settings()
        
        if not settings['hold_plots']:
            self.view.clear_plot()
        
        self._apply_settings(settings)
        
        if self.model.regenerate_spectrum():
            self._plot_spectrum()
            self._load_experiment_if_available(settings)
            self._apply_plot_formatting(settings)
    
    def on_peak_shape_changed(self, shape_type, enabled):
        """Handle peak shape checkbox changes"""
        if shape_type == 'gaussian':
            self.model.update_plot_setting('gaussian_enabled', enabled)
            self.model.update_plot_setting('lorentzian_enabled', False)
        elif shape_type == 'lorentzian':
            self.model.update_plot_setting('lorentzian_enabled', enabled)
            self.model.update_plot_setting('gaussian_enabled', False)
    
    def on_save_figure_clicked(self):
        """Save figure to file"""
        theory_file = self.model.get_theory_file_path()
        if not theory_file:
            self.view.show_error("No theory file selected")
            return
        
        try:
            import matplotlib.pyplot as plt
            output_path = theory_file.replace('.out', '.png').replace('.log', '.png')
            plt.savefig(output_path, dpi=500)
            self.view.display_message(f"Saved figure to {output_path}")
        except Exception as e:
            self.view.show_error(f"Failed to save figure: {str(e)}")
    
    def on_export_data_clicked(self):
        """Export spectrum data"""
        success, message = self.model.export_spectrum()
        if success:
            self.view.display_message(message)
        else:
            self.view.show_error(message)
    
    # Private helpers
    def _apply_settings(self, settings):
        """Update model with current view settings"""
        for key, value in settings.items():
            self.model.update_plot_setting(key, value)
    
    def _plot_spectrum(self):
        """Plot the spectrum on canvas"""
        from PlotManager import PlotManager
        
        spectrum = self.model.spectrum
        sticks = self.model.sticks
        plot_bars = self.model.get_plot_setting('plot_bars')
        
        PlotManager.plotTheoryData(self.ax, spectrum, sticks, int(plot_bars))
        self.view.update_canvas()
    
    def _load_experiment_if_available(self, settings):
        """Load experimental data if available"""
        exp_path = self.model.get_experiment_file_path()
        if exp_path:
            try:
                from PlotManager import PlotManager
                scale = settings['scale_exp']
                PlotManager.plotExperimentData(exp_path, scale)
            except Exception as e:
                self.view.display_message(f"Warning: Could not load experiment data: {str(e)}")
    
    def _apply_plot_formatting(self, settings):
        """Apply formatting to plot"""
        from matplotlib.ticker import AutoMinorLocator
        
        ymax = settings['ymax']
        if ymax != 'auto':
            self.ax.set_ylim(ymin=0, ymax=float(ymax))
        
        self.ax.set_xlim(xmin=settings['xmin'], xmax=settings['xmax'])
        self.figure.subplots_adjust(top=0.92, right=0.95, bottom=0.16, left=0.15)
        self.ax.set_xlabel('wavelength / nm', size=14)
        self.ax.set_ylabel('normalized absorbance', size=14)
        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.legend(edgecolor='none')
        self.ax.tick_params(axis='both', labelsize=14)
        self.figure.set(figheight=3.8, figwidth=7, dpi=100)
        
        self.view.update_canvas()
    
    # Observer interface
    def on_model_changed(self, event_type, data):
        """Handle model change notifications"""
        if event_type == 'error':
            self.view.show_error(data['message'])
        elif event_type == 'theory_file_loaded':
            self.view.display_message(f"Loaded {data['num_states']} excitation states")
        elif event_type == 'spectrum_exported':
            self.view.display_message(f"Spectrum exported successfully")