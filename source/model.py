"""Model: Encapsulates all data and business logic"""
import numpy as np
import re
from FileManager import FileManager
from DataManager import DataManager

class SpectrumModel:
    """Manages all spectral data and computations"""
    
    def __init__(self):
        self._file_manager = FileManager()
        self._data_manager = DataManager()
        
        # Application state
        self._theory_file_path = ""
        self._experiment_file_path = ""
        self._observers = []
        
        # Spectral data
        self.eV = None
        self.nm = None
        self.OS = None
        self.spectrum = None
        self.sticks = None
        
        # Plot settings
        self._plot_settings = {
            'xmin': 100,
            'xmax': 2000,
            'ymax': 'auto',
            'hwhm': 200,
            'gaussian_enabled': True,
            'lorentzian_enabled': False,
            'plot_bars': True,
            'hold_plots': False,
            'scale_exp': 1.0,
        }
    
    def attach_observer(self, observer):
        """Observer pattern for notifying View of changes"""
        self._observers.append(observer)
    
    def notify_observers(self, event_type, data=None):
        """Notify all observers of state changes"""
        for observer in self._observers:
            observer.on_model_changed(event_type, data)
    
    # File operations
    def set_theory_file(self, path):
        """Load theory file and extract spectral data"""
        self._theory_file_path = path
        try:
            lines = self._file_manager.openTheoryFile(path)
            self.eV, self.nm, self.OS = self._data_manager.retrieveDataFromLines(lines)
            self.notify_observers('theory_file_loaded', {
                'path': path,
                'num_states': len(self.eV) if self.eV else 0
            })
            return True
        except Exception as e:
            self.notify_observers('error', {'message': str(e)})
            return False
    
    def set_experiment_file(self, path):
        """Set experimental data file path"""
        self._experiment_file_path = path
        self.notify_observers('experiment_file_set', {'path': path})
    
    def get_theory_file_path(self):
        return self._theory_file_path
    
    def get_experiment_file_path(self):
        return self._experiment_file_path
    
    # Plot settings
    def update_plot_setting(self, key, value):
        """Update a plot setting and notify observers"""
        if key in self._plot_settings:
            self._plot_settings[key] = value
            self.notify_observers('setting_changed', {key: value})
            return True
        return False
    
    def get_plot_setting(self, key):
        return self._plot_settings.get(key)
    
    def get_all_plot_settings(self):
        return self._plot_settings.copy()
    
    # Spectrum generation
    def regenerate_spectrum(self):
        """Regenerate spectrum with current settings"""
        if self.eV is None or self.nm is None or self.OS is None:
            self.notify_observers('error', {'message': 'No theory file loaded'})
            return False
        
        try:
            hwhm = self._plot_settings['hwhm']
            gaussian = self._plot_settings['gaussian_enabled']
            
            self.spectrum, self.sticks = self._data_manager.makeSpectrumFromRetrievedData(
                eV=self.eV,
                nm=self.nm,
                OS=self.OS,
                hwhm=hwhm,
                peakShapeBool=int(gaussian)
            )
            self.notify_observers('spectrum_updated', {
                'spectrum': self.spectrum,
                'sticks': self.sticks
            })
            return True
        except Exception as e:
            self.notify_observers('error', {'message': f'Spectrum generation failed: {str(e)}'})
            return False
    
    def export_spectrum(self):
        """Export spectrum data to CSV files"""
        if self.spectrum is None or self.sticks is None:
            return False, "No spectrum to export"
        
        if not self._theory_file_path:
            return False, "No theory file loaded"
        
        try:
            base_filename = self._theory_file_path.replace('"', '').split('.')[0]
            np.savetxt(base_filename + '_spectrum.csv', self.spectrum.T, 
                      delimiter=',', header="nm,abs")
            np.savetxt(base_filename + '_sticks.csv', self.sticks.T, 
                      delimiter=',', header='nm,f')
            self.notify_observers('spectrum_exported', {'path': base_filename})
            return True, f"Exported to {base_filename}"
        except Exception as e:
            return False, str(e)