"""DaltonView 2.0 - Application entry point"""
import tkinter as tk
from model import SpectrumModel
from view import SpectrumView
from controller import SpectrumController

if __name__ == "__main__":
    # set MVC components
    model = SpectrumModel()
    view = SpectrumView(controller=None)  # placeholder controller
    controller = SpectrumController(model, view)
    view.controller = controller  # set controller after creation
    
    # run application
    view.mainloop()