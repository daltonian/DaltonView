class LogManager:
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