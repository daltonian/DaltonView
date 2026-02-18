class PlotManager:
    def __init__(self):
        self.counts = 0

    @staticmethod
    def plotTheoryData(ax,spectrum,sticks,barsBool):
        #... 
        ax.plot(spectrum[0],spectrum[1],linewidth=2,label="Theory Spectrum")
        if barsBool == 1:
            for x in range(len(sticks[0])):
                print(sticks[0,x],sticks[1,x])
            ax.bar(sticks[0],sticks[1],width=2)
        return
    
    @staticmethod
    def plotExperimentData(path_to_experimental_csv,scale_factor):
        experimental_data = np.loadtxt(path_to_experimental_csv,delimiter=',',skiprows=1)
        ax.plot(experimental_data[:,0],experimental_data[:,1]*float(scale_factor),linewidth=2,color='k',label='Experimental Spectrum')
        return