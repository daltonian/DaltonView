class DataManager:
    def __init__(self):
        self.counts = 0
    
    def retrieveDataFromLines(self,lines):
        """ Retrieves excitation data from quantum chemical output file 
            Returns list with eV, nm, and Oscillator Strength values """
        #...
        # Find lines which contain eV, oscialltor strength, and FC factor

        if 'Q-Chem' in lines[0]:
            try:
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
            except ValueError:
                print("This does not appear to be a TD-DFT file.\nDouble check selection")
                return
                           
        
        elif 'Gaussian' in lines[0]:
            try:
                # Find lines which contain eV, oscialltor strength, and FC factor
                start = [x for x in range(len(lines)) if 'Excitation energies and oscillator strengths:' in lines[x]][0]
                end = [x for x in range(len(lines)) if ' Population analysis using the SCF density.' in lines[x]][0]
                state_index = [x for x in range(start,end) if 'Excited State  ' in lines[x]]        # Find which lines have eV values

                eV = []
                nm = []
                OS = []
                for idx in state_index:
                    eV.append(float(re.findall('\d*\.?\d+',lines[idx])[1]))
                    nm.append(float(re.findall('\d*\.?\d+',lines[idx])[2]))
                    OS.append(float(re.findall('\d*\.?\d+',lines[idx])[3]))
            except ValueError:
                print("This does not appear to be a TD-DFT file.\nDouble check selection")
                return

        elif 'GAMESS' in lines[0]:
            try:
                # Find lines which contain eV, oscialltor strength, and FC factor
                start = [x for x in range(len(lines)) if 'R-TDDFT CALCULATION CONVERGED' in lines[x]][0]
                end = [x for x in range(len(lines)) if 'SUMMARY OF TDDFT RESULTS' in lines[x]][0]
                state_index = [x for x in range(start,end) if 'STATE # ' in lines[x]]        # Find which lines have eV values

                eV = []
                nm = []
                OS = []
                for idx in state_index:
                    eV.append(float(re.findall('\d*\.?\d+',lines[idx])[0]))
                    nm.append(1239.8/eV[-1])
                    OS.append(float(re.findall('\d*\.?\d+',lines[idx+1])[0]))
            except ValueError:
                print("This does not appear to be a TD-DFT file.\nDouble check selection")
        
        else:
            print("Could not determine input file type...\nVerify that the name of the program of origin is in the first line.")
            return

        return eV, nm, OS

    def makeSpectrumFromRetrievedData(self,eV,nm,OS,hwhm,peakShapeBool):
        # Define wavelength (or eV) region over which to create spectrum
        # # Position of wavelengths will effect the normalization of the
        # # spectrum, ie should ideally match experimental ref range
        min_eV = 0.4 # Minimum energy in eV (~3100 nm)
        max_eV = 15 # Maximum wavelength in eV (~80 nm)
        num_eV = int(3000*2+1)
        space = np.linspace(min_eV,max_eV,num_eV)

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