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