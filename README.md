# DaltonView
## TD-DFT spectrum viewing software

An executable version is available at https://daltonian.co/daltonview or this [direct download link](http://u.pc.cd/CjHotalK)

Just a few tips on how the file system is configured:
1. Output type is determined by the first line of the output file. Gaussian and Q-Chem print their respective names in the first line of an output file by default. If your implementation differs from this, you could just add 'Gaussian' or 'Q-Chem' to the first line of the output.
2. The experimental file should be configured as column data, with wavelength in the first column and absorbance as the second column. Headers are okay, the first row is skipped to allow for this.
3. Exporting will automatically add descriptors to the original output file name and save them in the same folder.



<img src="https://daltonian.co/images/daltonview_screenshot.png"
     alt="Program Screenshot"
     style="float: middle;" />
