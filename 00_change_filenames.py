
source activate CUTRUN

python
import os
import glob
files = glob.glob('*.fastq.gz') #the path could specified instead
for file in files:
parts = file.split('_') 
new_name = ''.join([parts[0], '_', parts[1], '_', parts[3], '_', parts[4]])
	ext = '.fastq.gz'
	new_file = ''.join([new_name, ext]) 
	print(new_file) # run this one just to check is correct before next step
os.rename(file, new_file)
