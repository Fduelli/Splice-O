import subprocess
#GTCTCACTGTGTTCTCTG
#GAGGGGGTGGGGGTGGGA
#hCLN3 ASO 6
#CCTGGAAGCTCTGCGGTC
#hCLN3 ASO 25
#TCCTCTTGGCCTTCACCT
# Run the Java JAR file and wait for it to finish
subprocess.run(['java', '-Xms512m', '-Xmx4g', '-jar', '/home/fduelli/customASO_project/ASO_walk_gen.jar'])

# If the above process completes, run the Python file in the same directory
filename = 'runspliceai_w_customASO.py'  # Replace with the name of your Python file
subprocess.run(['python', filename])

