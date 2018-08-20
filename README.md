1.  clone the env/environment.yml file via conda 
Instructions are here:
https://conda.io/docs/user-guide/tasks/manage-environments.html#cloning-an-environment  

2.  Adjust the settings in the config.ini file for:  
-  AA_unit = The repeat unit size in Amino acids (AA)
-  AA_skipNterm = number of AAs to skip at the START of the sequences  
-  AA_skipCterm = number of AAs to skip at the END of the sequences  

These AA units will be multiplied by 3 by the script to adjust the DNA
translaton.  

3.  Run the script specifying the infile in fasta format.  
```python repeat_plotter.py INFILE.fasta```  

4.  The resulting outfiles are:  
-  ```outfile.pdf``` is the figure as a pdf  
-  ```outfile.tab``` is a text output of the repeat pattern.  
