# PDA-Pred
Protein-DNA binding affinity prediction

PDA-Pred predicts protein-DNA binding affinity using structure-based features and classes based on structure, function and percentage of binding site of protein. Generally, contact potenitals, accessible surface area, interactions between various atoms and their energy contributions, base step parameters of the DNA, secondary strucutres and residue depth of protein are important to understand the binding affinity. PDA-Pred shows a correlation of 0.86 and a mean absolute error(MAE) of 0.76 in Self-consistency, correlation of 0.78 and a MAE of 0.98 kcal/mol in jack-knife test.

This program takes the input of protein-DNA complex and predicts the binding affinity.

Input options:
User can input PDB ID of protein-DNA complex or can provide the file in PDB format. *
User can enter the chain ID of protein. Multiple chain Ids can be entered with comma. (Eg. Chain ID: A,B) (Optional)
User can enter the chain ID of DNA. Multiple chain Ids can be entered with comma. (Optional)
User should select one of the DNA class for the prediction. *
1. Single strand
2. Double strand
User should select one of the strucutral class of protein for the prediction, if the DNA is double stranded *
1. all-α
2. all-β
3. αβ
4. others
User should select one of the Functional class of protein for the prediction. *
1. Regulatory
2. other (non-regulatory)

Dependencies:
1. Program is basically using python3, and demands few python packages to run.
2. Please ensure the packages such os,re, Bio, numpy, functools,sys, time, shutil, subprocess, math, uuid, cgitb, timeit, wget, glob, urllib, pandas, warnings are installed.
3. Make sure that the naccess, has the executable path, installed in the local system
4. Also it need the installation of 3vvv software, used for volume and surface area calculation (Source: https://github.com/vosslab/vossvolvox)
5. Foldx exectable file

Run:
1. Run "python3 pdpredict.py"
2. Enter option for PDB-ID and PDB-file
3. Provide the PDB ID or PDB-file with chain
4. Enter the DNA strand, strucutral and functional class
5. The program will be executed, and the result will be provided in result.txt
In the Form of PDB ID, Chain, Predicted affinity(kcal/mol), dissocation constant(M).
You can also access it through the webserver https://web.iitm.ac.in/bioinfo2/pdapred/index.html
