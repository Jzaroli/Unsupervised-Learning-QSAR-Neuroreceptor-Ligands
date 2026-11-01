# Unsupervised Learning: QSAR & Neuroreceptor Ligands
  
## Project Topic:  
The goal for today is to run unsupervised learning algorithms on a dataset comprised of small molecules known to bind to dopmainergi, serotonergioc, gabaergic, and histaminergic neuroes. Their chemical SMILES will be converted into various features and which will then be used to to classify the molecules back into their target receptor classes. The project uses both PCA and K-Means.

## Data Overview:

Database: ChEMBL is a database with a trove of cheminformatics related data:
  
https://www.ebi.ac.uk/chembl/ 
  
It self describes as "ChEMBL is a manually curated database of bioactive molecules with drug-like properties. It brings together chemical, bioactivity and genomic data to aid the translation of genomic information into effective new drugs."
  
Some chem domains and ML applications have low data and high dimensions. It can be incredibly useful to run clustering or PCA algorithms on the data to understand underlying patterns, find the right features, group chemicals based on underlying physical structure or bioactivity, and/or overall do prep work for downstream cheminformatics work.  

## Fetching data from ChEMBL:  
  
### Definitions (Wikipedia):    
* Ki, the dissociation constant, "describes the binding affinity of the inhibitor for the target (the concentration of free inhibitor needed to occupy 50% of the target's binding sites at equilibrium)."  
* IC50 is defined as "the concentration of an inhibitor required to reduce a specific biological or biochemical activity by 50% under defined assay conditions"
* The Simplified Molecular Input Line Entry System (SMILES) "is a specification in the form of a line notation for describing the structure of chemical species using short ASCII strings."
* In drug discovery and pharmacology, pChEMBL is a standardized, logarithmic scale used to measure how potently a small molecule binds to a biological target (such as a receptor or enzyme).
Because raw binding affinity values like $IC_{50}$ or $K_i$ span massive orders of magnitude (from millimolar down to picomolar), working with them directly can skew models. Taking the negative logarithm transforms them into an intuitive scale where higher numbers mean stronger potency.
  * The pChEMBL value is calculated as the negative base-10 logarithm of the molar concentration ($M$): $$\text{pChEMBL} = -\log_{10}(\text{Concentration in Moles/Liter})$$
  * In high-throughput screening and database curation, compounds requiring more than $1\,\mu M$ of concentration to show an effect are often considered weak binders, noise, or inactive background molecules.
  
| Concentration (Molar) | Concentration (Standard Units) | pChEMBL Value | Interpretation |
| :--- | :--- | :--- | :--- |
| $1.0\text{ M}$ | $1,000,000\,\mu M$ | **0** | Extremely weak / inactive |
| $1.0 \times 10^{-3}\text{ M}$ | $1,000\,\mu M$ ($1\text{ mM}$) | **3** | Very weak |
| $1.0 \times 10^{-6}\text{ M}$ | **$1\,\mu M$** | **6** | Active threshold |
| $1.0 \times 10^{-7}\text{ M}$ | $100\text{ nM}$ | **7** | Strong / Potent |
| $1.0 \times 10^{-8}\text{ M}$ | $10\text{ nM}$ | **8** | Very potent |
| $1.0 \times 10^{-9}\text{ M}$ | $1\text{ nM}$ | **9** | Extremely potent |
  
## Feature Engineering  
  
https://www.rdkit.org/docs/source/rdkit.Chem.html
https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
  
Method Descriptions:  
* MolFromSmiles: Construct a molecule from a SMILES string.
* MolWt: The exact molecular weight of the molecule.
* MolLogP: Atom-based calculation of LogP and MR using Crippen’s approach
* NumHDonors: The number of hydrogen bond donors in the molecule, typically calculated based on the number of heteroatoms attached to at least one hydrogen.
* NumHAcceptors:The number of hydrogen bond acceptors in the molecule, representing the count of heteroatoms (like nitrogen or oxygen) that can accept a hydrogen bond.
* TPSA: Topological Polar Surface Area, which estimates the surface area of the molecule's polar atoms (usually oxygen, nitrogen, and attached hydrogens) to predict drug transport and absorption properties.
* NumRotatableBonds: The number of single bonds that are not in a ring and can rotate freely, which is a key measure of molecular flexibility.
* RingCount: The total number of rings present in the molecule.
* FractionCSP3: The fraction of carbon atoms that are sp3 hybridized, used as an indicator of molecular three-dimensionality and complexity.
* MolMR: The molecular refractivity calculated using Crippen’s atom-based approach, which relates to the molecule's polarizability and bulk volume.

## Models:  

### PCA:  

#### Explained Variance Ratio Results:  
* Explained variance ratio sum (3 components): 0.8042
* Explained variance ratio sum (4 components): 0.8879
* Explained variance ratio sum (5 components): 0.9522
* Explained variance ratio sum (6 components): 0.9739
* Explained variance ratio sum (7 components): 0.9880
* Explained variance ratio sum (8 components): 0.9971
  
I proceeded with 4 principle components to keep dimensionality down for K-means while capturing 89% of the variance.  
  
#### PCA Results:  
  
| | MolWt | LogP | NumHDonors | NumHAcceptors | TPSA | NumRotatableBonds | RingCount | FractionCSP3 | MolRefractory |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **PC1** | 0.482867 | 0.220128 | 0.115861 | 0.367124 | 0.309662 | 0.401020 | 0.285417 | 0.082998 | 0.474468 |
| **PC2** | -0.098957 | -0.445129 | 0.540605 | 0.292315 | 0.543923 | -0.063309 | -0.240448 | -0.155059 | -0.180684 |
| **PC3** | 0.000883 | -0.435801 | -0.342324 | 0.230151 | -0.039181 | 0.242651 | -0.267760 | 0.711692 | -0.036140 |
| **PC4** | -0.023234 | -0.370437 | -0.183230 | 0.147726 | 0.095989 | -0.470683 | 0.754035 | 0.086357 | -0.007576 |

#### PCA Analysis:  
* PC1: This component is dominated by molecular weight and the molecule's polarizability and bulk volume, meaning large vs small molecules.
* PC2: This component is dominated by number of hydrogen donors and TPSA which estimates the surface area of the molecule's polar atoms, meaning hydrophilic vs lipophylic molecules
* PC3: This component is dominated by FractionCSP3, which is used as an indicator of molecular three-dimensionality and complexity. We can interpret this as complex three-dimensional molecules vs more flat ring-like structures. 
* PC4: This component is dominated by ring count vs number of rotatable bonds which is similar to PC3 in that on one end we have ring-based molecules but on the other we now have flexible configurations arising from the number of rotatable bonds

### K-Means Clustering:  

#### K-Means Evaluation Table:
  
| receptor_class | 0 | 1 | 2 | 3 |
| :--- | :---: | :---: | :---: | :---: |
| Dopamine_D2 | 806 | 71 | 26 | 92 |
| GABA_A | 63 | 586 | 1 | 350 |
| Histamine_H1 | 107 | 477 | 0 | 416 |
| Serotonin_5HT2A | 415 | 109 | 0 | 476 |
  
<b> Adjusted Rand Index (ARI): 0.1933 </b>  

## Observations and Conclusions:  
It wasn't possible to segment the data back into the four receptor classes based on the four PCA components. The k-means algorithm did not segment based on neuromolecular target but based on structural (dis)similarities.
  
* Cluster 0 represents the molecular weight/size component which shows that dopamine class segments well into this cluster. 
* The GABA class is more highly clusterd in PCA 2 and 4 for hydrophilicity vs lipophylicity and ring count vs rotatable bonds.
* The FractionCSP3 dominated cluster has very few candidates total and of the data present in this cluster, a vast majority of dopamine ligand molecules, suggesting that higher molecular weight/size correlates with three-dimensional complexity in this class.
* Finally, the fourth cluster which includes the ring count vs rotatable bonds component shows a strong presence from the GABA, histamine and serotonin class.
  
The Adjusted Rand Index (ARI) of 0.193 indicates a low degree of alignment between the unsupervised K-Means clusters and the true neuroreceptor classes. This occurs because the structural descriptor features group molecules by structural sub-groups rather than neuroreceptor classes.
  
## Future Directions:  
It would be interesting to review the generalities of the biology and the chemical structures of the ligands of each neuroreceptor class to then compare them to the results. It would also be interesting to run various algorithms using Morgan Fingerprints which can turn smaller molecular and chemical substructures into features!

## Additional Links and Documentation:
https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html  
https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html  
