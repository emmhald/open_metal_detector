# A Python program to analyze collections of Metal Organic Frameworks (MOFs) for open metal sites.
Author: Emmanuel Haldoupis <emmhald@gmail.com>

March 12th 2018

## Description

Given a set of CIF files this program can read the files, analyze the structures and
detect all the unique open metal sites present in each structure.

The first step is to create a collection containing all the desired CIF files. We can do
this by pointing to a folder containing the CIF files.

```
from omsdetector import MofCollection
mof_coll = MofCollection.from_folder(collection_folder="path to cif folder",
                                     analysis_folder="path to analysis folder")
```

Where **collection_folder** is the folder where the CIF files are located and **analysis_folder** is the folder where the results will be saved.


The analysis is run using the following command on the mof_coll object:

```
mof_coll.analyse_mofs()
```

Specifying a value for num_batches instructs the analysis to run in parallel in the specified number
of batches, each as a separate process.

Once the results have finished they can be summarized using the following methods:

```
mof_coll.summarize_results()
mof_coll.summarize_tfactors()
```

The summarize_results() method generates a table that summarizes the number of open metal sites found for each metal type.
The summarize_tfactors() method generates histograms (and stores them) for the distribution of the t-factors, which indicate
the degree of deviation from a closed coordination sphere for tetra, penta, and hexa-coordinated coordination spheres.

Finaly, a collection can be filtered to create a sub-collection using the following filters:

* "density": [min, max] (range of values)
* "oms_density": [min, max] (range of values)
* "uc_volume":  [min, max] (range of values)
* "metal_species": ["Cu", "Zn", ...] (list of metal species)
* "non_metal_species": ["C", "N", ...] (list of non metal species)
* "cif_okay": True (boolean value)
* "has_oms": True (boolean value)
* "mof_name":  [mof_name1, mof_name2] (string values)

For example:

```
co_oms = mof_coll.filter_collection(using_filter={"metal_species":["Co"], "has_oms":True})
```
See the example jupyter notebook for more details.

## Requirments
* python 3.6.3
* pymatgen  2018.2.13
* pandas 0.22.0
* numpy 1.14.1
* matplotlib 2.1.1

## Reference

Chung, Yongchul; Haldoupis, Emmanuel; Bucior, Benjamin; Haranczyk, Maciej ; Zhang, Hongda; Vogiatzis, Konstantinos; Milisavljevic, Marija; Ling, Sanliang; Camp, Jeffrey; Slater, Ben; Siepmann, J.; Sholl, David; Snurr, Randall *Computation-Ready, Experimental Metal-Organic Framework Database 2018: Additional Structures, Open Metal Sites, and Crystal Reconstruction* (submited to Chemistry of Materials)



# March 1, 2022 Update

Author: Grier M. Jones
University of Tennessee,Knoxville


Updated March 1, 2022
Updates:
- Updated pymatgen dependencies
- Updated reference
## Reference
Chung, Y. G.; Haldoupis, E.; Bucior, B. J.; Haranczyk, M.; Lee, S.; Zhang, H.; Vogiatzis, K. D.; Milisavljevic, M.; Ling, S.; Camp, J. S.; et al. Advances, Updates, and Analytics for the Computation-Ready, Experimental Metal-Organic Framework Database: CoRE MOF 2019. J. Chem. Eng. Data 2019, 64 (12), 5985–5998.



# Easy install for jupyter
Create conda environment, link it to a jupyter kernel, and install required packages
```
conda create --name omd
conda activate omd
ipython kernel install --user --name=omd
python3 -m pip install -r requirements.txt
python3 setup.py install
```
