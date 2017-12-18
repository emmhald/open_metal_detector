# Python program to analyze collections of MOF structures for open metal sites.
 Author: Emmanuel Haldoupis <emmhald@gmail.com>

 09/18/2017
##Description
Given a set of cif files this program can read all the files and analyze the structures and
detect all the unique open metal sites present in each structure.

This creates a collection, where collection_folder is where the cif files are and
output_folder where the results will be saved.
~~~
from mof_collection import MofCollection 
mof_coll= MofCollection(collection_folder=collection_folder, output_folder=output_folder)
~~~

The analysis is run using the following command on the mof_coll object
~~~
mof_coll.analyse_mofs(num_batches=3)
~~~
Specifying a value for num_batches instructs the analysis to run in parallel in the specified number
of batches, each as a separate process.

Once the results have finished they can be loaded and summarized:
~~~
mof_coll.load_results()
mof_coll.collect_statistics()
mof_coll.summarize_tfactors()
~~~

The collect_statistics generates a table that summarizes the number of open metal sites found for each metal atom.
The summarize_tfactors, generates histograms (and stores them) for the distribution of the t-factors, which indicate
the degree of deviation from a closed coordination sphere for tetra, penta and hexa-coordinated coordination spheres. 

See the example jupyter notebook for more details.

##Requirments
pymatgen
pandas
numpy
matplotlib
