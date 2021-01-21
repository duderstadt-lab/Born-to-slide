<p><img src="https://github.com/mjosch/Born-to-slide/blob/master/graphical_abstract.png" width=“800"></p>

Analysis software package used in:

### Born to slide: mobile origin licensing factors confer resistance to transcription conflicts
Matthias J Scherr, Syafiq Abd Wahab, Dirk Remus, Karl E Duderstadt  
published @xxx (doi:xxx)

Getting started:

1) Install Fiji (https://imagej.net/Fiji/Downloads) and add MARS to your update sites (http://sites.imagej.net/Mars/).  
A detailed installation guide can be found here: https://duderstadt-lab.github.io/mars-docs/install/. 
All MARS source code is publicly available in several repositories on Github at https://github.com/duderstadt-lab. 
The core library used for analysis and storage of data is contained in the mars-core repository. 
The graphical user interface is contained in the mars-fx repository. Documentation can be found at mars-docs.

2) Replace the following jars in your Fiji installation with the jars provided in https://github.com/mjosch/Born-to-slide/tree/master/MARS_jars
- mars-autocompletion
- mars-commands
- mars-core
- mars-fmt
- mars-fx
- mars-swing
- mars-trackmate

3) Create a local python environment with **born-to-slide.yml**

4) Download all datasets (organized in archives) from XXX (raw videos are deposited under xxx)

5) Download all python scripts and jupyter notebooks

6) Set the location of the dataset folder to the **filepath** variable in **awesome_data.py**

7) Run jupyter notebooks to check out the data in each figure panel
