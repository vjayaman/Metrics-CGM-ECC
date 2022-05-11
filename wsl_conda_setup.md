Set up wsl (windows subsystem for linux): 
https://docs.microsoft.com/en-us/windows/wsl/install

Set up miniconda: 
Go to https://docs.conda.io/en/latest/miniconda.html#linux-installers and hover over 
appropriate Linux installer link. Copy the link address `<miniconda linux installer>`int a Windows PowerShell terminal
	
	$ wget <miniconda linux installer>
	$ bash Miniconda[YOUR VERSION].sh

Restart your WSL terminal

In a WSL terminal:
	
	conda update -n base -c defaults conda

We are going to use the conda-forge channel (for most recent R version, etc.)

        conda config --add channels conda-forge
        conda config --set channel_priority strict

Then set up the environment *nb_env* using the YAML file in the documentation directory

	conda env create -f documentation/environment.yml
	conda activate nb_env

Run the following statement(s) (so we can use the `Rcpp` package)

        conda install -c conda-forge gxx

In an R environment, run: 

	install.packages("dint")

If going through the jupyter notebook documentation, run through the following in a WSL terminal:

	jupyter notebook --no-browser


Note that the whole process may take some time.  
