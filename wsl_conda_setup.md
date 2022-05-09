Set up wsl (windows subsystem for linux): 
https://docs.microsoft.com/en-us/windows/wsl/install

Set up miniconda: 
Go to https://docs.conda.io/en/latest/miniconda.html#linux-installers and hover over 
appropriate Linux installer link. Copy the link address `<miniconda linux installer>`
	
	$ wget <miniconda linux installer>
	$ bash Miniconda[YOUR VERSION].sh

Restart your WSL terminal

In a WSL terminal:
	
	conda update -n base -c defaults conda
	conda create --name docs_env
	conda activate docs_env
	conda install jupyter
	conda install -c r r-irkernel
	jupyter notebook --no-browser


