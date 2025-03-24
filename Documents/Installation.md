# Installation guide

All the required libraries are included in the [requirements.txt](../requirements.txt) file. Install them globaly or through an [enviroment](#creating-an-enviroment). 

If using the IAA server (recomended to be able to use images IDs) the use of an enviroment is mandatory. 

## Creating an enviroment. 

- Anywhere on your directories -> run the command: 
```shell
"python3 -m venv $myenv$"
```
Substituting "$myenv$" for whatever name you want to give it. 
- Activate the enviroment: run the command: 
```shell
source $myenv$/bin/activate
```
Again, substitute the name of "$myenv$". 

If the previous command raises an error try: 
```shell
source $myenv$/bin/activate.csh
```

## Install dependencies: 
Run in terminal:
```shell
pip install -r requirements.txt
```

You should be able to run any script after the installation finishes with the enviroment activated. 