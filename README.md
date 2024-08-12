# MPC2008-Public
Replication code for "Micro MPCs and Macro Counterfactuals: The Case of the 2008 Rebates" by Jacob Orchard, Valerie Ramey, and Johannes Wieland


## Permissions

You can use our code with proper attribution.

Please cite as:

Orchard, Jacob, Valerie A. Ramey, and Johannes F. Wieland. Micro MPCs and macro counterfactuals: the case of the 2008 rebates. The Quarterly Journal of Economics (Forthcoming)


# To Run Entire Project

You will need your own FREDKEY and BEA keys to download the source data. Place the FREDKEY in line 34 of `MPC/forecasting/code/build_forecast_data.do` and place the BEA key in line 14 of `MPC/downloaddata/code/pcefromBEA.py`. 

We used STATA version 16.1 and Python 3.11 to run this project.

### UNIX and MAC users
UNIX and MAC users can run the entire project using `make`. Simply type the following three commands in a terminal and the project will build from scratch. 

(1) `make install`

(2) `make venv`

(3) `make`

Once `make` executes successfully, the paper figures and tables are available in the folder `_finaltablesandfigures/output`.

The Monte Carlo simulations take XX hours to complete so we do not include them as part of the baseline replication flow. To run the Monte Carlo simulations execute this command:

(4) `make montecarlo`

The Monte Carlo figures will be produced in the folder `montecarlo/output`.

# To Run only the Empirical Regressions

An archived version of the data used to produce the main empirical regression tables (Tables 3-5) is stored in psmjregressions/output/archiveinteriew.zip. After unzipping the file, the main regression tables can be reproduced using Table_E1.do (table 3), Table_bias.do (Table 4), and Table_E2.do (Table 5). All of those do files are stored in psmjregressions/code.

For non-STATA users, a csv of the regression dataset is stored in psmjregressions/output/psmjsampleinteriew_wlabels.csv.

# Order of Tasks to Create Final Output

This project is divided into a series of subfolders that execute all of the tasks leading to final output beginning with downloaddata and ending with _finaltablesandfigures. Each subfolder contains both a code directory and, once-executed,  input and output directories. The `makefile` in the code folder documents how the inputs are converted in the outputs for the task. The input directory will have symbolic links to output from previous tasks, while the output directroy will include all of the output used by subsequent tasks. 

The makefile, "make" in the main folder shows the order of execution of the subfolders. The final output for the paper is mostly created in the forecasting, psmjregressions, model, and narrative subfolders. 


# How to Run MPC project on a Windows Computer

The replication code files are designed to be run on a server or computer that has access to the `make` command. We make no guarantees that the project will run smoothly on a Windows system, since `make` is not native to Windows. However, we have been able to run the project on Windows by installing a version of `make` using chocolatey. These are the steps that worked for us, however, we do not maintain any of the linked how-to sites or packages and cannot provide support for these steps. 

Alternatively, a windows user can run separately the code files they need in our replication files. 

## Preliminaries

Prior to running the project, you'll want to make sure that you have the latest versions of Anaconda and Git set up on your windows PC. Then install make via chocolatey. We'd recommend installing make last.

### Git and Git Bash (Required)

Follow the installation instructions [here to install Git](https://www.computerhope.com/issues/ch001927.htm). Git comes bundled with **Git Bash**, an alternative to windows command prompt. We will use Git Bash later on, as it can utilize some unix style commands. 

### Github Account and personal access token (Required)

Once Git is installed. You'll also need to create an account on Github. Once you create an account, [follow these steps](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) to create a personal access token. You'll use this personal access token when downloading and updating the project. 

### Anaconda (Required)

You may already have Anaconda installed on your system. If not, you can follow the guide [here](https://problemsolvingwithpython.com/01-Orientation/01.03-Installing-Anaconda-on-Windows/)

Once installed, we'd recommend[ adding Anaconda to your path.](https://www.geeksforgeeks.org/how-to-setup-anaconda-path-to-environment-variable/)


### Make for Windows (Required)

There are a number of ways to install make on windows. We've tried most and would highly recommend installing make via chocolatey, which is a package manager for windows (similar to anaconda, but for applications outside of python). 

[This article ](https://pakstech.com/blog/make-windows/) has step by step instructions on how to install chocolatey and then how to install make from chocolatey. Chocolatey should automatically add make to your path. 

### Microsoft Visual Studio Code (Optional)

This is completely optional, but Visual Studio Code (VScode) integrates well with git, git bash, and python. [You can download it here.](https://code.visualstudio.com/download)

There is a git extension on Visual Studio called [Git Lens](https://marketplace.visualstudio.com/items?itemName=eamodio.gitlens) that you can download via the extension manager. Git has a bit of a steep learning curve, and git lens allows you to see what git is doing in a GUI rather than having to push and pull via the command line. 

## Setting up the project

### Step 1: Ensure all applications are running

Open up **Git Bash**. Type the following commands into the terminal to make sure everything from the preliminary step is working correctly. Everything should be working fine, if you get a response from each command.

`conda --version`

`git --version`

`make --version`

### Step 2: Download project via Git

Navigate to the folder where you want to place the project (this should not be inside of dropbox since git and dropbox can clash)

Next, clone the project using the following line in **Git Bash**

`git clone url`

Where url is the project url you can find on Github. Git may prompt you to enter your username and password. Instead of using your Github username and password, [use the personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) you set up earlier. 

When this step is successful, you'll see all of the project code within the MPC folder, with the exception of the code for the did_imputation method. The did_imputation folder will be empty since this folder links to Boryusak's Github, not ours. In order to pull his most recent code run the following:

`git submodule init`
`git submodule update`

### Step 3: Setup the Conda Virtual Environment

There are a lot of versions of python and many many python packages with their own versions as well. In order to ensure that our code works as intended, it's important to make sure we are using the same version of these packages. That is what this step does.

Within **Git Bash** start by initiating Anaconda. 

`conda init`

You'll then need to close **Git Bash** and open it again.

Next, with **Git Bash** reopened, navigate to the MPC directory and type in the following:

`conda env create -f win_environment.yml`

This will create a _virtual environment_ that will have all of the same packages that we used when running the code on my windows machine. This environment is called MPC. It may take 10-20 minutes for conda to install all of the necessary packages (you also may be prompted to type **y** to agree to the setup). 

Activate the environment:

`conda activate MPC`

### Step 4: Run the project

Type the following into **Git Bash**

`make windows`

The entire project should run. 

FYI, the createconsumptionvariables step requires a large amount of RAM. You can make this step go quicker by editing

`MPC\globaloptions\hand\cexdownloadoptions.yml `so that the CEX start and end dates go from 2007-2009 rather than 1996-2019. Most of the code only needs the shorter time-span. 
