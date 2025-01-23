# particleSimulation
Simulate particle for Chemistry peoject

## Minimal preliminary process:

#### 1. Download Source Code 

- Goto [document git](https://github.com/MonkieDein/SingleLayerDiffusionModeling) > Code > Download ZIP and unzip folder.

#### 2. Download julia

- Install [Julia](https://julialang.org/downloads/) > Current stable release > select your OS installer

- Click the julia installer and install julia.

- When the installation completed select 
- [x] Run julia > Finish

#### 3. Run julia code

- Open the code folder (MonkieDein/SingleLayerDiffusionModeling) with vsCode you intended to run.

- Open terminal and run the following command:
    ```julia requirement.jl```
    ```julia code/experiment.jl```

## Optional preliminary process:

#### 1. Github Desktop : Automatic pull Source Code 

- Installation : [Github desktop](https://desktop.github.com/) > Download and run the installer

- Sign in : File > Options > Accounts > Sign in > Continue with browser

- Download code : File > Clone repository > URL >
URL : "https://github.com/MonkieDein/SingleLayerDiffusionModeling"
Local path : use default or choose desire %LOCAL_PATH%

#### 2. Add path to environment variable for window
[Available instruction in youtube](https://www.youtube.com/watch?v=42OXIbdc7bQ)

- open Window Start > type"julia" > right click the app icon > open file location > right click the julia icon in the folder > open file location 

- Now copy your julia file location : default as 
(C:\Users\%USERNAME%\AppData\Local\Programs\Julia-%version%\bin)

- open Window Start > type"Edit your environment variables" and select > Environment Variables... > at User variables for %USERNAME% find the Varible "*Path*" and click it > Edit... > New > paste the julia file location

#### 3. VS code : For developing and running code 

- [Available instruction in youtube](https://www.youtube.com/watch?v=oi5dZxPGNlk&t=204s)


## Directory Structure
- code : 
    - 3Dbasics: Basic physics code for 3d objects dynamics
    - reactChem: Chemistry code, chemical reaction simulation and dynamics
    - reactStatsPlots: Create plots or animation from the simulation
    - experiment: Running a MC simulation for each tg and wp
- plot : 
    - resulting plot from the simulation