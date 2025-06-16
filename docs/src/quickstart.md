```@meta
CurrentModule = KiteModels
```
# Quickstart

## Installation of Julia
For a quick test of this program, it is NOT needed to install VSCode, git or bash. Just installing Julia is sufficient, and that can be done in a few minutes. 

```@raw html
<details>
  <summary>Windows</summary>
```
    
### Windows
Please download and install Julia using `juliaup`. Launch the `Command Prompt` app and type:

```
winget install julia -s msstore
juliaup add 1.10
juliaup update
```
If that doesn't work, download [https://install.julialang.org/Julia.appinstaller](https://install.julialang.org/Julia.appinstaller) and double-click on the downloaded file to install it.

#### Optional
It is suggested to install [Windows Terminal](https://learn.microsoft.com/en-us/windows/terminal/install) . Copy and paste works better, unicode works much better and you can use it with `bash` or `Command Prompt`, whatever you prefer. It is suggested to set one of these two as default using the `Settings` menu of Windows Terminal.
  
```@raw html
</details>
```

```@raw html
<details>
  <summary>Linux</summary>
```

### Linux

Copy and past the following line to install julia:
```
curl -fsSL https://install.julialang.org | sh
```
Restart your terminal, and then execute:
```
juliaup add 1.10
juliaup update
```

It is suggested to add the following line to your ```.bashrc``` file:
```
alias jl='./bin/run_julia'
```
This makes it possible to run Julia with the shortcut `jl` later.

```@raw html
</details>
```

```@raw html
<details>
  <summary>Mac</summary>
```

### Mac
Please download and install `juliaup` as explained at https://github.com/JuliaLang/juliaup .

Restart your terminal, and then execute:
```
juliaup add 1.10
juliaup update
```

```@raw html
</details>
```

## Create a test project
Launch a command prompt and create a folder with the name "test":
```bash
mkdir test
cd test
julia --project="."
```
With the last command, we told Julia to create a new project in the current directory.

You can copy the examples to your project with:
```julia
using KiteModels
KiteModels.install_examples()
```

Your folder structure should now look like this:
```
shell> tree -d
├── data
├── examples
└── test
    ├── data
    └── examples

5 directories
```
You can access the operating system command line by typing the character ";", you then get a "shell" prompt and can enter operating system commands. To leave shell mode, type <BACKSLASH>. 
On windows you need to type ```tree /f``` instead of ```tree``` to see the files.

## Executing the first example
From the Julia prompt you can use the command "include" to execute a script:
```
include("examples/simulate_simple.jl")
```
On Windows you need to type "\\\\" instead of "/":
```
include("examples\\simulate_simple.jl")
```
You will see the 4-point kite fly for 30s. If you want to change the settings of the simulation, open the file `simulate_simple.jl` in your favorite text editor, modify the settings at the beginning of the file and execute the include command again.
You can use the <TAB> key for autocompletion, for example `include("ex<TAB>` completes to `include("examples\` which can save a lot of typing. If you type <TAB> again you get a list of files to choose from.

Try out changing the following default settings:
```
dt = 0.05
STEPS = 600
PLOT = true
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
```
Now you can quit Julia with the command ``` exit()```. If you want to launch Julia again, be sure to be in the correct folder and then type ```julia --project```. Without the parameter ```--project``` it will not load your project settings.

The first run of the script will be slow because Julia must compile the code. The second and any further run is very fast, but only as long as you do not leave your Julia session.

## Comparing the one-point and the four-point kite model
Start Julia in the project folder you created before:
```bash
cd test
julia --project
```
and then execute the command
```
using KiteModels
include("examples/compare_kps3_kps4.jl")
```
Use the command ```include("examples\\compare_kps3_kps4.jl")``` on Windows.

The last view of the animation should look like this:

![Initial State](kite.png)

You can save what you see with the command ```savefig("kite.png")```.

### Exercise
Modify the variable ```ALPHA_ZERO``` in line 11 of the script until the lift force of the 1 point model and the 4 point model match.

## Questions?
If you have any questions, please ask in the Julia Discourse forum in the section [modelling and simulation](https://discourse.julialang.org/c/domain/models) , or in in the section [First steps](https://discourse.julialang.org/c/first-steps) . The Julia community is very friendly and responsive.