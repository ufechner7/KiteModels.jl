```@meta
CurrentModule = KiteModels
```
# Advanced usage
For advanced users it is suggested to install git, bash and vscode or vscodium in addition to Julia. vscode and vscodium both have a very good plugin for Julia support, see [https://www.julia-vscode.org](https://www.julia-vscode.org/).
Installation instructions: [Julia and VSCode](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) .

Whe using vscode, I do NOT use the Julia terminal provided by vscode, but the normal bash terminal which is also available in vscode by selecting **Terminal->New Terminal** From this terminal I start Julia with ```julia --project``` or a different command as explained below. This makes it easier to understand what happens and is also faster when you need to restart.

## Creating a custom system image
To reduce the startup time it is suggested to use a custom system image that contains all the packages you use on a daily base in compiled form.

For end users, follow the instructions on [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl) .
Package developers can follow these instructions:

1. Go to the website https://github.com/ufechner7/KiteModels.jl and click on the **Fork** button at the top right.
2. clone the new repository which is owned by you with a command similar to this one: ```git clone https://github.com/aenarete/KiteModels.jl``` Your own git user name must appear in the URL, otherwise you will not be able to push your changes.

From a bash prompt you can create one using the following commands, assuming you created a folder "test" as explained on the Quickstart page:
```bash
cd KiteModels.jl
julia --project
```
```julia
using Pkg
Pkg.instantiate()
Pkg.precompile()
```
If you enter shell mode by pressing ";" and type the command ```tree```you should see the following files:
```
shell> tree
.
├── bin
│   ├── create_sys_image
│   ├── create_sys_image2
│   ├── run_julia
│   └── run_julia2
├── data
│   ├── settings.yaml
│   └── system.yaml
├── docs
│   ├── data
│   ├── make.jl
│   ├── Project.toml
│   └── src
│       ├── 4-point-kite.png
│       ├── advanced.md
│       ├── assets
│       │   └── logo.png
│       ├── examples_4p.md
│       ├── examples.md
│       ├── functions.md
│       ├── index.md
│       ├── initial_state_4p.png
│       ├── initial_state.png
│       ├── kite.png
│       ├── kite_power_tools.png
│       ├── kps4_hires.png
│       ├── kps4.png
│       ├── parameters.md
│       ├── quickstart.md
│       └── types.md
├── examples
│   ├── compare_kps3_kps4.jl
│   ├── plot2d.jl
│   ├── reel_out.jl
│   ├── simulate_ii.jl
│   └── simulate_simple.jl
├── LICENSE
├── Manifest-1.7.toml.default
├── Manifest-1.8.toml.default
├── Manifest.toml
├── Project.toml
├── README.md
├── src
│   ├── init.jl
│   ├── KiteModels.jl
│   ├── KPS3.jl
│   └── KPS4.jl
└── test
    ├── bench3.jl
    ├── bench4.jl
    ├── create_sys_image2.jl
    ├── create_sys_image.jl
    ├── plot2d.jl
    ├── plot_initial_state.jl
    ├── plot_kps3.jl
    ├── plot_kps4.jl
    ├── runtests.jl
    ├── test_for_precompile.jl
    ├── test_kps3.jl
    ├── test_kps4.jl
    ├── test_staticarrays.jl
    ├── test_steady_state.jl
    ├── test_sundials.jl
    └── update_packages.jl

9 directories, 55 files
```
Now leave Julia with the command ```exit()``` and then type:
```bash
cd bin
./create_sys_image --update
```
This will take about 6 min on a  Ryzen 7950X CPU. You should now see a new file in the bin folder:
```
~/repos/test/bin$ ls -lah kps*
-rwxrwxr-x 1 ufechner ufechner 723M apr 18 18:23 kps-image-1.10-main.so
```
You can launch julia such that it makes use of this system image with the commands:
```bash
cd ..
./bin/run_julia
```
If you now run any of the examples the time-to-first-plot (TTFP) should be less than 10s:
```julia
julia> @time include("examples/simulate_simple.jl")
lift, drag  [N]: 597.47, 129.29
Average number of callbacks per time step: 83.8866
  9.370429 seconds (17.95 M allocations: 1.359 GiB, 4.26% gc time, 50.76% compilation time: 26% of which was recompilation)

julia> 
```
A second run of this command needs about 3.6 s which means the startup time (load and compilation time of the package and the libraries) has been reduced to about 5.77s.

Without a system image the first time execution of the script "simulate_simple.jl" on the same computer is about 16.4 seconds
while the time for the second execution is the same (3.6s). So now about 7s of time are saved after each restart.

## Hints for Developers
### Coding style

- add the packages `TestEnv` and `Revise` to your global environment, not to any project

- avoid hard-coded numeric values like `9.81` in the code, instead define a global constant `G_EARTH` or read this value from a configuration file

- stick to a line length limit of 120 characters

- try to avoid dot operators unless you have to. 
Bad: `norm1        .~ norm(segment)`
Good: `norm1        ~ norm(segment)`

- if you need to refer to the settings you can use `se()` which will load the settings of the active project. To define the active project use a line like `set = se("system_3l.yaml")` at the beginning of your program.
- use the `\cdot` operator for the dot product for improved readability
- use a space after a comma, e.g. `force_eqs[j, i]`
- enclose operators like `+` and `*` in single spaces, like `0.5 * (s.pos[s.num_C] + s.pos[s.num_D])`;  
  exception: `mass_tether_particle[i-1]`
- try to align the equation signs for improved readability like this:
```julia
    tether_rhs        = [force_eqs[j, i].rhs for j in 1:3]
    kite_rhs          = [force_eqs[j, i+3].rhs for j in 1:3]
    f_xy              = dot(tether_rhs, e_z) * e_z
```

## Outlook

The next steps:
- re-implement the KPS4 model using ModelingToolkit
- add a Matlab/ Simulink wrapper similar to the Python wrapper [pykitemodels](https://github.com/ufechner7/pykitemodels)

