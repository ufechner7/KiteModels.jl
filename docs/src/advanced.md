```@meta
CurrentModule = KiteModels
```
# Advanced usage
For advanced users it is suggested to install git, bash and vscode or vscodium in addition to Julia. vscode and vscodium both have a very good plugin for Julia support, see [https://www.julia-vscode.org](https://www.julia-vscode.org/).
For Windows users, you can find some installation instructions here: [Julia on Windows](https://github.com/ufechner7/KiteViewer/blob/main/doc/Windows.md) .

Whe using vscode, I do NOT use the Julia terminal provided by vscode, but the normal bash terminal. From this terminal I start Julia with ```julia --project``` or a different command as explained below. This makes it easier to understand what happens and is also faster when you need to restart.

For Ubuntu Linux I use the following ppa to install vscode and to keep it up-to-date: [https://www.ubuntuupdates.org/ppa/vscode](https://www.ubuntuupdates.org/ppa/vscode) .

## Creating a custom system image
To reduce the startup time it is suggested to use a custom system image that contains all the packages you use on a daily base in compiled form.

From a bash prompt you can create one using the following commands, assuming you created a folder "test" as explained on the Quickstart page:
```bash
cd test
julia --project
```
```julia
using KiteModels
copy_bin()
```
If you enter shell mode by pressing ";" you should see the following files:
```
shell> tree
.
├── bin
│   ├── create_sys_image
│   └── run_julia
├── data
│   ├── settings.yaml
│   └── system.yaml
├── examples
│   ├── compare_kps3_kps4.jl
│   ├── plot2d.jl
│   └── simulate.jl
├── Manifest.toml
└── Project.toml

2 directories, 7 files
```
Now leave Julia with the command ```exit()``` and then type:
```bash
./create_sys_image --update
```
This will take about 6 min on a  i7-10510U CPU. You should now see a new file in the bin folder:
```
~/repos/test/bin$ ls -lah kps*
-rwxrwxr-x 1 ufechner ufechner 344M apr 18 18:23 kps-image-1.7.so
```
You can launch julia such that it makes use of this system image with the commands:
```
cd ..
./bin/run_julia
```
If you now run any of the examples the time-to-first-plot (TTFP) should be less than 1s:
```
julia> @time include("examples/simulate.jl")
lift, drag  [N]: 690.54, 136.91
Average number of callbacks per time step: 744.8983333333333
  4.566270 seconds (9.83 M allocations: 843.155 MiB, 5.78% gc time, 20.95% compilation time)

julia> 
```
A second run of this command needs about 3.5 s which means the startup time (load and compilation time of the libraries) has been reduced to about 1s.

Without a system image the first time execution of the script "simulate.jl" on the same computer needs about 52 seconds
while the time for the second execution is the same (3.5s). So we save now about 47s of time after each restart.

## Outlook

The next steps are:
- add the possibility to create log files of the simulation
- add winch model
- add controllers

Because the code for this already fully exists in Python and partially in Julia these steps will be in place soon.
