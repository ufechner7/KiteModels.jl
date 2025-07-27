# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Pkg

function globaldependencies()
    projectpath = Pkg.project().path
    basepath, _ = splitdir(projectpath)
    Pkg.activate()
    globaldependencies = keys(Pkg.project().dependencies)
    Pkg.activate(basepath)
    globaldependencies
end

if !("LiveServer" in globaldependencies())
    println("Installing LiveServer globally!")
    run(`julia -e 'using Pkg; Pkg.add("LiveServer")'`)
end
if !("LocalCoverage" in globaldependencies())
    println("Installing LocalCoverage globally!")
    run(`julia -e 'using Pkg; Pkg.add("LocalCoverage")'`)
end

using LocalCoverage
import LiveServer as LS

html_dir = tempdir()
coverage = generate_coverage("KiteModels"; run_test = true)
html_coverage(coverage; open = false, dir = html_dir)
LS.serve(launch_browser=true, dir=html_dir)
