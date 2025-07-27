# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# List all methods of KiteModels
using KiteModels

exported_names = names(KiteModels)
exported_functions = filter(n -> isdefined(KiteModels, n) && isa(getfield(KiteModels, n), Function), exported_names)
println("Exported methods:\n")
total = 0
for fun in exported_functions
    global total
    f = getfield(KiteModels, fun)
    mes = methods(f)
    for me in mes
        println(split(repr(me),'@')[1])
        total += 1
    end
end
println("\nTotal: $total")

function check_exported_docs(mod::Module)
    exported_symbols = names(mod, all=false)
    doc_status = Dict{Symbol,Bool}()
    for sym in exported_symbols
        #if isa(getfield(mod, sym), Function)
            doc_status[sym] = Base.Docs.hasdoc(mod, sym)
        #end
    end
    return doc_status
end

# Usage example:
results = check_exported_docs(KiteModels)
undocumented = filter(kv -> !kv[2], results)
println("\nUndocumented exported symbols: \n", keys(undocumented))
