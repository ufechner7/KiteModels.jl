# Exported Types

```@meta
CurrentModule = KiteModels
```

## Basic types
```@docs
SimFloat
KVec3
SVec3
AbstractKiteModel
AKM
```

## Struct KPS3 and KPS4 and SymbolicAWEModel
```@docs
KPS3
KPS4
SymbolicAWEModel
```
These structs store the state of the one point model and four point model. Only in unit tests
it is allowed to access the members directly, otherwise use the input and output functions.
