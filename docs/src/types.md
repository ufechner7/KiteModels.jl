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

## Struct KPS3 and KPS4 and KPS4_3L
```@docs
KPS3
KPS4
KPS4_3L
```
These structs store the state of the one point model and four point model. Only in unit tests
it is allowed to access the members directly, otherwise use the input and output functions.
The model `KPS4_3L` is not yet implemented, it is just a copy of `KPS4`.

