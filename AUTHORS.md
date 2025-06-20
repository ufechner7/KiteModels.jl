<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner
SPDX-License-Identifier: MIT
-->

# Project Authors
This project was started by Uwe Fechner in 2022 with the goal to make the results 
of his PhD thesis on the simulation and control of kite power systems, re-coded in a fast and modern programming language, available to the public.

Bart van de Lint joined the project in 2024. His major contribution is the first ModelingToolkit (MTK) based model, a ram-air kite model.

## Developers
- Uwe Fechner, Delft/ Den Haag, The Netherlands -  original author of the `KPS3` and `KPS4` kite models
- Bart van de Lint, Trondheim, Norway - author of the `SymbolicAWEModel` model

## Contributors
- Daan van Wolffelaar, Delft, The Netherlands - contributed the scripts 
  - `test_inertia_calculation.jl` 
  - `calculate_rotational_inertia.jl` and the related functions and made major contributions to 
  - `test_orientation.jl`
  
  He also tested the KPS4 model against real flight data and was able to find bugs and determine parameters to achieve a reasonable match between the model and a real system.
- Friso Broekhuizen, Delft, The Netherlands contributed
  - the code in the branch `kps5`, in particular the `KPS5` kite model
- Roland Schmehl, Delft contributed the 3D CAD file `kite.obj`.

