# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using KiteModels

kps4_3L::SymbolicAWEModel = SymbolicAWEModel(KCU(se()))

integrator = KiteModels.init!(kps4_3L; stiffness_factor=0.035, prn=false)
nothing