# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using KiteModels

set = load_settings("system_ram.yaml")
rak = SymbolicAWEModel(set)
KiteModels.init_sim!(rak; remake=false, reload=false)