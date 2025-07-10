# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using KiteModels

set = load_settings("system_ram.yaml")
sam = SymbolicAWEModel(set)
KiteModels.init!(sam; remake=false, reload=false)