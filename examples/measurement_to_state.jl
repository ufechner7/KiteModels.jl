using Revise, KiteModels, ModelingToolkit

s = KPS4_3L(KCU(se("system_3l.yaml")))
model!(s)
