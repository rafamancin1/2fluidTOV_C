import twofluidTOV

geom_to_Msun = 0.0006772199944005382

eos1 = twofluidTOV.EOS_Tabular("APR")
eos2 = twofluidTOV.EOS_SIDM(400, 3.1415)
sols = twofluidTOV.TOV_Family(eos1, eos2)
e1 = 0.6  # baryonic central energy density (GeV/fm^3)
F_chi = 0.1
res = sols.calc_lambda_and_mass_directly(e1, F_chi)
print(f"m_chi=400 MeV: Mass={res.M * geom_to_Msun} Lambda={res.lambda_param} F_chi={res.F_chi}")

# Swap the dark matter EOS without rebuilding the family
sols.set_analytic_eos(twofluidTOV.EOS_SIDM(220, 3.1415))
res = sols.calc_lambda_and_mass_directly(e1, F_chi)
print(f"m_chi=220 MeV: Mass={res.M * geom_to_Msun} Lambda={res.lambda_param} F_chi={res.F_chi}")
