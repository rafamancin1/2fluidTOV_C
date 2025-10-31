import twofluidTOV

geom_to_Msun = 0.0006772199944005382

eos1 = twofluidTOV.EOS_Tabular("APR")
#eos2 = twofluidTOV.EOS_Poly("Polytropic", 1e7, 2)
eos2 = twofluidTOV.EOS_SIDM(220, 3.1415)
sols = twofluidTOV.TOV_Family(eos1, eos2)
res = sols.calc_lambda_and_mass_directly(0.46, 0.5)
print(f"Mass={res.M * geom_to_Msun}")
print(f"Lambda={res.lambda_param}")


