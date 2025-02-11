import twofluidTOV

eos1 = twofluidTOV.EOS_Tabular("APR")
eos2 = twofluidTOV.EOS_Poly("Polytropic", 1e7, 2)
'''
sols = twofluidTOV.TOV_Family(eos1, eos2)
K = 1e7
F_chi = 0.1
mass = 1.3
print(sols.calc_lambda_parallel(K, F_chi, mass))
'''
e1 = 0.6
e2 = 0.8

integrator = twofluidTOV.TwoFluid_TOV(eos1, eos2)
result = integrator.integrate_two_fluid_tov(e1, e2)
print(result.F_chi)
print(result.lambda_param)