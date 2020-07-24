from models import *
import numpy as np
import os

if not os.path.exists("images"):
    os.mkdir("images")

model = Ex1_Model(1000, 500, 0, gamma = 1/5000, tmax = 1000)

print("R0 = {}".format(model.r0()))

(S0, E0, I0, N, fig1) = model.get_equilibre_endemique(graph = True)

# # Find Endemic equilibrium
# fig1.show()
fig1.write_image("../figures/nullclines.eps")
print(S0, E0, I0, N)


## Study  Stability for Endemic equilibrium
model.set_param(gamma  = 1/10000) 
print("R0 = {}".format(model.r0())) 
model.solve(tmax = 50, graph = False)

fig2 = model.graph_phase_EI()
# fig2.show()
fig2.write_image("../figures/study_Endemic_stabillity_1.eps")


model.set_param(S0 = 0, E0 = 50, I0= 4, gamma  = 1/1000, sigma = 3/13)
print("R0 = {}".format(model.r0())) 
graph = model.solve(tmax = 200)
# graph.show()
fig3 = model.graph_phase_EI()
# fig3.show()
fig3.write_image("../figures/study_Endemic_stabillity_2.eps")

model.set_param(S0 = 0, E0 = 50, I0= 4, alpha = 20/73, beta = 1/30)
print("R0 = {}".format(model.r0())) 
graph = model.solve(tmax = 50)
# graph.show()
fig4 = model.graph_phase_EI()
# fig4.show()
fig4.write_image("../figures/study_Endemic_stabillity_3.eps")



## Numerical study R0 <1 
model = Ex1_Model(300, 300, 300, gamma = 1/1000, tmax = 50)
print("R0 = {}".format(model.r0())) 
fig5 = model.solve()
# fig5.show()
fig5.write_image("../figures/numerical_study_R0_lt1_1.eps")
fig6 = model.graph_N()
# fig6.show()
fig6.write_image("../figures/numerical_study_R0_lt1_2.eps")
fig7 = model.graph_phase_SE()
# fig7.show()
fig7.write_image("../figures/numerical_study_R0_lt1_3.eps")
fig8 = model.graph_phase_SI()
# fig8.show()
fig8.write_image("../figures/numerical_study_R0_lt1_4.eps")


## Numerical study R0 > 1 
model = Ex1_Model(300, 300, 300, gamma = 1/2500, tmax = 50)
print("R0 = {}".format(model.r0())) 
fig9 = model.solve()
# fig9.show()
fig9.write_image("../figures/numerical_study_R0_gt1_1.eps")
fig10 = model.graph_N()
# fig10.show()
fig10.write_image("../figures/numerical_study_R0_gt1_2.eps")
fig11 = model.graph_phase_SE()
# fig11.show()
fig11.write_image("../figures/numerical_study_R0_gt1_3.eps")
fig12 = model.graph_phase_SI()
# fig12.show()
fig12.write_image("../figures/numerical_study_R0_gt1_4.eps")