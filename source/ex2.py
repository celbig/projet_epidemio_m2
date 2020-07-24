from models import *
import numpy as np
import os

if not os.path.exists("images"):
    os.mkdir("images")



# Simulations 
# Evolution pour r0 < 1 
model = Ex2_Model(200, 200, 200, 200, 
		p = 3/4,
		r = 8/10, 
		pi = 50,
		mu = 1/10,
		nu = 1/2000,
		nu_v = 1/4000,
		tmax = 200)
print("R0 = {}".format(model.r0())) # R0 = 0.8085653742524646
fig1 = model.solve()
# fig3.show()
fig1.write_image("../figures/ex2_r0_lt_1.eps")

# Evolution pour r0 > 1 
model = Ex2_Model(200, 200, 200, 200, 
		p = 1/4,
		r = 8/10, 
		pi = 50,
		mu = 1/10,
		nu = 1/2000,
		nu_v = 1/4000,
		tmax = 200) 
print("R0 = {}".format(model.r0())) # R0 = 1.9270676730509508
fig2 = model.solve()
# fig2.show()
fig2.write_image("../figures/ex2_r0_gt_1.eps")


# I vs S et Iv vs
model = Ex2_Model(200, 200, 200, 200, 
		p = 1/4,
		r = 8/10, 
		pi = 30,
		mu = 1/10,
		nu = 1/2000,
		nu_v = 1/4000,
		tmax = 200) 
print("R0 = {}".format(model.r0())) # R0 = 1.1562406038305706
model.solve()
fig3 = model.graph_phase_IS()
# fig3.show()
fig3.write_image("../figures/ex2_IS_1.eps")
fig4 = model.graph_phase_IvS()
# fig4.show()
fig4.write_image("../figures/ex2_IvS_1.eps")

model.set_param(p = 1/2, nu = 1/500, nu_v = 1/200)
print("R0 = {}".format(model.r0()))  # R0 = 7.261257124158084
model.solve()
fig5 = model.graph_phase_IS()
# fig5.show()
fig5.write_image("../figures/ex2_IS_2.eps")
fig4 = model.graph_phase_IvS()
# fig4.show()
fig4.write_image("../figures/ex2_IvS_2.eps")

model.set_param(p = 8/10, nu = 1/2000, nu_v = 1/5000)
print("R0 = {}".format(model.r0()))  # R0 = 0.6569124713425646
model.solve()
fig6 = model.graph_phase_IS()
# fig6.show()
fig6.write_image("../figures/ex2_IS_3.eps")
fig7 = model.graph_phase_IvS()
# fig7.show()
fig7.write_image("../figures/ex2_IvS_3.eps")


# Iv vs I et Sv vs S
model = Ex2_Model(200, 200, 200, 200, 
		p = 1/4,
		r = 8/10, 
		pi = 30,
		mu = 1/10,
		nu = 1/2000,
		nu_v = 1/4000,
		tmax = 500) 
print("R0 = {}".format(model.r0())) # R0 = 1.1562406038305706
model.solve().show()
fig8 = model.graph_phase_IvI()
# fig8.show()
fig8.write_image("../figures/ex2_IvI_1.eps")
fig9 = model.graph_phase_SvS()
# fig9.show()
fig9.write_image("../figures/ex2_SvS_1.eps")

model.set_param(S0 = 195, Sv0 = 0, I0 = 10, Iv0 = 10, tmax = 500)
model.solve().show()
fig10 = model.graph_phase_IvI()
# fig10.show()
fig10.write_image("../figures/ex2_IvI_2.eps")
fig11 = model.graph_phase_SvS()
# fig11.show()
fig11.write_image("../figures/ex2_SvS_2.eps")

model.set_param(S0 = 1995, Sv0 = 0, I0 = 10, Iv0 = 10, tmax = 500)
model.solve().show()
fig12 = model.graph_phase_IvI()
# fig12.show()
fig12.write_image("../figures/ex2_IvI_3.eps")
fig13 = model.graph_phase_SvS()
# fig13.show()
fig13.write_image("../figures/ex2_SvS_3.eps")
