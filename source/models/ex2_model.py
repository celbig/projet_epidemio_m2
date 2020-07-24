import numpy as np
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
from scipy.optimize import root


def ode_system_factory(p, r, pi, mu, nu, nu_v):
    return lambda t, x : np.array([
    	(1-p)*pi - mu *x[0] - (0.5 * (1-np.exp(-2*nu)) * x[2] + 0.5 * (1-np.exp(-2*nu_v)) * x[3])*x[0],
    	p * pi - mu * x[1] - (1-r) * (0.5 * (1-np.exp(-2*nu)) * x[2] + 0.5 * (1-np.exp(-2*nu_v)) * x[3]) * x[1],
    	(0.5 * (1 - np.exp(-2 * nu)) * x[2] + 0.5 * (1 - np.exp(-2 * nu_v)) * x[3]) * x[0] - (mu + nu) * x[2],
    	(1-r) * (0.5 * (1-np.exp(-2 * nu)) * x[2] + 0.5 * (1-np.exp(-2 * nu_v)) * x[3]) * x[1] - (mu + nu_v) * x[3]
        ])
def r0(p, r, pi, mu, nu, nu_v):
    return (((1-p)*0.5 * (1-np.exp(-2*nu))/(mu + nu)) + ( (p * (1-r) *0.5 * (1-np.exp(-2*nu_v))) / (mu + nu_v ))) * pi/mu 


class Ex2_Model():
	"""ODE system to model rabbit in a population of foxes"""
	def __init__(self, 
			S0,
			Sv0, 
			I0,
			Iv0,
			p,
			r,
			pi,
			mu,
			nu,
			nu_v,
			tmax = 50,
			Nt = 1000):
		self.p = p
		self.r = r 
		self.pi = pi 
		self.mu = mu
		self.nu = nu
		self.nu_v = nu_v
		self.S0 = S0
		self.Sv0 = Sv0
		self.I0 = I0
		self.Iv0 = Iv0
		self.tmax = tmax 
		self.Nt = Nt
		
		self.__init()

	def __init(self):
		self.S = None
		self.Sv = None
		self.I = None
		self.Iv = None
		self.t = None
		self.solver = None
		self.f = ode_system_factory(self.p, self.r, self.pi, self.mu, self.nu, self.nu_v)




	def r0(self):
		return r0(self.p, self.r, self.pi, self.mu, self.nu, self.nu_v)

	def set_param(self,
		p = None,
		r = None, 
		pi = None,
		mu = None,
		nu = None,
		nu_v = None,
		S0 = None,
		Sv0 = None,
		I0 = None,
		Iv0 = None,
		tmax = None,
		Nt = None):

		if p is not None:
			self.p = p
		if r is not None:
			self.r = r 
		if pi is not None:
			self.pi = pi 
		if mu is not None:
			self.mu = mu
		if nu is not None:
			self.nu = nu
		if nu_v is not None:
			self.nu_v = nu_v
		if S0 is not None:
			self.S0 = S0
		if Sv0 is not None:
			self.Sv0 = Sv0
		if I0 is not None:
			self.I0 = I0
		if Iv0 is not None:
			self.Iv0 = Iv0
		if tmax is not None:
			self.tmax = tmax 
		if Nt is not None:
			self.Nt = Nt
		
		self.__init()


	def solve(self, tmax = None, Nt = None, graph = True):
		if tmax is not  None :
			self.tmax = tmax 
		if Nt is not None  :
			self.Nt = Nt

		self.solver = solve_ivp(self.f, (0, self.tmax), [self.S0, self.Sv0, self.I0, self.Iv0],
			dense_output=True)
		self.t = np.linspace(0, self.tmax, self.Nt)

		z = self.solver.sol(self.t).T

		self.S = z[:,0]
		self.Sv = z[:,1]
		self.I = z[:,2]
		self.Iv = z[:,3]

		if graph:
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Solution of the system}$ "
            				# \\\\" +
		                	# ("a = {}, b={}, \\alpha = {}, \\beta = {},\\\\ \\gamma ={}, \\sigma = {}, R_0 = {}$"
		                	# .format(self.a, self.b, self.alpha, self.beta, self.gamma,
		                		# self.sigma, self.r0()))
	                	) 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.t,
					y = self.S,
					mode = "lines",
					name = "S" )
				)
			fig.add_trace(
				go.Scatter(
					x = self.t,
					y = self.Sv,
					mode = "lines",
					name = "Sv" )
				)
			fig.add_trace(
				go.Scatter(
					x = self.t,
					y = self.I,
					mode = "lines",
					name = "I" )
				)
			fig.add_trace(
				go.Scatter(
					x = self.t,
					y = self.Iv,
					mode = "lines",
					name = "Iv" )
				)
			fig.update_xaxes(title_text='$t$')
			return fig

	def graph_phase_IS(self):
		if (not self.S is None) and (not self.I is None):
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of I vs S}$") 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.S,
					y = self.I,
					mode = "lines",
					name = "Traj" )
				)
			# fig.add_trace(
			# 	go.Scatter(
			# 		x =[I_tilde],
			# 		y = [S_tilde],
			# 		mode = "markers",
			# 		name = "Equilibrium" )
			# 	)
			fig.update_xaxes(title_text='$S$')
			fig.update_yaxes(title_text='$I$')
			return fig

	def graph_phase_IvS(self):
		if (not self.S is None) and (not self.Iv is None):
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of Iv vs S}$") 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.S,
					y = self.Iv,
					mode = "lines",
					name = "Traj" )
				)
			# fig.add_trace(
			# 	go.Scatter(
			# 		x =[I_tilde],
			# 		y = [S_tilde],
			# 		mode = "markers",
			# 		name = "Equilibrium" )
			# 	)
			fig.update_xaxes(title_text='$S$')
			fig.update_yaxes(title_text='$I_v$')
			return fig

	def graph_phase_IIv(self):
		if (not self.I is None) and (not self.Iv is None):
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of I vs I_v}$") 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.Iv,
					y = self.I,
					mode = "lines",
					name = "Traj" )
				)
			# fig.add_trace(
			# 	go.Scatter(
			# 		x = [I_tilde],
			# 		y = [E_tilde],
			# 		mode = "markers",
			# 		name = "Equilibrium" )
			# 	)
			fig.update_xaxes(title_text='$I_v$')
			fig.update_yaxes(title_text='$I$')
			return fig

	def graph_phase_IvI(self):
		if (not self.Iv is None) and (not self.I is None):
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of I_v vs I}$") 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.I,
					y = self.Iv,
					mode = "lines",
					name = "Traj" )
				)
			# fig.add_trace(
			# 	go.Scatter(
			# 		x = [I_tilde],
			# 		y = [E_tilde],
			# 		mode = "markers",
			# 		name = "Equilibrium" )
			# 	)
			fig.update_xaxes(title_text='$I$')
			fig.update_yaxes(title_text='$I_v$')
			return fig

	def graph_phase_SvS(self):
		if (not self.S is None) and (not self.Sv is None):
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of S_v vs S}$") 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.S,
					y = self.Sv,
					mode = "lines",
					name = "Traj" )
				)
			# fig.add_trace(
			# 	go.Scatter(
			# 		x = [I_tilde],
			# 		y = [E_tilde],
			# 		mode = "markers",
			# 		name = "Equilibrium" )
			# 	)
			fig.update_xaxes(title_text='$S$')
			fig.update_yaxes(title_text='$S_v$')
			return fig

