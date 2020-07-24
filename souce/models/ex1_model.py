import numpy as np
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
from scipy.optimize import root


def ode_system_factory(a, b, alpha, beta, gamma, sigma):
    return lambda t, x : np.array([
    	a * np.sum(x) -b*x[0] - gamma * np.sum(x)*x[0] - beta * x[2] * x[0],
        beta * x[2] * x[0] - b * x[1] - gamma * np.sum(x) * x[1] - sigma * x[1],
        sigma * x[1] - b * x[2] - gamma * np.sum(x) * x[2] - alpha * x[2]
        ])
def r0(a, b, alpha, beta, gamma, sigma):
    return (a-b)*beta * sigma / (gamma *(a + sigma)*(a+alpha))


class Ex1_Model():
	"""ODE system to model rabbit in a population of foxes"""
	def __init__(self, S0, E0, I0,
			a = 1,
			b = 0.5,
			alpha = 1/73,
			beta = 1/79.69,
			sigma = 1/13,
			gamma = 1/10, 
			tmax = 50,
			Nt = 1000):
		self.a = a
		self.b = b 
		self.alpha = alpha 
		self.beta = beta
		self.sigma = sigma
		self.gamma = gamma
		self.S0 = S0
		self.E0 = E0
		self.I0 = I0
		self.tmax = tmax 
		self.Nt = Nt
		
		self.__init()

	def __init(self):
		self.S = None
		self.E = None
		self.I = None
		self.N = None
		self.t = None
		self.solver = None
		self.K = (self.a-self.b)/self.gamma 
		self.f = ode_system_factory(self.a, self.b, self.alpha, self.beta, self.gamma, self.sigma)
		self.iso_1 = lambda x:  1 + (self.alpha * x - self.alpha - self.a -self.sigma) * x / self.sigma
		self.iso_2 = lambda x: self.a / (self.beta * x * (self.a-self.b-self.alpha * x)/self.gamma + self.a - self.alpha * x)



	def r0(self):
		return r0(self.a, self.b, self.alpha, self.beta, self.gamma, self.sigma)

	def set_param(self,
		a = None,
		b = None, 
		alpha = None,
		beta = None,
		sigma = None,
		gamma = None,
		S0 = None,
		E0 = None,
		I0 = None,
		tmax = None,
		Nt = None):

		if a is not None:
			self.a = a
		if b is not None:
			self.b = b 
		if alpha is not None:
			self.alpha = alpha 
		if beta is not None:
			self.beta = beta
		if sigma is not None:
			self.sigma = sigma
		if gamma is not None:
			self.gamma = gamma
		if S0 is not None:
			self.S0 = S0
		if E0 is not None:
			self.E0 = E0
		if I0 is not None:
			self.I0 = I0
		if tmax is not None:
			self.tmax = tmax 
		if Nt is not None:
			self.Nt = Nt
		
		self.__init()


	def __convert_equilibre_endemique(self,S,I):
	    N = (self.a -self.b -self.alpha *I ) / self.gamma
	    E = N * (1 - S -I)
	    I = N * I
	    S = N*S
	    return (E, I, S, N)

	def solve(self, tmax = None, Nt = None, graph = True):
		if tmax is not  None :
			self.tmax = tmax 
		if Nt is not None  :
			self.Nt = Nt

		self.solver = solve_ivp(self.f, (0, self.tmax), [self.S0, self.E0, self.I0],
			dense_output=True)
		self.t = np.linspace(0, self.tmax, self.Nt)

		z = self.solver.sol(self.t).T

		self.N = np.sum(z, axis = 1)
		self.S = z[:,0]
		self.E = z[:,1]
		self.I = z[:,2]

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
					y = self.E,
					mode = "lines",
					name = "E" )
				)
			fig.add_trace(
				go.Scatter(
					x = self.t,
					y = self.I,
					mode = "lines",
					name = "I" )
				)
			fig.add_shape(
			        # Line Horizontal
			            type="line",
			            x0=0,
			            y0=self.K,
			            x1=self.tmax,
			            y1=self.K,
			            line=dict(
			                # color="LightSeaGreen",
			                # width=4,
			                dash="dashdot",
			            ),
			    )
			fig.update_xaxes(title_text='$t$')
			return fig

	def graph_N(self):
		if self.N is not None :
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Graph of N}$"
            			# \\\\}"+              	
	                	# "a = {}, b={}, \\alpha = {}, \\beta = {},\\\\ \\gamma ={}, \\sigma = {}, R_0 = {}$"
	                	# .format(self.a, self.b, self.alpha, self.beta, self.gamma,
		                		# self.sigma, self.r0())
	                	) 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.t,
					y = self.N,
					mode = "lines",
					name = "S" )
				)
			fig.update_xaxes(title_text='$t$')
			fig.update_yaxes(title_text='$N$')
			return fig

	def graph_phase_SE(self):
		if (not self.S is None) and (not self.E is None):
			(S_tilde, E_tilde, I_tilde, N) = self.get_equilibre_endemique()
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of S vs E}$"
            			# for}\\\\" +
	                	# "a = {}, b={}, \\alpha = {}, \\beta = {},\\\\ \\gamma ={}, \\sigma = {}, R_0 = {}$""".format(self.a, self.b, self.alpha, self.beta, self.gamma,
		                		# self.sigma, self.r0())
	                	) 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.E,
					y = self.S,
					mode = "lines",
					name = "Traj" )
				)
			fig.add_trace(
				go.Scatter(
					x =[E_tilde],
					y = [S_tilde],
					mode = "markers",
					name = "Equilibrium" )
				)
			fig.update_xaxes(title_text='$E$')
			fig.update_yaxes(title_text='$S$')
			return fig

	def graph_phase_SI(self):
		if (not self.S is None) and (not self.I is None):
			(S_tilde, E_tilde, I_tilde, N) = self.get_equilibre_endemique()
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of S vs I}$"
            			 # for} \\\\" + 
	                	# "a = {}, b={}, \\alpha = {}, \\beta = {},\\\\ \\gamma ={}, \\sigma = {}, R_0 = {}$""".format(self.a, self.b, self.alpha, self.beta, self.gamma,
		                		# self.sigma, self.r0())
	                	) 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.I,
					y = self.S,
					mode = "lines",
					name = "Traj" )
				)
			fig.add_trace(
				go.Scatter(
					x =[I_tilde],
					y = [S_tilde],
					mode = "markers",
					name = "Equilibrium" )
				)
			fig.update_xaxes(title_text='$I$')
			fig.update_yaxes(title_text='$S$')
			return fig

	def graph_phase_EI(self):
		if (not self.S is None) and (not self.I is None):
			(S_tilde, E_tilde, I_tilde, N) = self.get_equilibre_endemique()
			fig = go.Figure(data = [],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Phase graph of E vs I}$"
            			 # for} \\\\" + 
	                	# "a = {}, b={}, \\alpha = {}, \\beta = {},\\\\ \\gamma ={}, \\sigma = {}, R_0 = {}$""".format(self.a, self.b, self.alpha, self.beta, self.gamma,
		                		# self.sigma, self.r0())
	                	) 
            		)
                )

			fig.add_trace(
				go.Scatter(
					x = self.I,
					y = self.E,
					mode = "lines",
					name = "Traj" )
				)
			fig.add_trace(
				go.Scatter(
					x = [I_tilde],
					y = [E_tilde],
					mode = "markers",
					name = "Equilibrium" )
				)
			fig.update_xaxes(title_text='$I$')
			fig.update_yaxes(title_text='$E$')
			return fig


	def get_equilibre_endemique(self, graph = False):
		iso = lambda x: self.iso_1(x) - self.iso_2(x)
		I_tilde = root(iso, 0.5)
		if not I_tilde.success:
			return None
		I_tilde = I_tilde.x[0]
		if I_tilde == 0 or I_tilde >= 1 :
			return None
		S_tilde = self.iso_1(I_tilde)
		(E0, I0, S0, N) = self.__convert_equilibre_endemique(S_tilde,I_tilde)
		if not graph:
			return (S0, E0, I0,  N)

		x = np.linspace(0,1,500)
		graph = go.Figure(data =[],
                layout = go.Layout(
            		title = go.layout.Title(text="$\\text{Nullcines}$"
            			 # for } \\\\" + 
	                	# "a = {}, b={}, \\alpha = {}, \\beta = {},\\\\ \\gamma ={}, \\sigma = {}, R_0 = {}$""".format(self.a, self.b, self.alpha, self.beta, self.gamma,
		                		# self.sigma, self.r0())
	                	) 
            		)
                )
		graph.add_trace(go.Scatter(x=x, y=self.iso_1(x),
                    mode='lines',
                    name='$1 -\\frac{(a + \\sigma + \\alpha - \\alpha x) x}{\\sigma}$'))
		graph.add_trace(go.Scatter(x=x, y=self.iso_2(x),
                    mode='lines',
                    name='$1 -y = \\frac{a}{a - \\alpha x + \\beta N x}$'))
		graph.add_trace(go.Scatter(x=[I_tilde], y=[S_tilde],
                    mode='markers', name='Endemic equilibrium'))
		graph.update_xaxes(title_text='$\\bar I$', range=[0, 1])
		graph.update_yaxes(title_text='$\\bar S$', range=[0, 1])
		print(S0, E0, I0)
		return (S0, E0, I0, N, graph)



