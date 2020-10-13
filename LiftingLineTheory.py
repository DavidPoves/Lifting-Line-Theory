import numpy as np
import sympy as sp


class LiftingLineTheory(object):
	def __init__(self):
		"""
		Class implementing all the necessary methods obtain wing data using the Prandtl Lifting Line Theory. The process
		to solve the system, the user should specify the stations at which the equation should be evaluated and the
		terms to be considered for the Fourier series involved. The steps to solve are:
		1. Initialize this class. No inputs are required.
		2. Load all the geometric properties of the wing: aspect_ratio, taper_ratio, surface, chord_distribution and
		zero_lift_angle. The chord distribution must be an unique function of theta. The zero lift angle input can be
		a constant or a function of theta. Both of them should be (or return) radians.
		3. Call the solve method from this class to obtain the desired Fourier coefficients. The coefficients to be
		evaluated are the ones introduced by the user in considered_terms. Other inputs are angle_attack, which can be
		a constant or a function of theta (in radians) and the station angles.
		4. Once the coefficients are known, the rest of the parameters can be easily evaluated.
		"""
		# Preallocate data.
		self.AR = None
		self.surface = None
		self.taper_ratio = None
		self.chord_distribution = None
		self.zero_lift_angle = None
		self.b = None  # Wing span.
		self.considered_terms = None  # Considered terms of the Fourier series.
		self.station_angles = None  # Station angles used to define the system to be solved.
		self.system = None  # Array containing th system to be solved.
		self.angle_of_attack = None
		self.symbols = None
		self.Fourier_coefficients = dict()

	def load_geometric_parameters(self, aspect_ratio, taper_ratio, surface, chord_distribution, zero_lift_angle):
		"""
		Load geometric parameters and compute some other useful params.
		:param aspect_ratio: Aspect ratio of the entire wing.
		:param taper_ratio: Ratio between the tip and root chords.
		:param surface: Surface of the whole wing.
		:param chord_distribution: Distribution of the chord along the span direction as a function of theta c(theta),
		where theta is the polar angle in radians.
		:param zero_lift_angle:
		:return:
		"""
		self.AR = aspect_ratio
		self.taper_ratio = taper_ratio
		self.surface = surface
		self.chord_distribution = chord_distribution
		self.zero_lift_angle = zero_lift_angle

		# Compute auxiliary geometric params.
		self.b = np.sqrt(self.AR*self.surface)  # Wing span.

	@staticmethod
	def get_zero_lift_angle(inp, theta):
		"""
		Deal with the zero lift angle of attack definition.
		:param inp: Input of the user. It can be a function defining the zero lift angle of attack as a function of the
		station angle or a constant.
		:param theta: Angle at which the zero lift angle of attack must be evaluated.
		:return: Zero lift angle of attack. Its units will depend on the user input.
		"""
		if hasattr(inp, '__call__'):  # Check if the input is a function.
			return inp(theta)
		else:
			return inp

	@staticmethod
	def get_angle_of_attack(inp, theta):
		"""
		Deal with the angle of attack definition.
		:param inp: Input from the user of the angle of attack. It can be an evaluable function alpha(theta), where
		theta is a station angle, or it can be a constant in radians.
		:param theta: Station angle in radians.
		:return:
		"""

		if hasattr(inp, '__call__'):  # Check if the input is a function.
			return inp(theta)
		else:
			return inp

	def define_system(self):
		"""
		Define the system coming from the fundamental Prandtl lifting line coefficient considering the terms desired
		by the user.
		:return:
		"""

		# Define the Fourier terms symbols.
		terms_dict = dict()
		for i in self.considered_terms:
			terms_dict[f'A{i}'] = sp.Symbol(f'A{i}')

		self.system = np.array([])

		# Define the Lifting Line Equation.
		for theta in self.station_angles:
			# Compute the length of the chord at the given station.
			c_theta = self.chord_distribution(theta)

			# Define the Fourier terms
			# f1 = An*sin(n*theta), where theta is the station angle.
			f1 = sp.Add(*[terms_dict[f'A{str(int(n))}'] * sp.sin(n * theta) for n in np.linspace(
				self.considered_terms[0], self.considered_terms[-1], len(self.considered_terms))])

			# f2 = n*An*(sin(n*theta))/(sin(theta))
			f2 = sp.Add(*[n * terms_dict[f'A{str(int(n))}'] * (sp.sin(n * theta)) / (sp.sin(theta))
			              for n in np.linspace(self.considered_terms[0], self.considered_terms[-1],
			                                   len(self.considered_terms))])

			# Define the Sympy equations.
			self.system = np.append(self.system, sp.Eq((2*self.b)/(sp.pi*c_theta)*f1 +
			                                           LiftingLineTheory.get_zero_lift_angle(self.zero_lift_angle,
			                                                                                 theta) +
			                                           f2, LiftingLineTheory.get_angle_of_attack(self.angle_of_attack,
			                                                                                     theta)))
		self.symbols = tuple(list(terms_dict.values()))
		self.system = tuple(self.system)

	def solve(self, angle_of_attack, station_angles, considered_terms):
		"""
		Solve for the considered Fourier coefficients given an angle of attack and station angles, which come from the
		change of variable y = -b/2 * cos(theta).
		:param angle_of_attack: Angle of attack in radians.
		:param station_angles: List containing the station angles in radians.
		:param considered_terms: List containing the considered terms for the Fourier series. Must match in length with
		the station_angles variable.
		:return:
		"""

		# Do previous check.
		assert len(station_angles) == len(considered_terms), 'Lengths of station_angles and considered_terms variables do not match.'

		self.angle_of_attack = angle_of_attack
		self.station_angles = station_angles
		self.considered_terms = considered_terms

		# Define the Lifting line theory system with the considered terms.
		self.define_system()

		# Solve the system.
		ans = sp.solve(self.system, self.symbols)

		# Extract the coefficients.
		for key, value in ans.items():
			self.Fourier_coefficients[str(key)] = value

	def compute_lift_coefficient(self):
		"""
		Compute the lift coefficient of the whole wing.
		:return: Lift coefficient of the entire wing.
		"""
		return self.Fourier_coefficients[str(sp.Symbol('A1'))] * np.pi * self.AR

	def compute_induced_drag_parameter(self):
		"""
		Compute the induced drag parameter as given in Anderson's book page 449.
		:return: Induced drag parameter.
		"""
		# For this computation, the first term does not enter in the computation as an index of the summation.
		list_terms = self.considered_terms[1:]

		# Compute aux. term.
		term = sp.Add(*[n * (self.Fourier_coefficients[str(sp.Symbol(f'A{str(int(n))}'))] / self.Fourier_coefficients[
			str(sp.Symbol('A1'))]) ** 2 for n in np.linspace(list_terms[0], list_terms[-1], len(list_terms))]).doit()

		return term

	@staticmethod
	def compute_oswald_coefficient(induced_drag_parameter):
		"""
		Compute the Oswald coefficient from the induced drag parameter.
		:param induced_drag_parameter: The induced drag parameter.
		:return: Oswald Coefficient.
		"""
		return (1+induced_drag_parameter)**(-1)

	def compute_induced_drag_coefficient(self, oswald_coefficient=None, induced_drag_parameter=None):
		"""
		Compute the induced drag coefficient given the lift coefficient and the oswald coefficient/ind.drag param.
		:param oswald_coefficient: Oswald coefficient. Optional, default is None.
		:param induced_drag_parameter: Induced drag parameter. Optional, default is None.
		:return: The induced drag coefficient.
		"""
		C_L = self.Fourier_coefficients['A1'] * np.pi * self.AR
		if oswald_coefficient is not None:
			return C_L**2/(np.pi*oswald_coefficient*self.AR)
		elif induced_drag_parameter is not None:
			return C_L**2/(np.pi*self.AR) * (1+induced_drag_parameter)
