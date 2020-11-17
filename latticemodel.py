import numpy as np
import matplotlib.pyplot as plt

"""
FUNKSJONER
"""
# coloumbkraften
def coloumb(q1, q2, r1, r2): 
	R = r1-r2
	Rnorm = np.linalg.norm(R)
	return q1*q2*R/Rnorm**3

# friksjonskraft
def friction(v, k): 
	return -v*k

# kinetisk energi (rotasjonell og translasjonell)	
def KE(w, v, m, I): 
	return 0.5*I*w**2 + 0.5*m*v**2

# potensiell elektrisk energi
def PE(q1, q2, r1, r2): 
	R = r1-r2
	Rnorm = np.linalg.norm(R)
	return 0.5*q1*q2/Rnorm
	
# angulært momentum
def angular_momentum(I, w):
	return I*w

# torque
def torque(f, r): 
	return np.cross(f, r)

# arm
def arm(r, ms): 
	return r-ms
# massesenter	
def com(m1, m2, r1, r2): 
	M = m1+m2
	return (m1*r1 + m2*r2)/M 
	
# Integrator Eulers metode
def integrator(q, a_mom, thet, I1, torq, dt, rc):
	w = a_mom/I1
	a_mom_ny = a_mom+torq*dt
	w_ny = a_mom_ny/I1
	theta_ny = thet + dt*w
	if (q==1): 
		r1_ny = np.array([8/np.sqrt(65)*bl*np.cos(theta_ny), 8/np.sqrt(65)*bl*np.sin(theta_ny)])+rc
	else: 
		r1_ny = np.array([-1/np.sqrt(65)*bl*np.cos(theta_ny), -1/np.sqrt(65)*bl*np.sin(theta_ny)])+rc


	return a_mom_ny, theta_ny, r1_ny
 
"""
Setter opp systemet. Avstanden mellom massesenterene til molekylene er 1 til de fire nærmeste naboene
"""

bl = 0.5 # bindingslengde
mp = 1 # masse positv del av molekylet
mn = 8 # masse negativ ende av molekylet
N = 4 # lengde side
N2 = N**2 
I = (mp*mn)/(mp+mn)*bl**2

com_list = []

for i in range(N): 
	for j in range(N): 
		com_list.append(np.array([i, j]))

plist = []
nlist = []
theta_l = []
for i in range(N2): 
	random_int = np.random.randint(100)
	theta = 2*np.pi/100*random_int
	theta_l.append(theta)
	plist.append(np.array([com_list[i][0] + 8/np.sqrt(65)*bl*np.cos(theta), com_list[i][1] + 8/np.sqrt(65)*bl*np.sin(theta)]))
	nlist.append(np.array([com_list[i][0] - 1/np.sqrt(65)*bl*np.cos(theta), com_list[i][1] -1/np.sqrt(65)*bl*np.sin(theta)]))
print(plist)
fig1 = plt.figure(figsize=(8,8))

for i in range(N2): 
	plt.plot(plist[i][0], plist[i][1], 'or')
	plt.plot(nlist[i][0], nlist[i][1], 'ob')
plt.show()

a_mom = []
for i in range(N2): 
	a_mom.append(0)
v = []
for i in range(N2): 
	v.append(np.array([0,0]))

dt = 0.01
t = []
potential_Energy = []
kinetic_Energy = []
for i in range(1001): 
	t.append((i+1)*dt)
	pot_E = 0
	kin_E = 0
	for j in range(N2): 
		torq = 0 
		for k in range(N2): 
			if (j!=k): 
				# Beregner kraft i alle punktene
				ppforce = coloumb(1, 1, plist[i*N2+j], plist[i*N2+k])+friction(v[i*N2+j][0], 0.1)
				pnforce = coloumb(1, -1, plist[i*N2+j], nlist[i*N2+k])+friction(v[i*N2+j][0], 0.1)
				npforce = coloumb(-1, 1, nlist[i*N2+j], plist[i*N2+k])+friction(v[i*N2+j][1], 0.1)
				nnforce = coloumb(-1, -1, nlist[i*N2+j], nlist[i*N2+k])+friction(v[i*N2+j][1], 0.1)
				pforce = ppforce + pnforce
				nforce = npforce + nnforce
				# Beregner torque på molekylet
				a_lp = arm(plist[i*N2+j], com_list[j])
				a_ln = arm(nlist[i*N2+j], com_list[j])
				torqp = torque(pforce, a_lp)
				torqn = torque(nforce, a_ln)
				torq += torqp + torqn
				# Beregner potensiell energi for molekylet
				pot_Epp = PE(1, 1, plist[i*N2+j], plist[i*N2+k])
				pot_Epn = PE(1, -1, plist[i*N2+j], nlist[i*N2+k])
				pot_Enp = PE(-1, 1, nlist[i*N2+j], plist[i*N2+k])
				pot_Enn = PE(-1, -1, nlist[i*N2+j], nlist[i*N2+k])
				pot_E += pot_Epp+pot_Epn+pot_Enp+pot_Enn
		amp, tnp, rnp = integrator(1, a_mom[i*N+j], theta_l[i*N+j], I, torq, dt, com_list[j])
		amn, tnn, rnn = integrator(-1, a_mom[i*N+j], theta_l[i*N+j], I, torq, dt, com_list[j])
		a_mom.append(amp)
		theta_l.append(tnp)
		plist.append(rnp)
		nlist.append(rnn)
		v.append([amp*rnp/I, amn*rnn/I])
		kin_E += 0.5*amp**2/I
	potential_Energy.append(pot_E) 
	kinetic_Energy.append(kin_E)

fig2 = plt.figure(figsize=(8,8))	
for i in range(N2): 
	plt.plot(plist[1000*N+i][0], plist[10*N+i][1], 'or')
	plt.plot(nlist[1000*N+i][0], nlist[10*N+i][1], 'ob')
plt.show()

plt.figure(figsize=(8,8))
plt.plot(t, potential_Energy, 'r')
plt.plot(t, kinetic_Energy, 'b')
plt.legend(['Potential Energy', 'Kinetic Energy'])
plt.xlabel('Time in seconds')
plt.show()
		
