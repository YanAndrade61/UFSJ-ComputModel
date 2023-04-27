import numpy as np
from scipy.integrate import  solve_ivp
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import math

path = 'data/'

def odeSystem(t, u, alpha, beta, gama):

    S, I, R = u
    dS_dt = -beta*S*I + alpha*R
    dI_dt =  beta*S*I - gama*I
    dR_dt =  gama*I   - alpha*R
 
    return [dS_dt,dI_dt,dR_dt]  

def isReferenceTime(times, ct):
    for t in times: 
        if (abs(ct - t) <= 10**(-5)):
            return True 
    return False

def solve(x):
    global data, reference_times
    dt = 0.01
    tfinal = 50
    times = np.arange(0,tfinal+dt,dt)

    N = 1000
    S = 0.995*N
    I = 0.005*N
    R = N - S - I
    u = [S, I ,R]

    alpha, beta, gama = x
    params = (alpha,beta,gama)
    
    def solveOde(t, y):
        return odeSystem(t, y, *params)

    results = solve_ivp(solveOde,(0, tfinal), u, t_eval=times, method='Radau')    

    u = results.y[:3,:]
    errorS = errorI = errorR = 0
    sumS = sumI = sumR = 0
    i = j = 0  
    for t in times:
        if isReferenceTime(reference_times,t):
            
            Sdata = data[i][1]
            Idata = data[i][2]
            Rdata = N - Sdata - Idata
            
            errorS += (u[0][j] - Sdata)*(u[0][j] - Sdata) 
            errorI += (u[1][j] - Idata)*(u[1][j] - Idata) 
            errorR += (u[2][j] - Rdata)*(u[2][j] - Rdata)

            sumS += Sdata*Sdata
            sumR += Rdata*Rdata
            sumI += Idata*Idata
            
            i += 1
        j += 1

    errorS = math.sqrt(errorS/sumS)                 
    errorI = math.sqrt(errorI/sumI)                  
    errorR = math.sqrt(errorR/sumR)

    return errorS + errorI + errorR

def cb(x, convergence):
    global progress_err
    progress_err.append(solve(x))

if __name__ == "__main__":
    global data, reference_times, progress_err 
    data = np.loadtxt(path+'sir.csv', delimiter=',')
    reference_times = data[:,0]
    progress_err = []
    bounds = [
        (0.01, 1), (0.01, 1), (0.01,1)
    ]

    #chama evolução diferencial, result contém o melhor individuo
    solucao = differential_evolution(solve, 
                                     bounds, 
                                     strategy='best1bin',
                                    #  updating="deferred", 
                                     maxiter=30,
                                     popsize=100,
                                     atol=10**(-3),
                                     tol=10**(-3),
                                     mutation=0.8,
                                     recombination=0.5,
                                     disp=True, 
                                     workers=5,
                                     callback=cb)
    
    print(solucao.x)
    print(solucao.success)
    #saving the best offspring...
    np.savetxt('solucao_ajuste.txt',solucao.x, fmt='%.2f')        
    best=solucao.x
    error = solve(best)
    print("ERROR ", error)

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    ax.set(xlabel='time (it)', ylabel='error', title='Evolução do erro')
    ax.plot(range(len(progress_err)), progress_err)
    ax.grid()
    fig.savefig('atv4-erro.png', format='png')
    plt.show()