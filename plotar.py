import matplotlib.pyplot as plt

def plotar(array):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    ax.set(xlabel='time (it)', ylabel='error', title='Evolução do erro')
    ax.plot(range(len(array)), array)
    ax.grid()
    fig.savefig('atv4-erroAG.png', format='png')
    plt.show()

arr = [3.58072
    ,0.499111
    ,0.499111
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.458579
    ,0.37719
    ,0.37719
    ,0.37719
    ,0.37719
    ,0.37719
    ,0.37719
    ,0.37719
    ,0.37719]

plotar(arr)