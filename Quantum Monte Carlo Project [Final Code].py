import tqdm as tqdm 
import numpy as np
import matplotlib.pyplot as plt
import playsound as ps
import json

def S(x,L,T=1,m=1,k=1,type=1): # Define the function S(x) for a fixed T
    if L>1:
        if type==1:
            mat = np.asmatrix(np.identity(L)*(-m*L*T))
            mat = np.roll(mat,1,axis=1)+np.roll(mat,-1,axis=1)
            for i in range(L):
                mat[i,i]=2*m*L*T+k/(L*T)
            xs = np.asmatrix(x)
            out = 1/2 * xs*mat*np.transpose(xs)
        elif type==2:
            out=0
            for i in range(L):
                out += 1/2 * (m*(L*T)**2)*(x[(i+1)%L]**2 - 2*x[i]*x[(i+1)%L]) + 1/2 * (k+m*(L*T)**2)*x[i]**2
            out = out/(L*T)
    else:
        out = (1/(2*T))*x**2
    return out

# Defining a blocking function to perform 'blocking' on the provided dataset
# to account for dependency in consecutive measurements. Determines the block
# size of the dataset split into 10, then finds the block average of the 2nd to
# 10th blocks. Then the std dev is calculated with the block averages and the
# error is given by the standard deviation divided by the square root of size
# of the (dataset - 1).

def blocking(x):
    blocks = 10
    block = int(np.size(x,0)/blocks)
    # xavg = [np.mean(i) for i in x]
    xavg = list(np.average(x,axis=1))
    blockavg = np.array([np.mean(xavg[int(block*n):int(block*(n+1)-1)]) for n in range(1,10)])
    dev = np.std(blockavg)
    err = dev/np.sqrt(blocks-2)
    mean = np.mean(blockavg)
    return mean,err

def DeltaS(x1,x2,j,L,T=1,m=1,k=1):
    old = 1/2 * (m*L*T*(x1[(j+1)%L]**2 + 2*x1[j]**2 + x1[(j-1)%L]**2 - 2*(x1[j]*x1[(j+1)%L]+x1[j]*x1[(j-1)%L])) + k/(L*T) * x1[j]**2)
    new = 1/2 * (m*L*T*(x2[(j+1)%L]**2 + 2*x2[j]**2 + x2[(j-1)%L]**2 - 2*(x2[j]*x2[(j+1)%L]+x2[j]*x2[(j-1)%L])) + k/(L*T) * x2[j]**2)
    return new-old

def savedata(data,err,L,num_T,paths):
    with open('data'+str(L)+'-'+str(num_T)+'-'+str(paths)+'.json','w') as file:
        json.dump(data,file)
        file.close()
    with open('err'+str(L)+'-'+str(num_T)+'-'+str(paths)+'.json','w') as errfile:
        json.dump(err,errfile)
        errfile.close()

def readdata(L,num_T,paths):
    file = open('data'+str(L)+'-'+str(num_T)+'-'+str(paths)+'.json')
    data = json.load(file)
    file.close()
    errfile = open('err'+str(L)+'-'+str(num_T)+'-'+str(paths)+'.json')
    err = json.load(errfile)
    errfile.close()
    return data,err

L=3;k=1;T=1;onesig=0;twosig=0
x2avg=[];x2=[];x2errdata=[]
x2avgdata=[];x2errdata2=[]
paths = 1000000
iterations = 1
typ = 2
num_T = 30
T_range = np.linspace(2/num_T,2,num=num_T)
write = True
if write == True:
    for T in T_range:
        x2avg=[];x2=[];x2errdata=[]
        for loop in range(iterations):
            x_initial=np.zeros(L) # Setting initial parameters for an arbitrary initial state x=0 and T=1
            x_final=np.zeros(L)
            record=[list(x_initial)] # Recording all the mean values of x

            for i in tqdm.tqdm(range(paths)):
                for j in range(L):
                    # Generating a random shift, applying it to the position x_j and calculating
                    # the new value of S(x) and the difference compared to the previous value.

                    x_j = np.random.randint(L)
                    shift=(np.random.rand()-.5)
                    x_final[x_j] += shift
                    Delta_S = DeltaS(x_initial,x_final,x_j,L,T=T)

                    # Working out the probability that the shift is accepted.

                    if Delta_S <= 0:
                        acceptProb = 1
                    else:
                        acceptProb = np.exp(-Delta_S)

                    # If the shift is accepted then the updated value of x is appended
                    # to the list of values and the initial values of both S and x are
                    # updated for the next iteration. If it is rejected then the unaltered
                    # value of x is repeated and neither the x nor S values change.

                    acc = np.random.rand()
                    if  acc <= acceptProb:
                        # record.append(list(x_final))
                        x_initial[x_j] += shift
                    else:
                        # record.append(record[-1])
                        x_final[x_j] -= shift

                record.append(list(x_initial))

            # Plotting the recorded values of x against iterations in order
            # to graphically show the chain. Setting ylim as a function of
            # the recorded values ensures the data is always in a reasonable
            # frame.

            xmean,errmargin = blocking(record)
            x2 = list(np.square(record))
            x2mean,x2err = blocking(x2)
            x2avg.append(x2mean)
            x2errdata.append(x2err)

            print('Expectation of x:',xmean,'+-',errmargin,'\nExpectation of x^2:',x2mean,'+-',x2err)

        x2avgdata.append(x2mean)
        x2errdata2.append(x2err)

    savedata(x2avgdata,x2errdata2,L,num_T,paths)
else:
    x2avgdata,x2errdata2=readdata(L,num_T,paths)

plt.errorbar(T_range,x2avgdata,yerr=x2errdata2,fmt='.')
# plt.plot(T_range,T_range)
plt.plot(np.linspace(2/200,2,200),[1/2 * (np.exp(-1/i)+1)/(1-np.exp(-1/i)) for i in np.linspace(2/200,2,200)])
# plt.plot(np.linspace(2/200,2,200),[2*i*(8*i**2+1)/(1+16*i**2) for i in np.linspace(2/200,2,200)])
plt.xlabel('Temperature (K)')
plt.ylabel('Expectation of x^2')


plt.show()
