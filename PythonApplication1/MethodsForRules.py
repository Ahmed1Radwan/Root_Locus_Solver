import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
import tkinter as tkn
import math

poles = np.array([0,-25,-50+10j,-50-10j])
zeros = np.array([])
s = sym.Symbol('s')
k = sym.Symbol('k')

# step 2
def calculateNumberOfAsymptotes(numberOfPoles,numberOfZeros):
    print("Number of Asymptotes")
    return numberOfPoles - numberOfZeros

# Step 3
def AnglesOfAsymptotes(numberOfAsymptotes):
    Angles = list()
    if(numberOfAsymptotes <= 0):
        Angles.append(0)

    for x in range(numberOfAsymptotes):
        theta = ((2*x+1)/numberOfAsymptotes)*180
        Angles.append(theta)

    return Angles

# step 4
def calculateCentroidOfAsymptotes(poles,zeros):
    realPartOfPoles = poles.real
    realPartOfZeros = zeros.real
    return (sum(realPartOfPoles) - sum(realPartOfZeros)) / (poles.size-zeros.size)

# step 5 
def calculateBreakPoints(poles,zeros):
    charEquation = charachterEquation(poles,zeros)
    #print(charEquation)
    coff = charEquation.coeffs()
    coff.append(0)
    #print(coff)
    diffOfEquation = np.polyder(coff,1)
    diffOfEquation = diffOfEquation *-1
    BreakPoints = np.roots(diffOfEquation)
    return BreakPoints

# step 6
def calculateRouthTable(polynomial,epsilon,symbolic):
    dim=len(polynomial)
    coeff=dim
    columns=int(np.ceil(coeff/2))
    RA=sym.zeros(coeff,columns) # create the Routh table
    s=sym.Symbol('S')

    #assemble the first and second rows
    top=polynomial[0::2]
    down=polynomial[1::2]
    if len(top)!=len(down):
        down.append(0)

    for i in range(len(top)):
        RA[0,i]+=top[i]
        RA[1,i]+=down[i]

    rows=coeff-2 #number of rows that need determinants
    index=np.zeros(rows)

    for i in range(rows):
        index[rows-i-1]=np.ceil((i+1)/2)
        
    for i in range(2,coeff):  #go from the 3rd row to the last
        if np.sum(np.abs(RA[i-1,:]))==0: # while row is zero
            print('\n Row of zeros detected in row %d. \n Finding Auxillary Polynomial by differentiating row %d'%(i,i-1))
            order=coeff-i+1
            order_arr=np.arange(order,-1,-2)
            poly=0
            #get the polynomial for differentiation
            for k in range(len(order_arr)):
                poly+=(s**(order_arr[k]))*RA[i-2,k]
                diff=sp.diff(poly,s)
                a=sp.Poly(diff,s)
                c=a.coeffs()
            for l in range(columns-len(c)):
                c.append(0)
            for l in range(columns):
                RA[i-1,l]=c[l]
        
        elif RA[i-1:0]==0: # first element is zero
            print('\n First element is zero. Replacing with epsilon')
            RA[i-1:0]=epsilon
            
        for j in range(int(index[i-2])):
            if symbolic:
                RA[i,j]=-sym.det(sym.Matrix([[RA[i-2,0],RA[i-2,j+1]],[RA[i-1,0],RA[i-1,j+1]]]))/RA[i-1,0]
            else:
                RA[i,j]=-my_det(np.array([[RA[i-2,0],RA[i-2,j+1]],[RA[i-1,0],RA[i-1,j+1]]]).astype(float))/RA[i-1,0]
           
            
    return RA

def my_det(X):
    X = np.array(X, dtype='float64', copy=True)
    n = len(X)
    s = 0
    if n != len(X[0]):
      return ValueError
    for i in range(0, n):
      maxElement = abs(X[i, i])
      maxRow = i
      for k in range(i + 1, n):
          if abs(X[k, i]) > maxElement:
              maxElement = abs(X[k, i])
              maxRow = k
      if maxRow != i:
          s += 1
      for k in range(i, n):
          X[i, k], X[maxRow, k] = X[maxRow, k], X[i, k]
      for k in range(i + 1, n):
          c = -X[k, i] / X[i, i]
          for j in range(i, n):
              if i == j:
                  X[k, j] = 0
              else:
                  X[k, j] += c * X[i, j]
    det = (-1)**s
    for i in range(n):
      det *= X[i, i]
    return det


# form the charachter equation from poles and zeros
def charachterEquation(poles,zeros):
    charEquation = 1
    for x in poles:
        charEquation = charEquation *(s-x)

    charEquation = sym.Poly(charEquation,s)
    return charEquation



def Gui():


    root = tkn.Tk()
    root.title('GUI')
    root.geometry('700x700')

    polesLabel = tkn.Label(text = 'poles',font=30).pack()
    polexInfo = tkn.Text(root,width = 60,height=4)
    polexInfo.insert(tkn.INSERT,poles)
    polexInfo.pack()

    zerosLabel = tkn.Label(text = 'Zeros',font = 30).pack()
    zerosInfo = tkn.Text(root,width = 60,height = 4)
    zerosInfo.insert(tkn.INSERT,zeros)
    zerosInfo.pack()

    centroid_Asymptotes = calculateCentroidOfAsymptotes(poles,zeros)
    centroidLabel = tkn.Label(text = 'centroid of asymptotes',font = 20 ).pack()
    centroidInfo = tkn.Text(root,width = 60,height = 4)
    centroidInfo.insert(tkn.INSERT,centroid_Asymptotes)
    centroidInfo.pack()

    RA = calculateRouthTable([1,125,5100,65000,k],0.01,True)
    routhLabel = tkn.Label(text = 'Routh Table',font=20).pack()
    routhInfo = tkn.Text(root,width = 79,height = 4)
    routhInfo.insert(tkn.INSERT,RA)
    routhInfo.pack()

    sRow = RA[3,0]
    KMaxValue = sym.solve(sRow,k)
    s2Function = RA[2,0]*s**2+ KMaxValue[0]
    ImaginaryIntersectionPoints = sym.solve(s2Function,s)
    imagIntersectionLabel = tkn.Label(text='intersection points with Imaginary axis',font = 10).pack()
    pointsInfo = tkn.Text(root,width = 60,height = 4)
    pointsInfo.insert(tkn.INSERT,ImaginaryIntersectionPoints)
    pointsInfo.pack()

    angles_Asymptotes = AnglesOfAsymptotes(len(poles)-len(zeros))
    anglesLabel = tkn.Label(text = 'Angles of asymptotes',font = 20).pack()
    anglesInfo = tkn.Text(root,width=60,height = 4)
    anglesInfo.insert(tkn.INSERT,angles_Asymptotes)
    anglesInfo.pack()

    BreakPoints = calculateBreakPoints(poles,zeros)
    breakPointsLabel = tkn.Label(text='Break away Points',font=20).pack()
    bpInfo=tkn.Text(root,width=60,height=4)
    bpInfo.insert(tkn.INSERT,BreakPoints)
    bpInfo.pack()

    Graph()

    root.mainloop()



def Graph():

    plt.ylim(-150,150)
    plt.xlim(-200,100)
    plt.axvline(x=0,linestyle=':',color='black',linewidth=0.4)
    plt.axhline(y=0,linestyle=':',color='black',linewidth=0.4)

    plt.title('root locus')
    plt.xlabel('real part')
    plt.ylabel('imaginary part')

    # drawing poles
    print("Poles ")
    print(poles)
    print("zeros")
    print(zeros)
    plt.plot(0,0,'xg',linewidth=1.2)
    plt.plot(-25,0,'xg',linewidth=1.2)
    plt.plot(-50,-10,'xg',linewidth=1.2)
    plt.plot(-50,10,'xg',linewidth=1.2)

    # drawing centroid of asymptotes
    centroid_Asymptotes = calculateCentroidOfAsymptotes(poles,zeros)
    print("centroid of asymptotes")
    print(centroid_Asymptotes)
    plt.plot(centroid_Asymptotes,0,'+k',linewidth=1.5)

    # calculate routh array and get the intersection to Imaginary axis
    RA = calculateRouthTable([1,125,5100,65000,k],0.01,True)
    sRow = RA[3,0]
    KMaxValue = sym.solve(sRow,k)
    s2Function = RA[2,0]*s**2+ KMaxValue[0]
    ImaginaryIntersectionPoints = sym.solve(s2Function,s)
    print("Intersection Points with the imaginary axis")
    print(ImaginaryIntersectionPoints)

    plt.plot(0,520**(1/2),'oc',linewidth=1)
    plt.plot(0,-520**(1/2),'oc',linewidth=1)

    # drawing angles of asymptotes
    angles_Asymptotes = AnglesOfAsymptotes(len(poles)-len(zeros))
    print("Angles of Asymptotes")
    print(angles_Asymptotes)


    plt.plot([centroid_Asymptotes ,-1*centroid_Asymptotes],[0, (abs(centroid_Asymptotes)*2*math.tan(45))],linewidth=1)
    plt.plot([centroid_Asymptotes ,-1*centroid_Asymptotes],[0 ,-(abs(centroid_Asymptotes)*2*math.tan(45))],linewidth=1);
    distance = centroid_Asymptotes*2+centroid_Asymptotes;
    plt.plot([centroid_Asymptotes, distance],[0 ,(abs(centroid_Asymptotes)*2*math.tan(45))],linewidth=1);
    plt.plot([centroid_Asymptotes, distance],[0 ,-(abs(centroid_Asymptotes)*2*math.tan(45))],linewidth=1);


    # calculate calculateBreakPoints
    BreakPoints = calculateBreakPoints(poles,zeros)
    print("Break away Points ")
    print(BreakPoints)
    plt.plot(BreakPoints[2],0,'^k',lineWidth=4)


    # tracking the poles motion 
    charEquation = charachterEquation(poles,zeros)
    coff = charEquation.coeffs()
    coff.append(0)

    coff[4] = 100
    r = np.roots(coff)
    rr = r.real
    ri = r.imag
    #realPart = r[0].real
    #imagPart = r[0].imag


    i=1000
    while i < 100000000:

        coff[4]=i
        r=np.roots(coff)
        realPart = r.real
        imagPart = r.imag
        plt.plot(realPart[0],imagPart[0],'.r',linewidth=0.01)
        plt.plot(realPart[1],imagPart[1],'.m',linewidth=0.01)
        plt.plot(realPart[2],imagPart[2],'.g',linewidth=0.01)
        plt.plot(realPart[3],imagPart[3],'.y',linewidth=0.01)
        i = i + 300000



    plt.show()





def main():
    Gui()
    
    return


if __name__ == '__main__':
   main()
