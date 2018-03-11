# before running this script, please run myVerifications.m to generate the data

import ECOS
using MAT
println(ECOS.ver())

system_n=0
system_p=0

system_Y=zeros(system_n,system_p)
system_O=zeros(system_n,system_n,system_p)

file = matopen("./dataForECOSJulia/ECOSJulia_n10_p20.mat")
system_n=Int(read(file, "internal_n"))
system_p=Int(read(file, "internal_p"))
system_Y=read(file, "catenate_Y")
system_O=read(file, "catenate_O")
close(file)

n = system_p+system_n*system_p +system_n
m = 1*system_p+(1+system_n)*system_p
p = system_n*system_p
l = system_p 
ncones = system_p 
q = floor.(Int, reshape((1+system_n)*ones(1,system_p),system_p))

# Construct matrix G

Gpre1= [-1*eye(system_p) zeros(system_p, system_n*system_p+system_n)] # this matrix is used to extract -1*[t1, \ldots, tp] from the optimization variable
Gpre2= [ zeros(system_n*system_p, system_p) -1*eye(system_n*system_p) zeros(system_n*system_p, system_n)] # this matrix is used to extract -1*[a1, \ldots, ap] from the optimization variable
Gpre3=zeros((1+system_n)*system_p,n) # Gpre3 is used to extract -1*[t1, a1, t2, a2, \ldots, tp, ap] from the optimization variable, the method is by merging Grep1 and Grep2
for index=1:system_p
    Gpre3[(1+(index-1)*(system_n+1)):index*(system_n+1),:]=[reshape(Gpre1[index,:],1,n);
                      Gpre2[(1+(index-1)*system_n):index*system_n,:]]
end
G=[-1*eye(system_p) zeros(system_p, system_n*system_p+system_n); 
   Gpre3] 

# Construct matrix A
vec_O=zeros(system_p*system_n,system_n)
for index=1:system_p
    vec_O[1+(index-1)*system_n:index*system_n, :]=system_O[:,:,index]   
end

A=[zeros(system_n*system_p, system_p) eye(system_n*system_p) vec_O]

c = [ones(1, system_p) zeros(1, system_n*system_p) zeros(1, system_n)]

h = zeros(m ,1) 

Y=zeros(system_n*system_p,1)
for index=1:system_p
    Y[1+(index-1)*system_n:index*system_n,:]=system_Y[:,index]
end

b=Y

ecos_A = ECOS.ECOSMatrix(float(A))
ecos_G = ECOS.ECOSMatrix(float(G))
prob = ECOS.setup(n, m, p, l, ncones, q, 0,  ecos_G, ecos_A, reshape(c, length(c)), reshape(h, length(h)), reshape(b, length(b)))

exitflag = ECOS.solve(prob)
