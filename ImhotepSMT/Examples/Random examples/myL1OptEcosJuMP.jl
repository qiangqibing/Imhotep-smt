# before running this script, please run myVerifications.m to generate the data

using JuMP
using ECOS
using MAT

n=0
p=0

Y=zeros(n,p)
O=zeros(n,n,p)


file = matopen("./dataForECOSJulia/ECOSJulia_n100_p20.mat")
n=Int(read(file, "internal_n"))
p=Int(read(file, "internal_p"))
Y=read(file, "catenate_Y")
O=read(file, "catenate_O")
close(file)

m = Model(solver=ECOSSolver())

@variable(m, x[1:n])
@variable(m, a[1:n, 1:p]) 
@variable(m, t[1:p]>=0)

@objective(m, Min, sum(t[i] for i=1:p))

for i=1:p
@constraint(m, norm(a[:,i])<=t[i])
@constraint(m, Y[:,i].==O[:,:,i]*x+a[:,i] )
end

# print(m)
status = solve(m)

println("ecos_x = ", getvalue(x))
# println("ecos_a = ", getvalue(a))