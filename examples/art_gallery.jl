include("./plotting.jl")

using CGeometry
using CGeometry.DCELs

dcel = testradialpoly(10^2)
dcel2 = animmakemonotone(dcel, 0.1)
updatefaces!(dcel2)
dcels = explote(dcel2)
for dcel in dcels
    plotdcel(dcel)
    sleep(0.2)
end
