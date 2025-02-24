include("reactChem.jl")
include("reactStatsPlots.jl")
using StatsBase
using DelimitedFiles
using Distributions

wpInitsArray = [0.783,0.822,0.858,0.887,0.893,0.909,0.919,0.926,0.936]          # wp initialization for each layer
layerRadisArray = [0,51,57.19,62.11,66.21,69.76,72.89,75.72,78.31,80.68,82.88]  # radius of each layer
Tgvalues = [25.4,45.1,102.1]                                                    # Tg values
colorpalletes = palette(:cool, length(Tgvalues))
N = 1000                                                                        # N : NumberOfRadical
Vs = fill(0.0,(length(Tgvalues),length(wpInitsArray),N))                        # Depth Radius (Deepest)
for (ntg,tgv) in enumerate(Tgvalues)                                            # Iterate over TgValues
    for (nwp,wpInit) in enumerate(wpInitsArray)                                     # Iterate over wp : layer increase with wp
        # Initialize Multi Layer Particle variables
        wpEnd = 1.0                                                                 # initial Wp values
        startWpTime , endWpTime = 0.0 , 30*60.0                                     # linear Relationship Wp reaction end Time
        parRadius = layerRadisArray[(2+nwp)]                                        # particle radius
        reactionTemp = 70                                                           # T : Reaction Temparature (°C)
        layerRadius = layerRadisArray[[1]]                                          # Layer radius 
        Tg₀ = fill(tgv,size(layerRadius))                                           # initial Tg
        Tmon = 106                                                                  # monomer temp 106°C

        # --- Make particle observable function --- 
        layerR = Observable(layerRadius)                                            # layerR : layer Radius away from center)
        Wp = Observable(wpInit)                                                     # Wp : Weightage of polymer
        layerΔT = Observable(reactionTemp .- Tg₀)                                   # layerΔT = T - Tg # (@lift(reactionTemp .- updateTg($Wp,Tg₀,Tmon)) )
        layerD = @lift([10^logD($Wp,t,unit="nm") for t in $layerΔT])                # Layers Diffusion coefficient
        par = mLparticle(parRadius,Wp,layerR,layerΔT,layerD)                        # create particle with (Radius,Wp,Layers,Temp,Diffusion)        

        # Simulation propagation Time statistics
        lPropl = 297                                                                # Number of propagation steps before simulation end
        propTime = @lift(propagationTimeInterval($(par.Wp),reactionTemp))           # propTime : Time Interval for Monomer to propagate
        propStats = getPropStats(wpInit,reactionTemp,endWpTime,lPropl,              # get propagation statistics
                                wpEnd = wpEnd, startTime=startWpTime)                   
        T = propStats.T                                                             # T        : Total simulation time interval

        # initialize Radicles
        zmerInit = 4                                                                # Initial zmer length
        # maxStepLength = 50          # Alternative: minimum(diff(par.layerR[]))      # stepSizeLength: Smaller length will discretize time longer
        zmer = Observable(zmerInit)                                                 # Make zmer an observable
        τ =   Observable(T/1000) # Observable(MinTimeForStepsize(maxStepLength,par.layerD[],wpInit,zmerInit,confident=0.9)  ) #
        radRadius = 0.75                                                            # initialize radical Radius
        initLoc = par.obj.radius-radRadius*(1+1e-10)                                # initial radius location
        Rad = [Radicle([initLoc,0,0],par,r=radRadius,zmer = zmer,τ = τ) for n in 1:N] # Initialize N radicals
        for rad in Rad                                      
            random_position(rad.obj)                                                # start from different edge of the sphere
        end

        sims = simulate(propTime,Wp,zmer,collect(0:τ[]:T),par,Rad,wpInit,endWpTime, # run simulations
        wpEnd=wpEnd,startTime=startWpTime,wallCollideDo="random",eachStepDo ="random",seed=1)    
        anim3D(sims,par,radicalRadius=radRadius,fps=60,videoName=wdir("animation/tg-$tgv/")*"wp-$wpInit.mp4",radC=colorpalletes[[ntg]])
        anim2D(sims,par,τ,radicalRadius=radRadius,fps=60,videoName=wdir("animation/2D/tg-$tgv/")*"wp-$wpInit.mp4",axislim=85)
        info = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)             #
        v = info.sortedValues[size(info.sortedValues)[1],:]
        num_bins = max(1,Int(ceil(maximum(v) - minimum(v))))
        his = Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[length(layerRadisArray)],
        ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,layerRadisArray[length(layerRadisArray)]), xlabel="Distance (nm) from center", 
        ylabel="Probability", title="Histogram of closest depth reach by radicles", legend=false)
        Vs[ntg,nwp,:] = info.sortedValues[size(info.sortedValues)[1],:]
        # Add a rectangle (shaded region)
        x1,x2 = parRadius,ceil(layerRadisArray[end])
        y1,y2 = 0,1  
        Plots.plot!(his,[x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], fillrange=0, fillalpha=0.3, color=:black, label="Shaded Region")
        savefig(his,wdir("plots/DeepestHistogram/tg-$tgv/")*"wp-$wpInit.png")
    end
end
