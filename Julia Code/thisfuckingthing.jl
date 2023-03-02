cd("C:/Users/mcaldwe2/OneDrive - University of Wyoming/Documents/UWyo/PhD project/YNP-ungulate-niche-proj")
pwd()

#load required packages
using  Plots, Statistics, KernelDensity, DataFrames, CSV


#load data
time_df = CSV.read("./Code output/indiv_alldata_highact_times.csv", DataFrame)

#generate stupid id thing cause julia is dumb
jul_id1 = groupby(time_df, :id_yr_seas)
insertcols!(time_df, 1, :jul_id => jul_id1.groups)
#create unique id string
ids = unique(time_df[:, "id_yr_seas"])

#loop by individual
for i in 1:(length(ids))
    a = subset(time_df, :id_yr_seas => ids[i])

    a = time_df[(time_df.id_yr_seas => ids[i]), :]
#create kernel densities from data. It returns an object that you must extract the density from
a=kde(a)
b=kde(b)

#replace objects with the kernal density only
a=a.density
b=b.density

min_b=minimum(b)
max_b=maximum(b)
max_a=maximum(a)
min_a=minimum(a)


#initialize empty overlap vector
overlap = Vector{Float64}()

#collect values of overlap
for i in 1:(length(a))
    if a[i]>b[i] && a[i]<max_a
        append!(overlap,b[i])
    elseif a[i]<b[i] && a[i]<max_a
        append!(overlap,a[i])
    end
end

#plot histograms to check if overlap is accurate
plot(a,fill=(0,:green), opacity=[0.30],legend=false)
plot!(b,fill=(0,:blue), opacity=[0.30],legend=false)
plot!(overlap,fill=(0,:red), opacity=[0.3],legend=false)

#gen sums
sum_overlap=sum(overlap)
sum_a=sum(a)
sum_b=sum(b)

#gen ratios of overlap