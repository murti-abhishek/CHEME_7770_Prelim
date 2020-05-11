IPTG = [0, 5e-4, 0.005, 0.012, 0.053, 0.216, 1] # mM
mRNA = [19, 21, 41, 67, 86, 93, 93] # <n> number of mRNA copies / cell
mRNA_low = [18, 17, 37, 65, 84, 91, 92] # low
mRNA_high = [20, 26, 44, 69, 88, 95, 94] # high

using Plots
#plot(IPTG,mRNA)

# ----------------- Parameters -----------------

av_number = 6.02e23      # number/mol
mass_of_single_cell = 1e-12 # g Bionumbers : 101789  #2.8e-13
fraction_dry_cell = 0.27 # Bionumbers : 110086
RNAPII_copy_number = 4600 #1150 # copies/cell # Bionumbers : 108601
R_X_T = RNAPII_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

# ------------- Kinetic transcription limit -------------
gene_length = 1000 # nt #given
transcription_elongation_rate = 50 # nt/s # Bionumbers : 111871
k_E_x = transcription_elongation_rate*(1/gene_length)

characteristic_initiation_time = 42 # s #https://github.com/varnerlab/JuGRN-Generator/blob/master/src/distribution/Default.json
kI = (1/characteristic_initiation_time)
τₓ = k_E_x/kI
Kₓ = 0.24 #nmol/gDW # https://github.com/varnerlab/JuGRN-Generator/blob/master/src/distribution/Default.json

# ------------- (Hill function) -------------
K_D = 49.6e-3 # Bionumbers : 101976
n = 0.9 # Important parameter
f = zeros(length(IPTG))

# ------------- (Control function) -------------
W1 = 0.01 # Important parameter
W2 = 0.05 # Important parameter
u = zeros(length(IPTG))

# ------------- Dilution and degradation -------------

cell_doubling_time = 40 / 60 # hr
μ = log(2) / cell_doubling_time # hr^-1
mRNA_half_life = 5 / 60 # hr
θ = log(2) / mRNA_half_life # hr^-1

# ------------- Model definition ------------
gene_copy_number = 2
G = gene_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell)   # nmol/gDW

r_x = k_E_x * R_X_T * (G / (τₓ * Kₓ + (τₓ + 1) * G)) * 3600

m = zeros(length(IPTG))

for i = 1:length(IPTG)

    f[i] = IPTG[i]^n / (IPTG[i]^n + K_D^n)
    u[i] = (W1 + W2 * f[i]) / (1 + W1 + W2 * f[i])
    m[i] = (r_x * u[i]) / (μ + θ)

end


plot(IPTG[2:end], m[2:end],linewidth=3, xaxis = ("IPTG (mM)", :log),ylabel = "lacZ mRNA (nmol/gDW)",label="Model",legend=:bottomright)

mRNA_new = mRNA*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW
mRNA_low_new = mRNA_low*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW
mRNA_high_new = mRNA_high*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW
plot!(IPTG[2:end], mRNA_new[2:end],label="EXP")
plot!(IPTG[2:end], mRNA_low_new[2:end],label="low")
plot!(IPTG[2:end], mRNA_high_new[2:end],label="high")
savefig("Problem_1.png")
