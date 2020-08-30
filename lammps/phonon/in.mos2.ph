# 3D copper block simulation
boundary     p p p
units        metal
atom_style   atomic

# geometry
read_data	 datasw13.pos

mass            1 32.066
mass            2 95.94
mass            3 32.066
mass            4 32.066
mass            5 95.94
mass            6 32.066



#4 potentials
pair_style      hybrid/overlay sw sw lj/cut 10.0 &
kolmogorov/crespi/z 14.0  kolmogorov/crespi/z 14.0 kolmogorov/crespi/z 14.0 kolmogorov/crespi/z 14.0 


# SW for layer 1
pair_coeff      * *  sw 1 mos216.sw S Mo S NULL NULL NULL 
# SW for layer 2
pair_coeff      * *  sw 2 mos216.sw NULL NULL NULL S Mo S
# Set arbitrary interactions - 0, before overlaying with KC
pair_coeff * * lj/cut 0.0 3.4

#------------------------------------------
# S layer1: 1, 4, 7, 10=3
# Mo layer1: 2, 5, 8, 11=2
# S layer1: 3, 6, 9, 12=1
# S layer2: 15, 18, 21, 24=4
# Mo layer2: 14, 17, 20, 23=5
# S layer2: 13, 16, 19, 22=6

#------------------------------------------
# S-S KC interactions.
#------------------------------------------

# KC between LAYER 3 and LAYER 4
pair_coeff 3 4 kolmogorov/crespi/z 1 ILP.KC NULL NULL S S NULL NULL 

#------------------------------------------
# S-Mo KC interactions.
#------------------------------------------

# KC between LAYER 2 and LAYER 4
pair_coeff 2 4 kolmogorov/crespi/z 2  ILP.KC NULL Mo NULL S NULL NULL
# KC between LAYER 3 and LAYER 5
pair_coeff 3 5 kolmogorov/crespi/z 3  ILP.KC NULL NULL S NULL Mo NULL 


#------------------------------------------
# Mo-Mo KC interactions.
#------------------------------------------
# Mo layer1: 2, 5, 8, 11 = 2
# Mo layer2: 14, 17, 20, 23 = 5

# KC between LAYER 2 and LAYER 5
pair_coeff 2 5 kolmogorov/crespi/z 4  ILP.KC NULL Mo NULL NULL Mo NULL 
neighbor     2. nsq
neigh_modify every 1 delay 0 check yes

#Langevin random seed
variable     dt equal 2e-3
variable     r  equal 57085
variable     T  equal 300
variable     dT equal "v_dt * 100"

timestep ${dt}

# initialize
velocity     all create $T 28459 rot yes dist gaussian mom yes
reset_timestep 0

# fixes 
fix          1 all langevin $T $T ${dT} 73504 zero yes
fix          2 all nve
fix          3 all phonon 10 50000 1000000 map.in mos2 nasr 50

# output
#                    1    2    3  4  5     6   7
thermo_style custom step temp pe ke press pxx pyy
thermo       100

restart      1000000 restart.one restart.two

dump   1 all xyz 50000 dump.mos2.xyz

# execution
run 	 6000000
write_restart Restart.final
