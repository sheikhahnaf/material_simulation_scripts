# 3D copper block simulation
boundary     p p p
units        metal
atom_style   atomic


read_data	mos2x2.pos
mass            1 95.94
mass            2 32.066



pair_style sw 


# SW for layer 1
pair_coeff      * * mos216.sw Mo S 

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
fix          3 all phonon 10 50000 1000000 mos2x2.map.in mos2.1l nasr 50

# output
#                    1    2    3  4  5     6   7
thermo_style custom step temp pe ke press pxx pyy
thermo       100

restart      1000000 restart.one restart.two

dump   1 all xyz 50000 dump.mos2.1l.xyz

# execution
run 	 6000000
write_restart Restart.final
