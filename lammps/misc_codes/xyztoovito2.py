import numpy




ifile = input('name trajectory file of the file-\n')
fi  = open(ifile,'r')
timestep = 0
read_line = False
out=ifile[:-4]+".defo"

while True:
  line = fi.readline()
  if len(line)==0:
    break
  if line.find('281933') != -1:
    timestep += 1
    
    with open(out,"a") as fo:
             fo.writelines("ITEM: TIMESTEP \n")
             fo.writelines("{} \n".format(timestep-1))
             fo.writelines("ITEM: NUMBER OF ATOMS\n")
             fo.writelines("281933 \n")
             fo.writelines("ITEM: BOX BOUNDS pp pp pp \n")
             fo.writelines("-100 100 \n")
             fo.writelines("-100 100 \n")
             fo.writelines("{} {}\n".format(-(timestep*(0.2/40)*250+260),(timestep*(0.2/40)*250+260)))
             print("{} {}\n".format(-(timestep*(0.2/40)*250+250),(timestep*(0.2/40)*250+250)))
             print(read_line)
             fo.writelines("ITEM: ATOMS id type x y z\n")
    continue         
  if line.find('Atoms') != -1:
 
    count = 0
    print(read_line)
    continue
  else :
    a0 = int(line.split()[0])
    a1 = float(line.split()[1])
    a2 = float(line.split()[2])
    a3 = float(line.split()[3])
        
    atom_id=count+1
    #print(line)
    
    with open(out,"a") as fo:
      fo.writelines("{} {} {} {} {} \n".format(atom_id,a0,a1,a2,a3))
      
    count += 1          
    continue
fi.close()










