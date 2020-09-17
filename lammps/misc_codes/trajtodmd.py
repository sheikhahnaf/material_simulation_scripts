import numpy




ifile = input('name trajectory file of the file-\n')
fi  = open(ifile,'r')

m = 0
n = 0
read_time = False
read_natoms = False

while True:
  line = fi.readline()
  if len(line)==0:
    break
  if line.find('ITEM: TIMESTEP') != -1:
    read_time = True
    continue
  if read_time == True:
    n += 1
    a0 = int(line.split()[0])
    if n==1:
      step1 = a0
    if n==2:
      step2 = a0
      t = step2 - step1   #  Output interval
    nsteps = a0
    read_time = False
    continue
  if line.find('ITEM: NUMBER OF ATOMS') != -1:
    m += 1
    if m==1:
      read_natoms = True
    continue
  if read_natoms == True:
    a0 = int(line.split()[0])
    natoms = a0
    read_natoms = False
    continue
fi.close()

nmax = n+1    # NO. of data

colx = 2
coly = 3
colz = 4
colvx = 5
colvy = 6
colvz = 7





x  = numpy.array([[0.0 for i in range(nmax)] for j in range(natoms)])
y  = numpy.array([[0.0 for i in range(nmax)] for j in range(natoms)])
z  = numpy.array([[0.0 for i in range(nmax)] for j in range(natoms)])
vx  = numpy.array([[0.0 for i in range(nmax)] for j in range(natoms)])
vy  = numpy.array([[0.0 for i in range(nmax)] for j in range(natoms)])
vz  = numpy.array([[0.0 for i in range(nmax)] for j in range(natoms)])

fi = open(ifile,'r')
t = -1
read_data = False
while True:
  line = fi.readline()
  if len(line)==0:
    break
  if line.find('ITEM: TIMESTEP') != -1:
    t += 1 
    read_data = False
    continue
  if line.find('ITEM: ATOMS') != -1:
    read_data = True
    continue
  if read_data == True:
      a0 = int(line.split()[0])
      a1 = float(line.split()[colx])
      a2 = float(line.split()[coly])
      a3 = float(line.split()[colz])
      a4 = float(line.split()[colvx])
      a5 = float(line.split()[colvy])
      a6 = float(line.split()[colvz])
      id=a0-1
      x[id,t] = a1
      y[id,t]= a2
      z[id,t] = a3
      vx[id,t] = a4
      vy[id,t]= a5
      vz[id,t] = a6
      continue
A= numpy.vstack((x,y,z,vx,vy,vz))
numpy.savetxt("dmd.txt",A.T, delimiter=",")
print(A.shape)


fi.close()




