#!/usr/bin/python 
import os,re
from shutil import copyfile,copy
import sys
from numpy import *
from numpy.linalg import *
from optparse import OptionParser
import subprocess as sp

### Default parameters ###
pars = OptionParser()
pars.add_option("-i", default = "start.gro", dest = "iname", \
        help = "input file (.gro) [default: start.gro]")
pars.add_option("-o", default = "rot", dest = "oname", \
        help = "output file, just name without extension [default: rot]")
pars.add_option("-p", default = "temp.top", dest = "top", \
        help = "temporary topology file [default: temp.top]")
pars.add_option("--em", default = "em.mdp", dest = "em_mdp", 
        help = "protocol file for energy minimization [default: em.mdp]")
pars.add_option("--ru", default = "run.mdp", dest = "run_mdp", \
        help = "protocol file for md simulation [default: run.mdp]")
pars.add_option("--at", default = "dih.ndx", dest = "ind", \
    help = "index file with atom numbers of the dihedral [default: dih.ndx]")
pars.add_option("--np", default = 1, dest = "np", \
        help = "number of CPUs to use for mdrun [default: 1]")
pars.add_option("--ns", default = 36, dest = "ns", \
        help = "number of positions to scan [default: 36]")
pars.add_option("--dir", default = "./TMP/", dest = "drct", \
        help = "directory for temporary files [default: ./TMP]")
pars.add_option("--gr", default = "gmx", dest = "grom",\
        help = "path to the installed version of gromacs")

improper_barrier = 10000 # the value of the parameter for the improper potential
CONNECT_BARRIER = 0.16   # the maximal distance for attribution of atoms to the rotating 
                         # or staying group. CONNECT_BARRIER must be increased if "atom numbers do not cover the 
                         # molecule", and decreased if "trying to rotate a cycle molecule"

#############################################
############## Functions ####################
#############################################

class StdClass:
    pass

def CreateDir(_scrDir, _topFile):
  try:
    os.mkdir(_scrDir)
    tl,topFileName = os.path.split(_topFile)
    topRelPath = os.path.join(_scrDir, topFileName)
    copyfile(_topFile, topRelPath)
  except Exception as e:
    print '==> Problem in scratch directory creation <=='
    raise e
  return topRelPath
 

def ReadGro(gro):
  coord = []
  atom = []
  F = open(gro)
  line = F.readline() # read the first line of .gro file
  natom = int(float(F.readline())) # get the number of atoms in the molecule from the second line of the .gro file
  for p in range(natom): # create an array of atom cooridnates
    line = F.readline()
    coord.append([float(line[22:29]),float(line[30:37]),float(line[38:45])]) # check that positions of atom cooridnates
                                                                             # are given correctly!!!
    #coord.append([float(line[22:35]),float(line[38:51]),float(line[54:67])])
    atom.append(line[0:21]) # save all the information about atoms (except for atom coordinates)
  coord = array(coord)
  F.close()
  return atom, coord


def MoveCoord(coord, inds):
  """ Move second atom in opt.at to zero point  """
  coord = coord - coord[inds[1],:]
  ### Orient the roatation bond (between 2nd and 3rd atom) in the direction of the Oz axis ###
  # rotation OX #
  COSA = coord[inds[2],2]/norm(coord[inds[2],[1,2]])
  SINA = coord[inds[2],1]/norm(coord[inds[2],[1,2]])
  MX = matrix([[1,0,0],[0,COSA,SINA],[0,-SINA,COSA]])
  coord = array(coord*MX)
  if abs(coord[inds[2],1]) > 10**(-9):
    print '################\n#\n#'
    print '# !!!ERROR!!!'
    print '# WRONG ROTATION MATRIX OX!!!'
    sys.exit(1)

  # rotation OY #
  COSA = coord[inds[2],2]/norm(coord[inds[2],[0,2]])
  SINA = coord[inds[2],0]/norm(coord[inds[2],[0,2]])
  MY = matrix([[COSA,0,SINA],[0,1,0],[-SINA,0,COSA]])
  coord = array(coord*MY)
  if abs(coord[inds[2],0]) > 10**(-9):
    print '################\n#\n#'
    print '# !!!EROR!!!'
    print '# WRONG ROTATION MATRIX OY!!!'
    sys.exit(1)
  return coord


def DetermineRotPart(coord, CONNECT_BARRIER, inds):
  """ Determine the two parts of the molecule: 
  which one rotates and which one stays
  build a connection matrix """

  M1 = zeros([len(coord),len(coord)])
  M2 = zeros([len(coord),len(coord)])
  for p in range(len(coord)):
    for q in range(p,len(coord)):
      if ( sqrt(((coord[p,:]-coord[q,:])**2).sum()) < \
            CONNECT_BARRIER ) and \
         ( p != inds[1] or q != inds[2] ) and \
         ( p != inds[2] or q != inds[1] ):
        M1[p,q] = 1
        M1[q,p] = 1

  # detect rotating and staying groups #
  rot = list(M1[inds[1],:].nonzero()[0])
  stop = 0
  while stop != 1:
    stop = 1
    for p in rot:
      for q in list(M1[p,:].nonzero()[0]):
        if q not in rot:
          rot.append(q)
          stop = 0
  stop = 0
  stay = list(M1[inds[2],:].nonzero()[0])
  while stop != 1:
    stop = 1
    for p in stay:
      for q in list(M1[p,:].nonzero()[0]):
        if q not in stay:
          stay.append(q)
          stop = 0

  # check that all atoms are attributed either to the rotating or to staying group #
  if len(rot) == len(coord):   
    print '!!! ERROR !!!\n TRYING TO ROTATE A CYCLE MOLECULE !!! Decrease the value of connectivity barrier!!! '
    exit(1)
  elif len(rot) + len(stay) != len(coord):
    print '!!! ERROR !!!\n SOMETHING WRONG!!! ROTATING AND STAYING GROUPS DO NOT COVER THE MOLECULE!!! Increase the value of connectivity barrier!'
    exit(1)
  return rot


def WriteGro(_coords, _atomLines, _scrDir):
  # write a .gro file with correct coordinates #
  W = open(os.path.join(_scrDir,'TMP.gro'),'w') 
  W.write("File created by ProfileMD.py, TPP project \n%7d\n" % _coords.shape[0])
  for k in range(_coords.shape[0]):
    W.write(_atomLines[k]+ "%13.10f %13.10f %13.10f\n" % tuple(_coords[k,:]) )
  W.write('   6.00000   6.00000   6.00000') # set box
  W.close()
  return




def MinimizeEnergy(_grom, _TMP, _scrDir, _emmdp, _topfile, _nt):
  try:
    exlist = [_grom,'grompp','-f',_emmdp,'-c',\
          os.path.join(_scrDir,'TMP.gro'),\
          '-p',_topfile,'-o',_TMP.tpr,\
          '-po',os.path.join(_scrDir,'mdout_em.mdp')]
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE )
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running gmx grompp'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
    print sout
    raise Exception("Cannot create .tpr file for energy minimization! Check out GROMACS log above.")

  try:
    exlist = [_grom,'mdrun','-s', _TMP.tpr,'-c', _TMP.gro,\
      '-nt',str(_nt),'-g',os.path.join(_scrDir,'mdem.log'),\
      '-e',os.path.join(_scrDir,'emener.edr'),\
      '-o',os.path.join(_scrDir,'emtraj.trr')]
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE )
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running gmx mdrun'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
    print sout
    raise Exception("Cannot run energy minimization! Check out GROMACS log above.")

  return 
 

def MDrun(_grom, _TMP, _scrDir, _runmdp, _topFile, _nt):
  try:
    exlist = [_grom,'grompp','-f',_runmdp,'-c',_TMP.gro,\
          '-p',_topFile,'-o',_TMP.runtpr,\
          '-po',os.path.join(_scrDir,'mdout_ru.mdp')]
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE )
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running gmx grompp'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
    print sout
    raise Exception("Cannot create .tpr file for relaxation! Check out GROMACS log above.")

  try:
    exlist = [_grom,'mdrun','-s',_TMP.runtpr,'-e',_TMP.edr,'-deffnm',_TMP.nm,\
      '-nt', str(_nt),'-g',os.path.join(_scrDir,'mdrun.log')]
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE )
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running gmx mdrun'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
    print sout
    raise Exception("Cannot run relaxation! Check out GROMACS log above.")
  return


def GetEnergy(_grom, _TMP, extractTerms=False):
  try:
    exlist = [_grom,'energy','-f',_TMP.edr,'-o',_TMP.xvg,'-xvg',\
        'xmgrace' if extractTerms else 'none']
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE, stdin=sp.PIPE )
    zz.stdin.write('1 2 3 4 5 6 7 8 9 10 11 12 13 14 \n')
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running gmx energy'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
    print sout
    raise Exception("Cannot extract energy terms! Check out GROMACS log above.")
  zz.stdin.close() 
  if extractTerms:
    f = open(_TMP.xvg, 'r')
    terms = []
    eners = []
    # extracting terms
    for li in f:
      mm = re.match(r'\@ s[0-9]+ legend "([^"]+)"', li)
      if mm:
        terms.append(mm.group(1))
      else:
        if li[0] != '@':
          eners.append(li)
    f.close()
    # rewrite file without terms
    f = open(_TMP.xvg, 'w')
    f.writelines(eners)
    f.close()
    return terms
  return []


def CleanMemory(_scrDir):
  print "Clearing temporary files..",
  sys.stdout.flush()
  os.unlink(os.path.join(_scrDir,'emener.edr'))
  os.unlink(os.path.join(_scrDir,'mdem.log'))
  os.unlink(os.path.join(_scrDir,'mdrun.log'))
  os.unlink(os.path.join(_scrDir,'mdout_ru.mdp'))
  os.unlink(os.path.join(_scrDir,'mdout_em.mdp'))
  os.unlink(os.path.join(_scrDir,'emtraj.trr'))
  os.unlink(os.path.join(_scrDir,'angdist.xvg'))
  for f in os.listdir(_scrDir):
    if f[0] == '#':
        os.unlink(_scrDir+f)
  print "DONE"
  return


def SetAngle(_grom, _scrDir, _ind):
  """ calculate angles """
  #TODO: refactor
  try:
    exlist = [_grom,'angle','-f',\
          os.path.join(_scrDir,'TMP.gro'),\
          '-n',_ind,'-type','dihedral',\
          '-ov',os.path.join(_scrDir,'TMP_ang.xvg'),\
          '-od',os.path.join(_scrDir,'angdist.xvg'),\
          '-all','-xvg', 'none']
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE )
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running gmx angle'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
    print sout
    raise Exception("Cannot calculate angle value! Check out GROMACS log above.")
  angle = loadtxt(os.path.join(_scrDir,'TMP_ang.xvg'))[1]
  return angle
 

def TopologyUpdate(_svTop, _topFile, angle):
  # fix the angle in the topology file #
  FITP = open(_topFile, 'w')
  _svTop[9] = "#define mydih    2  %8.3f  %8.3f\n" % (angle, improper_barrier) # check that 4th line in the topology file was empty!
  # write the resulting topology file #
  for ln in _svTop:
    FITP.write(ln)
  FITP.close()
  return
  

def WriteFrame(_grom, p, _aniFile, _scrDir):
  """ Write current frame """
  try:
    exlist = [_grom,'editconf','-f',\
        os.path.join(_scrDir,'TMP.gro'),'-o',\
        os.path.join(_scrDir,'TMP.pdb')]
    zz = sp.Popen( exlist, stderr=sp.PIPE, stdout=sp.PIPE )
    serr,sout = zz.communicate()
  except OSError as e:
    print '** Error running editconf'
    print '** The command was: \n' + ' '.join(exlist) + '\n'
    raise e

  if zz.returncode != 0:
      print sout
      raise Exception("Cannot convert .gro to .pdb! Check out GROMACS log above.")
  T = open(os.path.join(_scrDir,'TMP.pdb'),'r')
  for i in range(3):
      _aniFile.write(T.readline())
  _aniFile.write('MODEL        %5d\n' % (p))
  T.readline()
  while True:
      ln = T.readline()
      if not ln:
          break
      _aniFile.write(ln)
  T.close()      
  return


def WriteEnergy(p, angle, _TMP, _E):  
  # calculate average potential energy value #
  en = loadtxt(_TMP.xvg)[:,1:] # read only energy terms
  n_ener = len(en) # the number of energy terms written in .edr file
  energy = mean(en[int(n_ener*0.5):, :],0) # average energy over the second part of the trajectory
  _E.write('%15.3f ' % angle)
  for i in range(14):
    _E.write('%15.3f ' % energy[i])
  _E.write('\n')
  return


###########################################################
################## Main algorithm #########################
###########################################################

if __name__ == "__main__":
  ### Set necessary parameters ###
  opt, args = pars.parse_args()

  print "Initial structure: ", opt.iname
  print "GMX topology: ", opt.top
  print "GMX emin protocol: ", opt.em_mdp
  print "GMX run protocol: ", opt.run_mdp
  print "Index file: ", opt.ind
  print "Number of threads to run GMX: ", opt.np
  print "Number of points to scan: ", opt.ns
  print "Directory of saving simulations: ", opt.drct
  print "GROMACS run: ", opt.grom

  params = StdClass()
  params.iname = opt.iname
  params.oname = opt.oname
  params.top = opt.top
  params.em_mdp = opt.em_mdp
  params.run_mdp = opt.run_mdp
  params.ind = opt.ind
  params.np = int(opt.np)
  params.ns = int(opt.ns)
  params.scrDir = opt.drct
  params.grom = opt.grom

  del opt

  print "================================== "
  print ""

  print "Loading Index ..",
  sys.stdout.flush()
  at = loadtxt(params.ind, skiprows=1, dtype=int64) 
  if len(at) != 4:
    print '!!! ERROR!!! NUMBER OF ATOMS IS NOT 4!!!!'
    sys.exit(1)
  else:
    rotIdx = int32(at)-1 
    # attribute positions in the coordinate array to atom nubers:
    # first atom has position 0, second --- 1, ... etc
  print "DONE"
  print "DIHEDRAL QUAD: ",at

  ### Create a rotation matrix  ###
  MZ = matrix([[cos(2*pi/params.ns),sin(2*pi/params.ns),0],\
          [-sin(2*pi/params.ns),cos(2*pi/params.ns),0],[0,0,1]])

  ### Read the topology file ###
  print "Try reading the topology ..",
  sys.stdout.flush()
  F_temp = open(params.top, 'r')
  savedTopology = F_temp.readlines()
  F_temp.close()
  print "DONE"
    
  ### Prepare pdb-file with the result of rotation ###
  print "Preparing animated PDB",
  sys.stdout.flush()
  S = open('_'+params.oname+'-rotated.pdb','w')
  S.write('REMARK  File created by ProfileMD.py, TPP project \n')
  print "DONE"

  ### Prepare file to write energies ###
  print "Preparing energy file",
  sys.stdout.flush()
  E = open('_energies.xvg','w')
  E.write('# Here we will write all energy terms\n')
  print "DONE"

  ### Parsing GRO file
  print "Reading GRO file..",
  atomLines, atomCoords = ReadGro(params.iname)
  print atomCoords.shape[0], " atoms"

  atomCoords = MoveCoord(atomCoords, rotIdx)
  rot = DetermineRotPart(atomCoords, CONNECT_BARRIER, rotIdx)

  ### Creating temporary directory
  print "Creating scratch directory..",
  sys.stdout.flush()
  if params.scrDir[-1] != '/':
      params.scrDir += '/'
  try:
    topFile = CreateDir(params.scrDir, params.top)
  except OSError as e:
    if e.errno == 17:
      print '** ERROR: directory exists!'
      print '** May be you did not remove scratch directory %s after previous run?' \
          % (params.scrDir)
    raise e
  print "DONE"

  for p in range(params.ns):
    print " ==> Processing point #", p
    WriteGro(atomCoords, atomLines, params.scrDir)
    WriteFrame(params.grom, p, S, params.scrDir)
    angle = SetAngle(params.grom, params.scrDir, params.ind)
    print "Current angle is: ", angle
    TopologyUpdate(savedTopology, topFile, angle)

    TMP = StdClass()
    TMP.tpr = os.path.join(params.scrDir,'TMP'+str(p)+'.tpr')
    TMP.gro = os.path.join(params.scrDir,'TMP'+str(p)+'.gro')
    TMP.edr = os.path.join(params.scrDir,'TMP'+str(p)+'.edr')
    TMP.nm =  os.path.join(params.scrDir,'TMP'+str(p)+'nm')
    TMP.runtpr = os.path.join(params.scrDir,'TMP'+str(p)+'run.tpr')
    TMP.xvg = os.path.join(params.scrDir,'TMP'+str(p)+'.xvg')

    print "Minimizing energy..",
    sys.stdout.flush()
    MinimizeEnergy(params.grom, TMP, params.scrDir, params.em_mdp, topFile, params.np)
    print "DONE"

    print "Making small dynamics..",
    sys.stdout.flush()
    MDrun(params.grom, TMP, params.scrDir,params.run_mdp,topFile, params.np)
    print "DONE"

    print "Extracting energy terms..",
    sys.stdout.flush()
    if p == 0:
      termLabels = GetEnergy(params.grom, TMP, extractTerms=True)
      termLabels = ['Dih.value'] + termLabels
      E.write("#" + "|".join([ "%14s " % i for i in termLabels]) + "\n")
    else:
      GetEnergy(params.grom, TMP, extractTerms=False)
    print "DONE"

    skips = WriteEnergy(p,angle, TMP, E)
  
    # rotate molecule and pass to the next angle position #  
    atomCoords[rot] = array(atomCoords[rot]*MZ)
    CleanMemory(params.scrDir)
  S.close()
  E.close()

