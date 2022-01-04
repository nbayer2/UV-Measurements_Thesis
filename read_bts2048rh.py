import os
import numpy as np
from datetime import datetime
from datetime import timedelta
# import matplotlib.pyplot as plt

def fileex(file,info=0,ver=0): #Test if a file exists
  vex=-2
  if not os.path.isfile(file):
    if info == 1:
      print("%s does not exist" % (file))
    vex=-2
  else:
    if ver >= 1:
      try:
        ftxt=open(file,"r")
        vex=0
        if ver>=2:
          try:
            llines=ftxt.readlines()
            nlin=len(llines)
            if nlin >= 1:
              vex=1
            else:
              vex=0
              if info == 1:
                print("% contains no line" % (file))
          except:
            if info == 1:
              print("%s: Not possible to read any line" % (file))
            vex=0
        else: #ver >=2
          vex=1
        ftxt.close()
      except:
        vex=-1
        if info == 1:
          print("%s may not be a txt file" % (file))
    else:# ver >=1
      vex=1
  return vex



def readtxt(ftxt_to_read):
  """
  This routine read a text file and put the information in a python dictionary
  that is returned.
  """
  dictfile={}
  dictfile["fpath"] = ftxt_to_read
  nlin=0
  try:
    ftxt=open(ftxt_to_read,"r")
    llist = ftxt.readlines()
    nlin=len(llist)
    dictfile["nlin"]=nlin
    dictfile["lineslist"]=llist
    ftxt.close()
    if nlin > 0:
      dictfile["error"] = "none"
    else:
      dictfile["error"] = "no_lines"
  except:
    dictfile["error"] = "not_readable"
    dictfile["nlin"] = nlin
    dictfile["lineslist"] = []
  return dictfile



def readmatrix(ftxt,ntop=0,nbot=0):
    dictxt=readtxt(ftxt)
    nmat=int(dictxt["nlin"]-ntop-nbot)
    lintxt=dictxt["lineslist"]
    if nmat >= 1:
      line0=lintxt[ntop].split()
      ncol=len(line0)
      if ncol >= 1:
        mat=np.zeros((nmat,ncol),dtype=float)
        #print ftxt, mat.shape
        for i in range(nmat):
          listlin=lintxt[ntop+i].split()
          #print i, listlin[0]
          ncoli= len(listlin)
          if ncoli == ncol:
            for j in range(ncol):
              mat[i,j]=float(listlin[j])
      else:
        mat=np.zeros((nmat,0),dtype=float)
        print("problem: No columns")
    else:
      print("problem, no lines or negative number of lines! for %s" % (ftxt))
      mat=np.zeros((0,0),dtype=float)
    return mat



def integband(w0,w1,mwvl,mrad,modus="norm",pos="y"):
  intband=0.
  ae0=0.094
  ae1=0.015
  be0=298.
  be1=139.
  lambe0=298.
  lambe1=328.
  lambe2=400.
  nw=mwvl.shape[0]
  for i in range(nw-1):
    mi0=mwvl[i]
    mi1=mwvl[i+1]
    ri0=mrad[i]
    ri1=mrad[i+1]
    if pos=="y":
      if ri0 < 0:
        ri0=0
      if ri1 < 0:
        ri1=0
    wvi=0.5*(mi0+mi1)
    radi=0.5*(ri0+ri1)
    if (mi0 >= w0) and (mi1 <= w1) and (ri0 >= 0) and (ri1 >=0):
      if modus == "uvind":
        if ((wvi < lambe2) and (wvi >= lambe0)):
          if wvi < lambe1:
            ae=ae0
            be=be0
          else:
            ae=ae1
            be=be1
          inti= radi*(mi1-mi0)*10.**(ae*(be-wvi))*40.
        else:
          inti=0.
      elif modus=="par":
        inti=0
      else:
        inti=radi*(mi1-mi0)
    else:
      inti=0
    intband=intband+inti
  return intband


# nhscan= maybe number of wavelength measurements per line
# nlinf= total number of lines in data file
def bts1scan(for0,iline,nlinf,nhscan,nhe):
  nmscan=221 # number of lines per wavelength measurements data block
  dscan={}
  ntop=iline+nhscan
  nbot=nlinf-ntop-nmscan+1
  m2scan=readmatrix(for0,ntop,nbot)
  m2scan1=readmatrix(for0,ntop+nmscan-1,nbot-1)
  m400=m2scan1[0][0]
  nlinb=m2scan.shape[0]
  ncolb=m2scan.shape[1]
  nwvl=int(nlinb*ncolb)+1
  mscan=np.zeros((nwvl),dtype=float)
  minfo=np.zeros((4),dtype=float)
  try:
    ntop1 = ntop-nhe-1
    nbot1 = nbot-1+nmscan+nhe
    ninfo1 = nlinf-ntop1-nbot1
    ntop2 = ntop-nhe
    nbot2 = nbot-1+nmscan+nhe-1
    ninfo2 = nlinf-ntop2-nbot2
    minf1=readmatrix(for0,ntop1,nbot1)
    minf2=readmatrix(for0,ntop2,nbot2)
    minfo[3]=minf1[0][0]
    minfo[0]=minf2[0][1]
    minfo[1]=minf1[0][1]
  except:
    dtxt=readtxt(for0)
    ltxt=dtxt["lineslist"]
    lh1=ltxt[iline+1].split()
    lh2=ltxt[iline+2].split()
    minfo[3]=int(lh1[0])
    minfo[0]=float(lh2[1])
    minfo[1]=float(lh1[1])
  msensib=1.e-8
  listscan=[]
  if nwvl >= 1:
    for ilin in range(nlinb):
      listilin=list(m2scan[ilin,:])
      listscan=listscan+listilin
    listscan.append(m400)
    mscan[:]=listscan[:]
    msup0=np.where(mscan[:] > 0.)
    msup00=msup0[0]
    nsup0=msup00.shape[0]
    if nsup0 > 0:
      msensib=min(mscan[msup00])
    if msensib <= 1.e-7:
      msensib = 0.5e-7
    if msensib <= 0.:
      msensib = 1.e-8
  minfo[2]=msensib
  dscan["mat"]=mscan
  dscan["n"]=nwvl
  dscan["info"]=minfo
  dscan["kexs"]=["n","mat","info"]
  return dscan

# for0 = filename Or0
# scans all data blocks, iteration via nscan
# nlinehead = number of header lines (in first line declared)
# nscan = number of timesteps or data blocks
# iscan = running index of nscan
# nlinscan = 226 fixed value (number of data lines per nscan) 
def bts_scans(for0,nscan,nlinhead,nlinscan,nhscan,nwvl,nhe,pref,inst="ben"):
  dscans={}
  miscans=np.zeros((nscan,2),dtype=int)
  mscans=np.zeros((nscan,nwvl),dtype=float)
  minfo=np.zeros((nscan,4),dtype=float)
  d_or0=readtxt(for0)
  lor0=d_or0["lineslist"]
  nlinf=d_or0["nlin"] # total number of lines in data file
  for iscan in range(nscan):
    iline=nlinhead+iscan*nlinscan  # first line of data block
    try:
      utcseci=int(lor0[iline])
      miscans[iscan,0]=iline # number of line of current data block
      miscans[iscan,1]=utcseci # time of current data block
      d1scan=bts1scan(for0,iline,nlinf,nhscan,nhe)
      mscani=d1scan["mat"]
      minfoi=d1scan["info"]
      mscans[iscan,:]=mscani[:]
    except:
      print("%s: Scan %i: Something wrong" % (pref,iscan))
    minfo[iscan,:]=minfoi
  dscans["info"]=miscans
  dscans["minfo"]=minfo
  dscans["mat"]=mscans
  dscans["keys"]=["info","mat","minfo"]
  return dscans



def bts_wvl(for0, nlinf):
  nheadwvl = 10
  nlwvl = 138
  dwvl = {}
  ntop = nheadwvl
  nbot = nlinf - ntop - nlwvl +1
  m2wvlb = readmatrix(for0, ntop, nbot)
  m2wvlb1 = readmatrix(for0, ntop+nlwvl-1, nbot-1)
  nlinb = m2wvlb.shape[0]
  ncolb = m2wvlb.shape[1]
  nwvlbts = int(nlinb * ncolb) +5
  mwvlb = np.zeros((nwvlbts), dtype=float)
  listwvl = []
  if nwvlbts >= 1:
    for ilin in range(nlinb):
      listilin = list(m2wvlb[ilin, :])
      listwvl = listwvl + listilin
    for i in range(m2wvlb1.shape[1]):
      listwvl.append(m2wvlb1[0][i])
    mwvlb[:] = listwvl[:]
  dwvl["mat"] = mwvlb
  dwvl["n"] = nwvlbts
  return dwvl



def bts_or0(for0,idate,pref):
  d_scans={}
  d_scans["keys"] =["listn","mwvl","mscans","miscans","minfos"]
  dicf=readtxt(for0)
  nlinf=dicf["nlin"]
  linesf=dicf["lineslist"]
  nlinmscan=121
  nlinhead=int(linesf[0].split()[0].split("=")[1])#182 # find number of header lines (first line of or0-file)
  nlinscan=226
  nhscan=5 # maybe number of wavelength measurements per line
  nwvl0=1101 # number of wavelength
  nhe = 3
  if nlinf >= 182:
    dwvl=bts_wvl(for0,nlinf)
    mwvl=dwvl["mat"]
    nwvl=dwvl["n"]
    nscan=int((nlinf-nlinhead)/nlinscan) # maybe number of timesteps of measurements
    #print(str(nscan) + " nscan")
    #print "mol4_bentham",nhscan
    if nscan >=1:
      dscans=bts_scans(for0,nscan,nlinhead,nlinscan,nhscan,nwvl,nhe, pref, inst="bts")
      miscans=dscans["info"]
      minfo=dscans["minfo"]
      mscans=dscans["mat"]
      #print("scans yes", nscan,  mscans.shape)#,mscans[int(nscan/2),:]
      #print(mscans)
    else:
      mscans=np.zeros((0,nwvl0),dtype=float)
      miscans=miscans=np.zeros((0,2),dtype=int)
    listn=[1,nlinf,nscan,nlinhead,nlinscan,nlinmscan]
    #print listn
  else:
    print("sorry, file %s is not completed (nlines < 182)" % (for0))
    mwvl=np.zeros((nwvl0),dtype=float)
    mscans=np.zeros((0,nwvl0),dtype=float)
    miscans=np.zeros((0,2),dtype=int)
    minfo=np.zeros((0,4),dtype=float)
    listn=[1,nlinf,0,0,0,nlinmscan]
  d_scans["listn"]=listn
  d_scans["mwvl"]=mwvl
  d_scans["mscans"]=mscans
  d_scans["miscans"]=miscans
  d_scans["minfo"]=minfo
  #print(d_scans)
  return d_scans




def read_oro_bts(fbtsday,keybts,i8date):
  dicbts={}
  vex=fileex(fbtsday)
  nmes=0
  nwvl=0
  nwvl0=1101#120 # number of wavelength
  muvint=np.zeros((0),dtype=float)
  muvind=np.zeros((0),dtype=float)
  if vex == 1:
    dicfbts=readtxt(fbtsday)
    nlinbts=dicfbts["nlin"]
    if nlinbts >= 1:
      nmes=int((nlinbts-1))
    else:
      vex=0
  if vex == 1:
    doro=bts_or0(fbtsday,i8date, "BTS")
    mwvl=doro["mwvl"]
    mspect=doro["mscans"]
    lnoro=doro["listn"]
    nmes=lnoro[2]
    mtime=np.zeros((nmes),dtype=float)
    mseconds=np.zeros((nmes),dtype=None)
    miscans=doro["miscans"]
    mtime[:]=miscans[:,1]/3600.
    mseconds[:]=miscans[:,1]
    nwvl=mwvl.shape[0]
    minfo=doro["minfo"]
  else:
    mspect=np.zeros((0,nwvl0),dtype=float)
    mtime=np.zeros((0),dtype=float)
    mwvl=np.zeros((nwvl0),dtype=float)
    minfo=np.zeros((4,0),dtype=float)
  if nmes > 0 and nwvl == nwvl0:
    dicbts["dataflag"]=1
    muvint=np.zeros((nmes),dtype=float)
    muvind=np.zeros((nmes),dtype=float)
    muva=np.zeros((nmes),dtype=float)
    muvb=np.zeros((nmes),dtype=float)
    for imes in range(nmes):
      # uva uvb definitions see also https://www.who.int/uv/health/solaruvradfull_180706.pdf
      mspecti=np.zeros((nwvl),dtype=float)
      mspecti[:]=mspect[imes,:]
      muvint[imes]=integband(280,400,mwvl,mspecti)
      muvb[imes] = integband(280, 315, mwvl, mspecti)
      muva[imes] = integband(315, 400, mwvl, mspecti)
      muvind[imes]=integband(280,400,mwvl,mspecti,modus="uvind")
  else:
    dicbts["dataflag"]=0
  dicbts["i8date"]=i8date
  dicbts["time"]=mtime
  dicbts["seconds"]=mseconds
  dicbts["spect"]=mspect
  dicbts["uvint"]=muvint
  dicbts["uvind"]=muvind
  dicbts["uvb"]=muvb
  dicbts["uva"]=muva
  dicbts["wvl"]=mwvl
  dicbts["n"]=[nmes,nwvl]
  dicbts["num"]= keybts
  dicbts["info"]=minfo

  # generate array of datetimes
  # source https://www.w3resource.com/python-exercises/numpy/python-numpy-datetime-exercise-4.php
  start = datetime.strptime(str(i8date), "%Y%m%d")
  mdatetimes = np.array([start + timedelta(seconds=i) for i in mseconds])
  dicbts["datetime"]=mdatetimes
  #print(mdatetimes)
  #print("WVL      " + str(dicbts["wvl"].size))
  #print("SPECT    " + str(dicbts["spect"].size))
  #print("UVA      "+str(dicbts["uva"].size))
  #print("DATETIME "+str(dicbts["datetime"].size))
  return dicbts

    

if __name__ == "__main__":
  #main()
  
  #p_bts = "/rao/x21758/level0/bts2048/global/2019/"
  p_bts = ""
  nam_bts = "MP190623rh.OR0"
  methodbts = "global"
  i8date = 20190623
  f_bts = p_bts+nam_bts
  print("Hello world!" + f_bts)
  
  d_bts1day = read_oro_bts(f_bts,methodbts,i8date)
  print(d_bts1day["i8date"])
  print(d_bts1day.keys())
  m_spect = d_bts1day["spect"]
  print(m_spect.shape)
  print(d_bts1day["n"])
  #print("No, hour, UTC, UVA, UVB, UV, UVind")
 # for i in range(d_bts1day["n"][0]):
 #    print("%3i, %2.2f, %s %.4f, %.6f, %.3f %.2f" % (i, d_bts1day["time"][i],d_bts1day["datetime"][i],d_bts1day["uva"][i],d_bts1day["uvb"][i],d_bts1day["uvint"][i],d_bts1day["uvind"][i]))
     
 # plt.plot(d_bts1day["seconds"],d_bts1day["uva"])
 # plt.show()
