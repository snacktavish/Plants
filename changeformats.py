import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import gzip
from subprocess import call
class Chromosome:
  def __init__(self,infii):
       self.infile=infii
       fi=open(infii)
       header=fi.readline().split()
       self.numinds=len(header[3:])
       self.indivs={index:elem for index, elem in enumerate(header[3:])}
       self.genos=dict.fromkeys(range(len(header[3:])),None)
       self.multihit=[]
       for item in self.genos: 
          self.genos[item]={}
       self.snplist=[]
       self.nsnp=0
       for num,lii in enumerate(fi):
          lin=lii.strip().split('\t')
          snp=lin[1]
          if len({x.strip('/') for x in lin[2:] if x not in ['./','./.']})==2:
            if num==0:
                self.chrm=lin[0]
            self.snplist.append(snp)
            ref=lin[2]
            for ii,item in enumerate(lin[3:]):
                self.genos[ii][snp]=self.translate(ref,item.split('/')[0])
            self.nsnp+=1
          if len({x.strip('/') for x in lin[2:] if x != './'})>2:
#            print("Snp %s has more than 2 nucleotides" %snp) 
            self.multihit.append(snp)  
          if num%1000000==0:
            print(num)
 #      self.tripletest()
  def translate(self,ref,val): #this is setting allele in refernce as 0
    assert ref in ['A','T','G','C'], "ref is not a base: %s" %ref
    assert val in ['A','T','G','C','.'], "val is not an acceptable value: %s" %val #do I need a test for  
    if val=='.':
      return('?')
    else:
      return(val)
#    if ref==val:
#        return('0')
#    else:
#        return('1')
  def export_fphase(self,outfile):
    import math
    fi=open(outfile,'w')
    fi.write(str(int(math.ceil(self.numinds/2.0)))+'\n')
    fi.write(str(self.nsnp)+'\n')
 #   fi.write(snp_locs)
    for ind in range(len(self.genos)):
#       fi.write("# id %s \n" %ind)
       for snp in self.snplist:
          fi.write(self.genos[ind][snp])
       fi.write('\n')
    if self.numinds%2!=0:
      for snp in self.snplist:
          fi.write('?')
      fi.write('\n')
  def tripletest(self):
    for snp in self.snplist:
      snpset=set()
      for ind in self.indivs:
        snpset.add(self.genos[self.indivs[ind]][snp])
    snpset.discard('?')
    if len(snpset) >2:
        print("Snp %i has more than 2 nucleotides" %snp)   
  def import_fphase(self,haplofile):
       self.impgenos=dict.fromkeys(range(len(self.genos)),None)
       for item in self.impgenos: 
          self.impgenos[item]={}
       fi=open(haplofile)
       for i,lin in enumerate(fi):
          if lin.startswith('#'):
            pass
          else:
            assert len(lin.strip())==len(self.snplist)
            indiv=i/2
            if indiv == self.numinds:
              break
            for snp, nuc in enumerate(lin.strip()):
                 if self.genos[indiv][self.snplist[snp]] in [nuc,'?']:
                      self.impgenos[indiv][self.snplist[snp]]=nuc
                 else:
                      print("%i indiv %s snp %s nuc %s raw %s"  %(i, indiv, snp, nuc, self.genos[indiv][self.snplist[snp]]))
                      print(lin)
  def export_chromop(self,outfile):
    import math
    outfi=open(outfile,'w')
    outfi.write('0\n')
    outfi.write(str(self.numinds)+'\n')
    outfi.write(str(len(self.snplist))+'\n')
    outfi.write('P')
    for item in self.snplist:
       outfi.write(' '+str(item))
    outfi.write('\n')
    outfi.write('S'*len(self.snplist)+'\n')
    for ind in self.indivs:
      for snp in self.impgenos[ind]:
        outfi.write(self.impgenos[ind][snp])
      outfi.write('\n')
    outfi.close()
  def import_copyprobs(self,infile):
    try: 
      call(["gunzip",infile])
    except: 
      pass
    self.copyprob_dict={}
    fi=open(infile[:-3])
    for i,lin in enumerate(fi):
      if i==0: pass
      elif lin.startswith('HAP'):
         hap=int(lin.split()[-1])-1
         self.copyprob_dict[hap]=[]
      else:
         lii=lin.split()
         indice=lii.index(max(lii[1:]))
         self.copyprob_dict[hap].append((lii[0],indice))
  def colors(self,val):
      from pylab import get_cmap
      NUM_COLORS = 30
      cm = get_cmap('Accent')
      color = cm(1.*val/NUM_COLORS)  # color will now be an RGBA tuple
      return(color)
  def painter_prep(self):
    self.linesections={}
    for hap in range(0,self.numinds):
        self.linesections[hap]=[]
        y=[int(hap),int(hap)]
        dat=self.copyprob_dict[hap]
        col=self.colors(dat[0][1]) #gets first hap
        loc=int(self.snplist[-1])
        count=0
        for indix,tup in enumerate(dat):
           oldcol=col
           col=self.colors(tup[1])
           count+=1
           if col!= oldcol:
              x=[loc,int(self.snplist[-indix])]
              self.linesections[hap].append([x,y,oldcol,count])
              loc=int(self.snplist[-indix])
              count=0
        x=[loc,int(self.snplist[0])]
        self.linesections[hap].append([x,y,oldcol,count])
  def paint(self,outfile):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_ylim([0,self.numinds+1])
        ax.set_xlim([-1500000,int(self.snplist[-1])])
        for hap in range(0,self.numinds):
            for item in self.linesections[hap]:
              if item[-1]>20:
                     ax.plot(item[0],item[1],color=item[2], linewidth=4)
            ax.plot([-1500000,0],item[1],color=self.colors(hap), linewidth=4)
            ax.text(-1500000,hap,self.indivs[hap])
        fig.savefig("%sh%i.pdf"%(outfile,hap))





def runner(infile,prefix):
   phasefi="%s.inp" %prefix
   c=Chromosome(infile)
   c.export_fphase(phasefi)
   call(["./fastPHASE_MacOSX-Darwin", "-n", "-B", "-T10", "-o%s" %prefix, phasefi])
   c.import_fphase('%s_haplotypes.out' %prefix)
   c.export_chromop("%s.chrominp" %prefix)
   call(["perl", "makeuniformrecfile.pl", "%s.chrominp" %prefix, "%s.recombfile" %prefix])
#   call(["./chromopainter-0.0.4/chromopainter", "-g", "%s.chrominp" %prefix, "-r", "%s.recombfile" %prefix, "-a", "0", "0", "-in"])
   call(["./chromopainter-0.0.4/chromopainter", "-g", "%s.chrominp" %prefix, "-r", "%s.recombfile" %prefix, "-a", "0", "0", "-in","-iM","-i","10", "-j", "-b"])
   c.import_copyprobs("%s.chrominp.copyprobsperlocus.out.gz" %prefix)
   c.painter_prep()
   c.paint(prefix)
   return(c)
  

c=runner("demo_data.txt", "demo_mini")


b=Chromosome('new_data.txt')
prefix="chrm1"
b.import_copyprobs("%s.chrominp.copyprobsperlocus.out.gz" %prefix)
b.painter_prep()
b.paint(prefix)


#b=Chromosome('new_data.txt')
#b.export_fphase('bigtest.inp')
#a.export_fphase("test.inp")
#a.tripletest()




c=runner("demo_data.txt", "demo_mini")



c=Chromosome("")
c.import_copyprobs("%s.chrominp.copyprobsperlocus.out.gz" %prefix)
c.painter("demo2")






copyprob_dict={}
for chrm in range(1,5):
   copyprob_dict[chrm]={}
   fi=gzip.open("Gdonor_chrm%i.copyprobsperlocus.out.gz"%chrm).readlines()
   i=1
   hap=1
   while i< len(fi):
    lin=fi[i]
    if lin.startswith('HAP'):
        hap=lin.split()[-1]
        i+=1
        copyprob_dict[chrm][hap]=[]
    else:
        copyprob_dict[chrm][hap].append(lin)
        i+=1
        yco=0


linesections={}
for chrm in copyprob_dict:
   linesections[chrm]={}

for chrm in [3,4]:
  for hap in range(1,1427):
    linesections[chrm][hap]=[]
    if hap%2==0:
        yco=hap-0.25
    else:
      yco=hap+0.25
    y=[yco,yco]
    dat=copyprob_dict[chrm][str(hap)]
    col=colors(dat[0].split()[1])
    i=dat[0].split()[0]
    count=0
    for lin in dat:
      lii=lin.split()
      oldcol=col
      col=colors(lii[1])
      count+=1
      if col!= oldcol:
         x=[i,lii[0]]
         linesections[chrm][hap].append([x,y,oldcol,count])
         i=lii[0]
         count=0
    x=[i,dat[-1].split()[0]]
    linesections[chrm][hap].append([x,y,oldcol,count])

def painter(a,b,chrm,outfile):
    fig = plt.figure()
    pylab.ylim([a,b])
    for hap in range(a,b):
            for item in linesections[chrm][hap]:
                plt.plot(item[0],item[1],item[2][0], linewidth=1.5)
    plt.savefig("../%s"%outfile)

