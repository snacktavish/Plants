import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import gzip
from subprocess import call
class Chromosome:
  def __init__(self,infii,prefix):
       self.prefix=prefix
       self.infile=infii
       fi=open(infii)
       header=fi.readline().split()
       self.numinds=len(header[3:])
       self.indivs={index:elem for index, elem in enumerate(header[3:])}
       self.rev_indivs={elem:index for index, elem in enumerate(header[3:])}
       self.indlist=self.rev_indivs.keys()
       self.genos=dict.fromkeys(range(len(header[3:])),None)
       self.multihit=[]
       for item in self.genos: 
          self.genos[item]={}
       self.snplist=[]
       self.ref_dict={}
       for num,lii in enumerate(fi):
          lin=lii.strip().split('\t')
          snp=lin[1]
          if num==0:
                self.chrm=lin[0].split('_')[-1]
          if len({x.strip('/') for x in lin[2:] if x not in ['./','./.']})==2:
            self.snplist.append(snp)
            ref=lin[2]
            self.ref_dict[snp]=ref
            for ii,item in enumerate(lin[3:]):
                self.genos[ii][snp]=self.translate(ref,item.split('/')[0])
          if len({x.strip('/') for x in lin[2:] if x != './'})>2:
#            print("Snp %s has more than 2 nucleotides" %snp) 
            self.multihit.append(snp)  
          if num%1000000==0:
            print(num)
       self.tripletest()
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
    fi.write(str(int(math.ceil(len(self.indivs)/2.0)))+'\n')
    fi.write(str(len(self.snplist))+'\n')
 #   fi.write(snp_locs)
    for ind in self.indorder:
#       fi.write("# id %s \n" %ind)
       for snp in self.snplist:
          fi.write(self.genos[ind][snp])
       fi.write('\n')
    if len(self.indivs)%2!=0:
      for snp in self.snplist:
          fi.write('?')
      fi.write('\n')
  def tripletest(self):
    self.triples=[]
    for snp in self.snplist:
      snpset=set()
      for ind in self.indivs:
        snpset.add(self.genos[ind][snp])
    snpset.discard('?')
    if len(snpset) >2:
        self.triples.append(snp)
        print("Snp %i has more than 2 nucleotides" %snp)
    for ind in self.indivs:
        for snp in self.triples:
          if snp in genos[ind]:
              del self.genos[ind][snp]
          snplist.remove(snp)
  def import_fphase(self,haplofile):
       self.impgenos=dict.fromkeys(self.genos.keys(),None)
       for item in self.impgenos: 
          self.impgenos[item]={}
       fi=open(haplofile).readlines()
       imp=[lin for lin in fi if lin!= '# \n']
       for i,ind in enumerate(self.indorder):
            lin=imp[i]          
            assert len(lin.strip())==len(self.snplist)
            for snp, nuc in enumerate(lin.strip()):
                 if self.genos[ind][self.snplist[snp]] in [nuc,'?']:
                      self.impgenos[ind][self.snplist[snp]]=nuc
                 else:
                      print("%i indiv %s snp %s nuc %s raw %s"  %(i, ind, snp, nuc, self.genos[indiv][self.snplist[snp]]))
                      print(lin)

  def export_nexus(self,filename,imp=True):
      assert(self.impgenos)
      nexfi=open(filename,'w')
      nexfi.write("#nexus\nbegin data;\n")
      nexfi.write("dimensions ntax={ntax} nchar={nchar};\n".format(ntax=self.numinds,nchar=len(self.snplist)))
      nexfi.write("format datatype=dna interleave=no gap=- missing=?;\n")
      nexfi.write("Matrix\n")
      if imp==True:
        dat=self.impgenos
      else:
       dat=self.genos
      for item in dat:
        nexfi.write("{0: <16}".format(self.indivs[item]))
        for snp in self.snplist:
          nexfi.write(dat[item][snp])
        nexfi.write('\n')
      nexfi.write(';\nend;\n')
  def trans(self,snp,val): #Same as ref = 0, alt=2
        assert(set([self.ref_dict[snp],val]).issubset(set(['A','T','G','C'])))
        if self.ref_dict[snp]==val:
                   return('0')
        else:
                    return('2')
  def export_EIG(self,outstr,group='None'):
        snpfi=open("eig/"+outstr+'.snp','w')
        goufi=open("eig/"+outstr+'.geno','w')
        indfi=open("eig/"+outstr+'.ind','w')
        for ind in self.impgenos:
          for snp in self.impgenos[ind]:
                gen=self.trans(snp,self.impgenos[ind][snp])
                goufi.write(snp+' '+'ind'+str(ind)+' '+gen+'\n')
        goufi.close()
        for ind in self.indivs:
            if group=='None':
              indfi.write('ind'+str(ind)+'\tF\t'+self.indivs[ind]+'\n')#by region      
        indfi.close()
        for snp in self.snplist:
            chrm=self.chrm
            snpfi.write(" ".join([snp,str(chrm),'0.0',snp])+'\n') #pulls outsnp number, chrm, dummy recombination distance, lcoaction
        snpfi.close()
        par=open("eig/par.example",'r')
        opar=open('eig/par.'+outstr,'w')
        for lin in par.readlines():
             opar.write(lin.replace('example',outstr))
        opar.close()
  def export_chromop(self,prefix,donor=0):
    import math
    outfi=open(prefix+chrominp,'w')
    outfi.write('0\n')
    outfi.write(str(len(self.indivs)+'\n')
    outfi.write(str(len(self.snplist))+'\n')
    outfi.write('P')
    if donor:
      assert type(donor)==dict
      dfi=open(prefix+"dlist",'w')
      new_indlist=[]
      for pop in donor:
          dfi.write(" ".join([pop,len(donor[pop])]))
          for ind in pop:
            new_indlist.append(ind)
      for item in self.indlist:
          if item not in new_indlist:
            new_indlist.append(item)
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
        line=dat[0][1] #gets first hap
        loc=int(self.snplist[-1])
        count=0
        for indix,tup in enumerate(dat):
           oldline=line
           line=tup[1]
           count+=1
           if line!= oldline:
              x=[loc,int(self.snplist[-indix])]
              self.linesections[hap].append([x,y,oldline,count])
              loc=int(self.snplist[-indix])
              count=0
        x=[loc,int(self.snplist[0])]
        self.linesections[hap].append([x,y,oldline,count])
  def paint(self,outfile):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_ylim([0,self.numinds+1])
        ax.set_xlim([-1500000,int(self.snplist[-1])])
        for hap in range(0,self.numinds):
            for item in self.linesections[hap]:
              if item[-1]>20:
                     ax.plot(item[0],item[1],color=self.colors(item[2]), linewidth=4)
            ax.plot([-1500000,0],item[1],color=self.colors(hap), linewidth=4)
            ax.text(-1500000,hap,self.indivs[hap])
        fig.savefig("%sh%i.pdf"%(outfile,hap))
        print("DOUBLE CHECK THESE COLOU?RS")







def strip_missing(self,snpcutoff1,indcutoff,snpcutoff2):
    snpcount={snp:0 for snp in self.snplist}
    for ind in self.genos:
         for snp in self.genos[ind]:
          if self.genos[ind][snp]!='?':
           snpcount[snp]+=1
    snpdel={snp for snp in snpcount if snpcount[snp]<snpcutoff1*len(self.indivs)}
    pruned_genos={ind:{} for ind in self.genos}
    for ind in self.genos:
      print(ind)
      for snp in self.genos[ind]:
        if snp not in snpdel:
          pruned_genos[ind][snp]=self.genos[ind][snp]
    indcount={ind:0 for ind in self.genos}
    for ind in pruned_genos:
      for snp in pruned_genos[ind]:
        if pruned_genos[ind][snp]!='?':
         indcount[ind]+=1
    inddel={ind for ind in indcount if indcount[ind]<indcutoff*(len(self.snplist)-len(snpdel))}
    print("inddel is",inddel)
    for ind in inddel:
      del pruned_genos[ind]
      del self.rev_indivs[self.indivs[ind]]
      del self.indivs[ind]
    snpcount2={snp:0 for snp in set(self.snplist)-set(snpdel)}
    for ind in pruned_genos:
      for snp in pruned_genos[ind]:
        if pruned_genos[ind][snp]!='?':
         snpcount2[snp]+=1
    snpdel2={snp for snp in snpcount2 if snpcount2[snp]<snpcutoff2*len(self.indivs)}
    pruned_genos2={ind:{} for ind in pruned_genos}
    for ind in pruned_genos:
      print(ind)
      for snp in pruned_genos[ind]:
        if snp not in snpdel2:
          pruned_genos2[ind][snp]=self.genos[ind][snp]
    self.genos=pruned_genos2
    self.snplist=list(set(self.snplist)-set(snpdel)-set(snpdel2))
    self.numinds=len(pruned_genos2)
    self.indorder=self.indivs.keys()
    self.indorder.sort()







def  parse_partitions(chrom, partfile,outfile):
     fi=open(partfile).readlines()
     for i,lin in enumerate(fi):
       if lin.startswith("Ntax"):
         end=lin.split("Nchar=")[-1].split()[0]
       if lin.startswith("MDLscore"):
         num=i+1
     dat=fi[num].split()
     dat.append(end)
     breaks=[]
     for item in dat[2:]:
         breaks.append(chrom.snplist[int(item)-1])
     ofi=open(outfile,'w')
     ofi.write('\n'.join(breaks))
     ofi.close()




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



markerslis= ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','1','2','3','4','8','^',r'$\alpha$',r'$\beta$',r'$\clubsuit$',r'$\spadesuit$',r'$\star$']
colorslis=['red','blue','green','black']

for item in 

import matplotlib.pyplot as plt

A=["1054.bam",]
B=["SRR072711.bam", "SRR072712.bam"]


def slider(chrom, prefix, window):
  chrom.export_EIG(prefix)
  markerslis= [r'$\spadesuit$',r'$\star$', r'$\spadesuit$', r'$\clubsuit$','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','1','2','3','4','8','^']
  singmarks=[r'$\gamma$', r'$\sigma$',r'$\infty$',r'$\alpha$',r'$\beta$',r'$\Theta$', r'$\Xi$', r'$\Phi$',r'$\$$', r'$\#$', r'$\%$', r'$\S$',r'$A$',r'$B$',r'$C$',r'$D$',r'$E$',r'$F$']
  colorslis=['red','green','black','orange','pink']
  symdict={}
  for ind,sym in zip(chrom.rev_indivs,zip(markerslis+singmarks+markerslis,colorslis*15)):
     symdict[ind]=sym
  start=0
  end=0
  n=(len(chrom.snplist)/window)+2
  fulldict={}
  for item in chrom.rev_indivs:
    fulldict[item]=[]
  fulldict["ends"]=[]
  fig=plt.figure(figsize=(5,n*5))
  i=0
  while end < len(chrom.snplist):
      print(end < len(chrom.snplist))
      i+=1
      print("end is")
      print(end)
      fi=open("eig/%s.snp"%prefix).readlines()
      end=start+window
      del fi[start:end]
      start=end
      try:
        endloc=chrom.snplist[end]
      except:
        endloc=chrom.snplist[-1]
      bad=open("eig/bad%s.snp"%prefix,"w")
      for lin in fi:
        bad.write(lin)
      bad.close()
      call(["../EIG5.0.1/bin/smartpca", "-p", "eig/par.%s"%prefix])
      outdict={}
      outs=open("eig/%s.evec"%prefix).readlines()
      for lin in outs[1:]:
        outdict[lin.split()[-1]]=(float(lin.split()[1]),float(lin.split()[2]))
      ax=fig.add_subplot(n,1,i+1)
      ax.set_title=(str(start))
      ax.set_xlim(-1,1)
      ax.set_ylim(-1,1)
      ax.text(-0.5,0.45,endloc)
      for item in outdict:
        x,y=outdict[item]
        ax.plot(x,y,marker=symdict[item][0],color=symdict[item][1],alpha=0.6)
  keyz=fig.add_subplot(n,1,1)
  keyz.set_xlim(0.5,10)
  keyz.set_ylim(-5,5+(len(outdict.keys()))/3)
  for e,item in enumerate(outdict):
      if e<=(len(outdict.keys())/3):
            keyz.plot(1,e,marker=symdict[item][0],color=symdict[item][1])
            keyz.text(1.5,e,item, fontsize=6)
      elif (len(outdict.keys())/3)<e<=2*(len(outdict.keys())/3):
            keyz.plot(3,0.5+(e-(len(outdict.keys())/3)),marker=symdict[item][0],color=symdict[item][1])
            keyz.text(3.5,0.5+(e-(len(outdict.keys())/3)),item, fontsize=6)
      else:
            keyz.plot(6,0.5+(e-(2*(len(outdict.keys())/3))),marker=symdict[item][0],color=symdict[item][1])
            keyz.text(6.5,0.5+(e-(2*(len(outdict.keys())/3))),item, fontsize=6)
  fig.tight_layout()
  fig.savefig("%s.pdf"%prefix)


import copy

def subset_inds(chrom,keeplis):
    x=copy.deepcopy(chrom)
    for item in x.rev_indivs.keys():
      if item not in keeplis:
        indnum=x.rev_indivs[item]
        del x.genos[indnum]
        del x.impgenos[indnum]
        del x.indivs[indnum]
        del x.rev_indivs[item]
    return(x)


'''      avals=[]
      bvals=[]
      for item in A:
        avals.append(outdict[item])
      for item in B:
        bvals.append(outdict[item])
      aavg=(sum(avals)/len(avals))
      bavg=(sum(bvals)/len(bvals))
      if aavg<bavg:
        aavg=aavg*-1
        bavg=bavg*-1
      for item in outdict:
        if outdict[item]<bavg or outdict[item]>aavg:
             fulldict[item].append("n/a")
        else:
          fulldict[item].append((outdict[item]-bavg)/(aavg-bavg))
      fulldict["ends"].append(end)
      return(fulldict)
'''


def post_proc(c,prefix):
  phasefi="%s.inp" %prefix
  export_fphase(c,phasefi)
  call(["./fastPHASE_MacOSX-Darwin", "-n", "-B", "-T10", "-o%s" %prefix, phasefi])
  import_fphase(c,'%s_haplotypes.out' %prefix)
  slider(c,prefix,len(c.snplist)+5)


c=Chromosome("big_demo.txt")
c.import_fphase('%s_haplotypes.out' %prefix)
c.export_nexus('test.nex')

c.export_EIG(prefix)


def runner(infile,prefix):
   phasefi="%s.inp" %prefix
   c=Chromosome(infile)
   c.export_fphase(phasefi)
   call(["./fastPHASE_MacOSX-Darwin", "-n", "-B", "-T10", "-o%s" %prefix, phasefi])
   c.import_fphase('%s_haplotypes.out' %prefix)
   c.export_chromop(prefix)
   call(["perl", "makeuniformrecfile.pl", "%s.chrominp" %prefix, "%s.recombfile" %prefix])
   call(["./chromopainter-0.0.4/chromopainter", "-g", "%s.chrominp" %prefix, "-r", "%s.recombfile" %prefix, "-a", "0", "0", "-in","-iM","-i","10", "-j", "-b"])
   c.import_copyprobs("%s.chrominp.copyprobsperlocus.out.gz" %prefix)
   c.painter_prep()
   c.paint(prefix)
   return(c)

def eigrun(chrom, prefix):
   chrom.export_EIG(prefix)
   call(["../EIG5.0.1/bin/smartpca", "-p", "eig/par.%s"%prefix])
   call(["perl", "../EIG5.0.1/bin/ploteig", "-i", "eig/%s.evec"%prefix, "-p", "eig/namesoutg.txt", "-x", "-k", "-o", "eig/%s.xtxt"%prefix])


 
    
c=runner("big_demo.txt", "demo_big")

c.export_EIG(prefix)

c=Chromosome("big_demo.txt")
  
b=Chromosome('new_data.txt')
prefix="chrm1"
b.import_fphase('%s_haplotypes.out' %prefix)
b.export_nexus("chrm1.nex")

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

