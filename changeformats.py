#import sys
#from __future__ import print_function
class Chromosome:
  def __init__(self,infii):
       self.infile=infii
       fi=open(infii)
       header=fi.readline().split()
       self.numinds=len(header[3:])
       self.indivs={index:elem for index, elem in enumerate(header[3:])}
       self.genos=dict.fromkeys(header[3::],None)
       self.multihit=[]
       for item in header[3:]: 
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
                self.genos[self.indivs[ii]][snp]=self.translate(ref,item.split('/')[0])
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
    for ind in self.genos:
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
       self.impgenos=dict.fromkeys(header[3::],None)
       for item in header[3:]: 
          self.genos[item]={}
       fi=open(haplofile)
       for lin in fi:
          if lin.startswith('#'):
            pass
          else:
            pass

a=Chromosome('testdata.txt')



b=Chromosome('new_data.txt')

b.export_fphase('bigtest.inp')

#a.export_fphase("test.inp")

#a.tripletest()
