class RawData:
  def __init__(self,infii):
       self.infile=infii
       fi=open(infii)
       header=fi.readline().split('\t')
       indivs={index:elem for index, elem in enumerate(header[3:])}
       self.genos=dict.fromkeys(header[3::],None)
       for item in header[3:]:
          genos[item]={}
       i=1
       for lii in fi:
          lin=lii.split('\t')
          if i==1:
              self.chrm=lin[0]
          snp=lin[1]
          ref=lin[2]
          for ii,item in enumerate(lin[3:]):
                self.genos[indivs[ii]][snp]=self.translate(ref,item.split('/')[0])
          i+=1
       self.nsnp=i
  def translate(self,ref,val): #this is setting allele in refernce as 0
    assert ref in ['A','T','G','C'], "ref is not a base: %s" %ref
    assert val in ['A','T','G','C','.'], "val is not an acceptable value: %s" %val #do I need a test for  
    if val=='.':
      return('?')
    if ref==val:
        return('0')
    else:
        return('1')
  def export_fphase(outfile):
    #
    pass






a=RawData('testdata.txt')





def gen(X,Y): #Decode's AB in genotypes
    if X=='A' and Y =='A': return '1'
    if X=='B' and Y =='B': return '2'
    if X=='A' and Y =='B': return '3'
    if X=='B' and Y =='A': return '3'


def ABtrans(A,B):
  if A and B in ['1','2']:
    if A==B=='1': return '1'
    if A==B=='2': return '2'
    if (A,B)==('1','2'): return '3'
    if (A,B)==('2','1'): return '3'
  else: return "error"


def phase_inp(geno,mark_dict,filtered_snps,outstr):
      for chrm in mark_dict:
        print("Chromosome %s" %chrm)
        oufi=open("%s%s.inp"%(outstr,chrm),'w')
        ii=0
        ids=[]
        for mark in mark_dict[chrm]:
              if mark[0] in filtered_snps:
                     ids.append(str(mark[1]))
                     ii+=1
        oufi.write(str(len(geno))+"\n")
        oufi.write(str(ii)+"\n"+"P ")
        oufi.write(' '.join(ids)+'\n')
        oufi.write("S" * ii +"\n")#WTF ARE THESE SSSSs for?
        for ind in geno:
                oufi.write("# id %s"%ind+"\n")
                hap1=[]
                hap2=[]
                for mark in mark_dict[chrm]:
                  if mark[0] in filtered_snps:
                    try:
                      genot=geno[ind][mark[0]]
                    except KeyError:
                      genot='?'
                    if genot=='1':
                       hap1.append('1')
                       hap2.append('1')
                    if genot=='2':
                       hap1.append('2')
                       hap2.append('2')
                    if genot=='3':
                       hap1.append('1')
                       hap2.append('2')
                    if genot=='?' or gen=='10':
                       hap1.append('?')
                       hap2.append('?')
                    if genot not in ['1','2','3','?','10']:
                       hap1.append('?')
                       hap2.append('?')
                assert(len(hap1)==len(hap2))
                oufi.write("".join(hap1)+'\n')
                oufi.write("".join(hap2)+'\n')
        oufi.close()





def phased_to_raw(phased_pref,outstr,chrms):
    for chrm in chrms:
        oufi=open('%s%s.txt'%(outstr,chrm),'w')
        snps=[mark_num_dict[x[0]] for x in mark_dict[chrm]] #gets the number for each snp from the name in mark dict
        hg=open("phase_outputSPAN/chrm%s_hapguess_switch.out"%(chrm)).readlines() #opens fastphase outputfile  
        assert len(snps)==len(hg[22].split())#makes sure number of snps in head == number in lines
        for x,lin in enumerate(hg):#skipps all the headers and footers
           if lin.startswith('# id'):
                idi=lin.split()[2]
                for i,snp in enumerate(snps):
                    A=hg[x+1].split()[i]
                    B=hg[x+2].split()[i]
                    geno=ABtrans(A,B)
                    oufi.write("["+snp+','+idi+','+geno+']'+'\n')
           else:pass       


def raw_to_eig(infiles,outstr,group='region'):#infiles is a list of files    
    try:
        len(mark_nam_dict)
    except:
      print("Create marker dictionary!")
    def trans(x):
        if x=='1':return '0'
        if x=='2':return '2'
        if x=='3':return '1'
        if x=='10':return '10'
    def ind_info(ind):
        ind_info=open('ids_fullSPAN.csv','r').readlines()
        inddict={}
        for lin in ind_info:
            lii=lin.split(',')
            inddict[lii[0]]=lii[1:]
        return inddict[ind]
    if type(infiles) != list:
       print("infiles needs to be alist")
    snps=set()
    ids=set()
    snpfi=open("eigSPAN/"+outstr+'.snp','w')
    goufi=open("eigSPAN/"+outstr+'.geno','w')
    indfi=open("eigSPAN/"+outstr+'.ind','w')
    for infile in infiles:
        goufi=open("eigSPAN/"+outstr+'.geno','a')
        indfi=open("eigSPAN/"+outstr+'.ind','a')
        print(infile)
        fi=open(infile)
        rawdat=[]
        for lin in fi:
          try:
                stri=str(str(lin)[1:-2]) #strips off brackets
                stri2=(str(stri)).split(',')
                snpnum=stri2[0].strip()
                snpnam=mark_nam_dict[snpnum]
                lst=[snpnum,stri2[1].strip(),trans(stri2[2].strip())]
                if stri2[1].strip() not in ind_skip:
                  ids.add(stri2[1].strip())
                  snps.add(snpnum)
                  if (lst[2]!='10'):
                                rawdat.append(lst)
          except:
                print(lin)
        for lin in rawdat:
                goufi.write(lin[0]+' '+lin[1]+' '+lin[2]+'\n')
        goufi.close()
    for item in ids:
        if group=='region':
               indfi.write(item+'\tF\t'+ind_info(item)[4].strip('"')+'\n')#by region      
        if group =='breed':
               indfi.write(item+'\tF\t'+ind_info(item)[2].strip('"')+'\n')#by breedcode
        if group =='country':
               indfi.write(item+'\tF\t'+ind_info(item)[5].strip('"')+'\n')
    indfi.close()
    for snp in snps:
            nam=mark_nam_dict[snp] 
            (num,chrm, pos)=mark_chrm_dict[nam]
            snpfi.write(" ".join([snp,chrm,'0.0',str(pos)])+'\n') #pulls outsnp number, chrm, dummy recombination distance, lcoaction
    snpfi.close()
    par=open("eig/par.example",'r')
    opar=open('eigSPAN/par.'+outstr,'w')
    for lin in par.readlines():
           opar.write(lin.replace('example',outstr))
    opar.close()
