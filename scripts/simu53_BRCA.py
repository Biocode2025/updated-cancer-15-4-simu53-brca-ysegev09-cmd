import random
def  Comp_seq(old,new):
  lenth_o=len(old)
  lenth_n=len(new)
  difranz=0
  i=0
  if lenth_o==lenth_n:
    while i<lenth_n:
      if old[i]!=new[i]:
        difranz=difranz+1
      i = i + 1  
  elif lenth_o>lenth_n:
    while i<lenth_n:
      if old[i]!=new[i]:
        difranz=difranz+1 
      difranz=lenth_o-lenth_n
      i = i + 1
  elif lenth_o<lenth_n:
    while i<lenth_o:
      if old[i]!=new[i]:
        difranz=difranz+1
      difranz=lenth_n-lenth_o
      i = i + 1
    

  return difranz
    
    
def DNA_RNA_Cod(line) : # התרגום לRNAוהגדלה
  big_leter=line.upper()
  hachlafa=big_leter.replace("T","U")
  return(hachlafa)

def Read_dict() :
  global RNA_codon_table

  RNA_codon_table={}
  milon = open('data/codon_AA.txt', 'r')
  for line in (milon):
    line = line.rstrip("\r\n")
    num=line[0:3]
    RNA_codon_table[num] = line[-1] 

def RNA_prot (seq):
  pro=""
  t=3
  i=0
  while t <= len(seq):
      
    sd=seq[i:t]
    pro=pro+RNA_codon_table.get(sd,"not found")
    
    i = i + 3
    t = t + 3
    if pro[-1]=="*":
      break
  
  return pro


def Mutate_DNA(seq):
  list=["A", "C", "G", "U"]
  num=len(seq)
  
  m=random.randrange(0,num)
  n=random.choice(list)
  seq=seq[0:m]+n+seq[(m+1):]
  return(seq)

def  Delete_DNA(seq):
  cama=random.randrange(0,3)
  for i in range(cama):
    num=len(seq)
    m=random.randrange(0,num)
    seq=seq[0:m]+seq[(m+1):]
  return(seq)

def Insert_DNA(seq):
  list=["A", "C", "G", "U"]
  num=len(seq)
  cama=random.randrange(0,3)
  for i in range(cama):
    m=random.randrange(0,num)
    n=random.choice(list)
    seq=seq[0:m+1]+n+seq[(m+1):]
  return(seq)
  
has_BRCA = input("Does the Female has a BRCA1,2 mutation? (Y=Yes, N=No) ")
if has_BRCA == 'Y':
    mutations_needed = 1
else:
    mutations_needed = 2
dop_until_chn=0
seq=""  
total_replications = 0
rep=0
file = open('data/pck_coli.txt', 'r')
make_mil=Read_dict()
for line in(file):
  if line[0]!=">":
    line = line.rstrip("\r\n")
    call3=DNA_RNA_Cod(line)
    seq=seq+call3
mata=seq
gan=0
avrag_run=0
avrag_run_c=0

for avrag_run_c in range(1000):
  
  while True: 
    avrag_run=avrag_run+1
    chan_mata=random.randrange(1,10001)
    if chan_mata==5000:
      ran_mata=random.randrange(1,101)
      if ran_mata<98: 
        
        dop_until_chn=dop_until_chn+1
        mata=Mutate_DNA(mata)
        protin=RNA_prot(seq) 
        protin_mata=RNA_prot(mata)
        ASVAHA=Comp_seq(protin,protin_mata)
        if ASVAHA==mutations_needed:
          break
      elif ran_mata==99:
        
        dop_until_chn=dop_until_chn+1
        mata=Insert_DNA(mata)
        break
      else:
        
        dop_until_chn=dop_until_chn+1
        mata=Delete_DNA(mata)
        break
  
  
avrag_run=(avrag_run/10)
avg_years=(avrag_run)/365

protin=RNA_prot(seq) 
protin_mata=RNA_prot(mata)
ASVAHA=Comp_seq(protin,protin_mata)
file_a=open('results/results_can.fasta', 'w')
file_a.write("original: "+protin+ "\n")
file_a.write("mutation: "+protin_mata+"\n" )
file_a.write("the num of replications is:   "+ str(dop_until_chn) +"\n")
file_a.write("The The average number of years required:   "+str(avg_years)+"\n")
file_a.close()
file_a=open('results/results_can.fasta', 'r')

for line in (file_a):
  print(line)
  