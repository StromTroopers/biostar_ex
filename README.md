# biostar_ex
```
import pandas as pd

tab=pd.read_csv("/beegfs/data/user/these/Transciptomic_analysis2/Table_RNA_reads3.txt",sep=";")


Species_name_file="/beegfs/data/user/these/Species_genome_names.txt"
list_of_names1=[]
for names in open(Species_name_file,"r"):
 	list_of_names1.append(names)

list_of_names2=[]
for names in list_of_names1:
 	list_of_names2.append(names.replace("\n", ""))

list_of_names2=['Anoplophora_glabripennis']

import os
assembly="Assembly1"
def lookupSRR(assembly):
  list_of_SRRs=[]
  subtab=tab.loc[tab["ScientificName"].str.contains(assembly)]
  for SRR in subtab["Run"]:
    list_of_SRRs.append(SRR)
  return list_of_SRRs

import os
def getSRRxargs(assembly):
  list_of_SRRs_1=[]
  list_of_SRRs_2=[]
  list_of_SRRs_single=[]
  for file in os.listdir("/Users/user/Desktop/these/snakemake_test_mapping/"+assembly+"/Mapping/"):
            if file.endswith("_1.tgz"):
              list_of_SRRs_1.append(file)
            if file.endswith("_2.tgz"):
              list_of_SRRs_2.append(file)
            elif file.endswith(".tgz"):
              if '_1' in file :
                continue
              else :
                list_of_SRRs_single.append(file)
  argument_single="-R ",",".join(list_of_SRRs_single)," "
  argument_paired="-1 ",",".join(list_of_SRRs_1)," -2 ",",".join(list_of_SRRs_2)
  if len(list_of_SRRs_single) >= 1 and len(list_of_SRRs_1) >= 1:
    arguments="".join(argument_single+argument_paired)
    return arguments
  if len(list_of_SRRs_single) >= 1 and len(list_of_SRRs_1) == 0:
    arguments=argument_single
    return "".join(argumentq)
  if len(list_of_SRRs_single) == 0 and len(list_of_SRRs_1) >= 1:
    arguments=argument_paired
    return "".join(arguments)

def Assembly_to_choose(species):
  if os.path.exists("/beegfs/data/user/these/Genomes/"+species+"/"+species+"_corrected2.fa"):
    assembly_file="/beegfs/data/user/these/Genomes/"+str(species)+"/"+str(species)+"_corrected2.fa"
  elif os.path.exists("/beegfs/data/user/these/Genomes/"+species+"/"+species+"_corrected.fa"):
    assembly_file="/beegfs/data/user/these/Genomes/"+str(species)+"/"+str(species)+"_corrected.fa"
  elif os.path.exists("/beegfs/data/user/these/Genomes/"+species+"/"+species+"_bis.fa"):
    assembly_file="/beegfs/data/user/these/Genomes/"+str(species)+"/"+str(species)+"_bis.fa"
  return assembly_file


#Write recipes that are able to render those targets


rule all:
     input:
      expand("/beegfs/data/user/these/Genomes/{species}/Mapping/Download_run.log", species = list_of_names2),
      expand("/beegfs/data/user/these/Genomes/Transcriptomic_reads/Reads_{read}", read = lookupSRR(config["list_of_names2"]))
      #expand("TEMPDIR/{setA}_{setB}_", setA = config["setA"], setB = config["setB"])

rule Download_reads:
     output: 
      outfile="/beegfs/data/user/these/Genomes/{species}/Mapping/Download_run.log",
      outfile2=expand("/beegfs/data/user/these/Genomes/Transcriptomic_reads/Reads_{read}", read = lookupSRR('{species}'.format(species={species})))
     shell:
      """
      rm -rf /beegfs/data/user/sra_reads/sra/{wildcards.read}*
      /beegfs/data/user/TOOLS/sratoolkit.2.11.0-ubuntu64/bin/prefetch --max-size 100000000 {wildcards.read} && /beegfs/data/user/TOOLS/sratoolkit.2.11.0-ubuntu64/bin/fasterq-dump -t /beegfs/data/user/these/Genomes/{wildcards.species}/Mapping/Transcriptomic_reads/ --threads 10 -f {wildcards.read} -O /beegfs/data/user/these/Genomes/{wildcards.species}/Mapping/Transcriptomic_reads/
      /beegfs/data/user/user_conda/bin/pigz --best /beegfs/data/user/these/Genomes/{wildcards.species}/Mapping/Transcriptomic_reads/{wildcards.read}*\n")
      touch {output.outfile}
      touch {output.outfile2}
      rm -rf /beegfs/data/user/sra_reads/sra/{params.read}*
      """
```




