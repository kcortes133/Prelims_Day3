# Characterization of Viral Composition in Metagenomics Sample

## Goal
Given a metagenomics sample give information about the composition of the viral community.

## Description
This program takes fastq files as input and outputs information about the viral community in the sample. 

Due to the large size of metagenomic samples, reads are limited to the first 150 bases and only 1/4 of randomly
subsetted data was used. The reads are then selected for likelihood of viral origin using kmer matching.Kmer
matching uses subsequences of the reads, of length k, and compares them to the subsequences of the viral genomes.
Reads with a high number of matching kmers are then aligned to the viral genomes. The viral genome with the highest score
during alignment is what the read is labeled as. 

Information about the viral population is then output to a csv file. Pie charts can be made to look at the relative 
composition of the sample, the viruses and virus host types represented. 

Example files are a subset of the SRR12464727 run from NCBI. 

## Install
- numpy~=1.22.4
- matplotlib~=3.5.2
- pandas~=1.4.2

## Usage
### Python
```python 
def main():
    readsF = 'Example/Input/SRR12464727.fastq'
    run = 'SRR12464727'
    
    kLen = 10

    # unzip viral genome files
    # will do automatically if genomes file doesnt exist and
    # 14,858 viral genomes
    # unzip the viral genomes in compressed folder if not unzipped
    if not os.path.exists('Outputs/'):
        os.makedirs('Outputs')

    # unzip viral genome databaset
    if not os.path.exists('viralDB/genomes'):
        os.makedirs('viralDB/genomes')
        vDir = 'viralDB/compressed/viral.'
        uvDir = 'viralDB/genomes/viral.'
        for i in range(1,5):
            vfile = vDir + str(i) + '.1.genomic .fna.gz'
            ovFile = uvDir + str(i) + '.genomic.fna'
            database.extract(vfile, ovFile)
            
    # make dictionary of viral database
    vGDB = database.makeVDB()
    s = time.time()
    # get kmers from viral genomes
    kmers = database.getkmers(vGDB, kLen)
    virusOutFile = 'Outputs/virusCount' + run +'_all.csv'

    # get reads from fastq file
    reads = dataExploration.getReads(readsFile)
    # randomly subsample reads
    reads = random.sample(reads, int(len(reads)/4))
    print(time.time() - s)
    # kmer matching
    top, totHits = kmerBinning.subsetReads(reads, kLen, kmers)
    print('Time to subset reads: ',time.time()-s)
    # align reads to viral genomes
    virusesFound, viralNames = kmerBinning.getVirusMatches(top, totHits, reads, vGDB)

    print('Possible Viruses found: ', len(top))
    print('Total time:  ', time.time()-s)
    viralHostsDB = database.getViralHosts()

    outFiles.writeViruses(virusOutFile, virusesFound, viralHostsDB, viralNames)
    plotGenerator.pieFromFile(virusOutFile, 10)
    plotGenerator.pieFromFile_Hosts(virusOutFile, 10)


main()
```
### Command Line
unzip the viral genome database
```
$ python .\main.py 

get viral community information from example dataset
$ python .\main.py --kBinning=True --readsFile='Example/Input/SRR12464727_example.fastq' --run='SRR12464727_example' 

get visualizations for example dataset
$ python .\main.py --visualize=True --viralSampleFile='Example/Output/virusCountSRR12464727_example_all.csv'         
```

## Input
- metagenomic fastq files

## Output
- csv file with viral community information
- virusID,virusName,virusCount,virusHost
- pie chart of 9 viruses and hosts with highest count in sample

Example pie charts for all of SRR12464727 reads, Example folder has results for subset of this data
![hosts](https://user-images.githubusercontent.com/22487858/173917029-1eac17ff-42cb-49eb-a074-354481afd36d.png)

![SRR12464727 pie chart](https://user-images.githubusercontent.com/22487858/173917086-afc0039c-347b-48c7-8fd5-3b6dea2a3fd1.png)

