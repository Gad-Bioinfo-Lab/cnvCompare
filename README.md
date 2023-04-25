# cnvCompare

cnvCompare is a bioinformatics tool used to compare and counts Copy Number Variations (Del or Dup) from several samples
The tool algorithm is described in : 


## Author
Yannis Duffourd (Inserm U1231 - Team GAD / CHU Dijon)

## Contributors
    Emilie Tisserant (Inserm U1231 - Team GAD)
    Anthony Auclair (Inserm U1231 - Team GAD / CHU Dijon)
    Marine Bergot (Inserm U1231 - Team GAD / CHU Dijon)
    Philippine Garret (Inserm U1231 - Team GAD)

## Requirements 
- cmake 3.10 or superior
- Boost 1.46.1 or superior
- gcc 11.x or superior (supporting ISO C++20)


## Install
git clone <pppp>
cd cnvCompare 
cmake . 
make 

## Usage 
cnvCompare -i </path/to/the/file/list> --vcf/--bed
    Allowed options : 
        -h, --help : displaying help 
        -i, --input <string>: List of input file(s) containing detected CNV from samples (prefer absolute path)
            These files will be all rewritten using the suffix provided including the counts for each DEL/DUP events
        -c, --control <string>: List of control file(s) containing detected CNV from samples (prefer absolute path)
            These files wont be rewritten but events contained will be counted
        -f, --filter <int>: apply a filter on the CNV size to be counted (default:0)
            Currently not yet implemented
        -w, --whole : Using the whole mode, counting on the whole genome. 
            Faster but need large amounts of RAM, depending on the size of the events, and the number of samples
            By default, if this option is not activated, the counts are realised chromosome by chromosome
        -d, --dict <string>: Path to the dictionnary file used to populate the chomosome list
            Mainly used for non-human genomes, if not provided will used chromosomes : [1-22],X,Y,MT/M
        --vcf : Input files are in VCF format 
        --bed : Input files are in BED format
        -s, --suffix <string>: The suffix used to name the ouput files
            The output files are written in the same path of the input file, adding this suffix in their names before the extension (.vcf or .bed or anyhting else) (default : "count")
        
You need to provide only 2 mandatory options : the input file list and the file format

examples :
    cnvCompare -i inputFile.list --vcf
    cnvCompare -i inputFile.list -c controlFile.list -w --vcf -d /path/to/genome.dict -s "counted" 

The input file format is one path to a single file on a line, the control file format is the same.
example :  
    /path/to/the/file/one.vcf
    /path/to/the/file/two.vcf
    /path/to/the/file/three.vcf



## Found a Bug ? 
Please search for any existing issues to avoid duplicate submission. 
If you can't find an issue corresponding, please feel free to open an issue in the tracker. 


## Contributing Funds 
You can donate to the [ARGAD association](https://www.helloasso.com/associations/association-pour-la-recherche-genetique-des-anomalies-du-developpement-argad), a non-profit association

## License 
Copyright 2023 : U1231 GAD - CHU Dijon 
Licensed under the AGPLv3: https://opensource.org/licenses/agpl-3.0

## Citing 
NYI