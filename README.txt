# TaIGa - TAxonomy Information GAtherer

This is a simple script that interacts with various utilities from NCBI's Entrez api in order to retrieve relevant taxonomic information for a collection of organisms. As of now, TaIGa is able to handle multiple types of Genbank format genome files, as well as a text file format list of organism names, separated by lines. TaIGa recieves a file as input, an output folder path, a valid user e-mail and one optional argument to identify the type of input file. TaIGa then uses Entrez to retrieve the TaxID, Genome ID and all taxonomic information of all taxa up to the organism name provided. Then, it builds a DataFrame and outputs it to a .csv file, so the user can visualize it as a table. TaIGa's goal is to make easier for researchers to gather mass taxonomical metadata for their projects. Therefore, TaIGa is best used when you have a big list of organisms or a big collection of genomes in a file. TaIGa is also a very cute anime character from the japanese romance animation ToraDora. You should watch it. 

# TaIGa was developed and is maintained by Maycon Douglas de Oliveira - 2019

## Dependencies:
    
    -Python 3.x, preferably 3.6
    -Biopython version 1.73 or newer (and all of its dependancies, obviously)
    -An internet connection

## How to run:

To run taiga, on your command line, from the TaIGa folder, do as following:
    $ python3 TaIGa.py [input file] [output path] [valid e-mail] --[one optional argument]

If you want further information on how to run TaIGa and what are the required and optional arguments it expects, you may run:
    $ python3 TaIGa.py -h

## Positional (required) Arguments:

[input file]: This is the full path to the file you will use as an input for TaIGa. By default, TaIGa expects it to be a Genbank format genome file with multiple records for multiple different organisms (repetitions allowed, but as of this version TaIGa will not check for them - this checking is to be implemented very soon). You can change this behaviour so TaIGa would expect: a Genbank format genome file with multiple records, all from the same organism; a Genbank format genome file with only one record; or a text file with organism names separated by line. Organism names refer to any valid taxonomic level that is available on NCBI's Taxonomy database.

[output path]: This is the full path to the output folder. This is where TaIGa will automatically create the output file, discussed below, and the missing file (also discussed below) if there is need for one. This folder must be a valid path on your system, but it doesn't need to be pre-created. TaIGa will check if the folder exists and, if it doesn't, it creates it at the provided path.

[valid e-mail]: This is just a valid e-mail of yours. Nothing will be sent to this e-mail, an neither TaIGa itself neither me will ever have use it for anything other than running TaIGa (in fact, I will never have access to this information. You may check the code yourself to confirm this). TaIGa only requires this field because it is usual to pass on this information when sending requests to Entrez. This is all TaIGa will use the e-mail for. You may pass on gibberish, if you so want, but I advise you not to. TaIGa will run fine anyways, as long as you provide something to this argument field.

## Optional Arguments (use only one per TaIGa run):

--single: This changes TaIGa's behaviour to, instead of expecting a Genbank format genome file with a collection of records of any sort, to expect only one record on the said file. This is, TaIGa, when run with this option on, will only accept your input file if it has only and no more than one record. Check your file closely before using this option. TaIGa will know.

--same: This changes TaIGa's behaviour to, instead of expecting a Genbank format genome file with multiple records from multiple, different organisms, to expect multiple records but all from the same organism. Do take care though: if your file has multiple records from different organisms and it also happens to have more than one record for the same organism, this mode won't work as expected. This will only take the first organism on the genome file and ignore all other records. If your file falls into this situation, run TaIGa in default mode (it is advisible to edit your input file before doing so, removing the duplicates). Eg.: this mode is useful if your genome file contains three 'Apis mellifera' records (each one with the DNA sequence of a cromosome). It is not useful, though, if your genome file contains three 'Apis mellifera' records and two 'Homo sapiens' records. It will only consider the first record it finds.

--name: This changes TaIGa's behaviour to, instead of expecting a Genbank format genome file of any sort, to expect a simple text file with a collection of organism names, all separated by line (linebreaks). Organism names, in this context, refer to any valid taxonomic level available on NCBI's Taxonomy database.

## TaIGa's Outputs:

If TaIGa runs successfuly, it will print a message to the screen to inform so. There is no '-v' or '--verbose' option, so there's no way (yet) to disable TaIGa's chit-chat. One way of working around this is to run TaIGa and send its output to a placeholder file, like so:
    
    $ python3 TaIGa.py [input file] [output path] [e-mail] > log.txt

This will create a 'log' file with everything TaIGa would normally print on the screen, and TaIGa will print no more.

After running successfuly, TaIGa will create at least one file at the informed output folder. If the output folder doesn't exist (or its parent folders), TaIGa will check it and create them for you. Do note that you still need to provide a valid path for TaIGa to run successfuly. Check it twice before running TaIGa. The created file will be named 'TaIGa_result.csv'. This is default and can only be changed on the source code, which you can surely do if you know what you're doing. It will be a .csv format file. To better visualize the results, import it to any spreadsheet viewer of yours.

Sometimes, TaIGa will also create a file named 'TaIGa_missing.txt'. This will be created if TaIGa wasn't able to find a valid corrected name or TaxID for any of the input organisms. On that file, you can check which organisms are those and go after that information manually.

On the subject of the 'corrected names': Entrez has a utility called 'correct spelling', which is a usefull way of checking some misspelling on some organism names before submiting them to TaIGa's search functions. The problem is, Entrez's search utility is very sensitive to anything it thinks is a misspelling. On top of that, Biopython's parsing functions can, sometimes, mess with the original organism names (way before the 'correct spelling' function is called). This usually happens with organism names containing special characters or things such as dates or researchers' names. This isn't something TaIGa handles very well as of now, so expect some perfectly normal organism names to, sometimes, appear on the 'missing' file.

## Extra Info:

As TaIGa depends on Entrez to function, it means that it also depends on NCBI's servers and, for that sake, its runtime may vary highly depending on date and time. For work hours in the US, it is common that NCBI's servers get cluttered and responses form the server get very slow. Still on that, TaIGa requires an internet connection to work and its running time will vary according to the speed of your connection. Those variables aren't something I can easily go around, so have that in mind when using TaIGa. And, last but not least, TaIGa will obviously take longer to run if you have a lot of organisms.

TaIGa is free to use, free to distribute and free to modify. TaIGa is a rather simple script and didn't take me much time to build it. It also depens heavily on pre-built modules to work properly. Being so, I won't ever complain if someone would take TaIGa and build, on top of it, a more sofisticated or efficient version of it and call it theirs. That said, I was very happy and accomplished when I finally had a stable version of TaIGa at hand, and I feel very proud of it. Because of that, I kindly ask that, if you want to use TaIGa for your work or want to build something on top of it, please, cite the original script and my name as its creator.

TaIGa's repository has an 'examples' folder on it. On that folder, you'll find an input file named 'names.txt' and both an output file and a missing information file. Those were run with the '--name' option on.

As said in the introduction, another inspiration for TaIGa's name is the cute romance anime character Taiga, from the japanese animation ToraDora. I highly recommend it.

Share love and knowledge and, on top of all, respect people.