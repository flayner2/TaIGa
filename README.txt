# TaIGa - TAxonomy Information GAtherer

This is a simple script that interacts with various utilities from NCBI's Entrez api in order to retrieve relevant taxonomic 
information for a collection of organisms. As of now, TaIGa is able to handle multiple types of Genbank format genome files, as 
well as a text file format list of organism names, separated by lines. TaIGa recieves a file as input, an output folder path, a 
valid user e-mail and one optional argument to identify the type of input file. TaIGa then uses Entrez to retrieve the TaxID, 
Genome ID and all taxonomic information of all taxa up to the organism name provided. Then, it builds a DataFrame and outputs it 
to a .csv file, so the user can visualize it as a table. TaIGa's goal is to make easier for researchers to gather mass 
taxonomical metadata for their projects. Therefore, TaIGa is best used when you have a big list of organisms or a big collection 
of genomes in a file. TaIGa is also a very cute anime character from the japanese romance animation ToraDora. You should watch 
it. 

# TaIGa was developed and is maintained by Maycon Douglas de Oliveira - 2019

## Dependencies:
    
    -Python 3.6x
    -Biopython version 1.73 or newer (and all of its dependencies, obviously)
    -Pandas version 0.24.2 or newer (and all of its dependencies, obviously)
    -An internet connection

## How to run:

All of the following examples (as well as the "Example run" section) take in account that you have multiple Python versions
installed on your machine and that one of them is python3.6. It might be the case that your machine only has python3.6 or 
greater (eg. python3.7) installed or that the command "python" is an alias for any python version >= 3.6. In those cases, 
running TaIGa with only "python" instead of "python3.6" is totally fine.

To run taiga, on your command line, from the TaIGa folder, do as following:
    $ python3.6 TaIGa.py [input file] [output path] [valid e-mail] --[optional arguments]
or:
    $ path/to/python3.6 TaIGa.py [input file] [output path] [valid e-mail] --[optional arguments]

If you want further information on how to run TaIGa and what are the required and optional arguments it expects, you may run:
    $ python3.6 TaIGa.py -h
or:
    $ path/to/python3.6 TaIGa.py -h

## Example run:
    $ python3.6 TaIGa.py example/example_input.txt example flayner5@gmail.com -c

Explaining: this will run TaIGa on a list of valid taxon names, stored as line separated values on a file named
'example_input.txt', defining the output folder to be 'example' and using my own e-mail. The TaIGa name correcting function is 
disasbled with the option '-c', as this usually yields better results and avoid certain bad behaviours. After this run, if it is 
successfull, there should be a 'TaIGa_result.csv' file inside the 'example' folder. The folder may also contain a 
'TaIGa_missing.txt' file if there's any missing data (in this example run, there is). If the option '-v' was used, you should 
also expect a 'TaIGa_run.log' file inside the main 'TaIGa.py' script folder. You can check the inputs and outputs of this 
example run inside the distributed 'example' folder.

## Positional (required) Arguments:

[input file]: This is the full path to the file you will use as an input for TaIGa. By default, TaIGa expects it to be a 
list of organism names separated by line in a text-like file. You can change this behaviour so TaIGa would expect: a line 
separated text file with a collection of taxon IDs; a Genbank format genome file with multiple records, all from the same 
organism; a Genbank format genome file with only one record; or a Genbank format genome file with multiple records from 
multiple organisms. Organism names refer to any valid taxonomic level that is available on NCBI's Taxonomy database.

[output path]: This is the full path to the output folder. This is where TaIGa will automatically create the output file, 
discussed below, and the missing file (also discussed below) if there is need for one. This folder must be a valid path on your 
system, but it doesn't need to be pre-created. TaIGa will check if the folder exists and, if it doesn't, it creates it at the 
provided path.

[valid e-mail]: This is just a valid e-mail of yours. Nothing will be sent to this e-mail, an neither TaIGa itself neither me 
will ever have use it for anything other than running TaIGa (in fact, I will never have access to this information. You may 
check the code yourself to confirm this). TaIGa only requires this field because it is usual to pass on this information when 
sending requests to Entrez. This is all TaIGa will use the e-mail for. You may pass on gibberish, if you so want, but I advise 
you not to. TaIGa will run fine anyways, as long as you provide something to this argument field.

## Optional Arguments:

### Input Modes (use only one per TaIGa run):

--single: This changes TaIGa's behaviour to, instead of expecting a list of names, to expect only one record on the said file. 
This is, TaIGa, when run with this option on, will only accept your input file if it has only and no more than one record. Check 
your file closely before using this option. TaIGa will know.

--same: This changes TaIGa's behaviour to, instead of expecting a list of names to expect multiple records but all from the same 
organism. Do take care though: if your file has multiple records from different organisms and it also happens to have more than 
one record for the same organism, this mode won't work as expected. This will only take the first organism on the genome file 
and ignore all other records. If your file falls into this situation, run TaIGa in default mode. Eg.: this mode is useful if 
your genome file contains three 'Apis mellifera' records (each one with the DNA sequence of a cromosome). It is not useful, 
though, if your genome file contains three 'Apis mellifera' records and two 'Homo sapiens' records. It will only consider the 
first record it finds.

--multi: This changes TaIGa's behaviour to, instead of expecting a list of names, to expect a Genbank format genome file with 
multiple records from multiple, different organisms. TaIGa does check for duplicate names and ignores them.

--tid: This changes TaIGa's behaviour to, instead of expecting any sort of name-based input, to expect a text file with a list 
of valid TaxIDs for a collection of organisms (or taxon levels). This is incompatible with the '-c' option, as TaIGa already 
skips the spelling correction when run with TaxIDs. 

### Other Options:

-c: This disables TaIGa's name correcting functionality. The usefulness of this is discussed below. This is incompatible with 
'--tid'. See '--tid' above.

-t: This sets the maximum number of retries TaIGa will do when fetching for taxonomic information for an organism. This can be 
very useful as Entrez will many times return broken responses. Default: 5.

-v: This disables TaIGa's standard verbose mode. If you use this option, TaIGa will print to the screen no more. Instead, it 
will automatically generate a log file called 'TaIGa_run.log' inside the folder TaIGa is in. This log file will contain all 
information about that particular TaIGa run.

## TaIGa's Outputs:

If TaIGa runs successfuly, it will print a message to the screen to inform so (or log it to a file if you have the '-v' option 
turned on).

After running successfuly, TaIGa will create at least one file at the informed output folder. If the output folder doesn't exist 
(or its parent folders), TaIGa will check it and create them for you. Do note that you still need to provide a valid path for 
TaIGa to run successfuly. Check it twice before running TaIGa. The created file will be named 'TaIGa_result.csv'. This is 
default and can only be changed on the source code, which you can surely do if you know what you're doing. It will be a .csv 
format file. To better visualize the results, import it to any spreadsheet viewer of yours.

Sometimes, TaIGa will also create a file named 'TaIGa_missing.txt'. This will be created if TaIGa wasn't able to find a valid 
corrected name or TaxID for any of the input organisms. On that file, you can check which organisms are those and go after that 
information manually.

On the subject of the 'corrected names': Entrez has a utility called 'correct spelling', which is a usefull way of checking some 
misspelling on some organism names before submiting them to TaIGa's search functions. The problem is, Entrez's search utility is 
very sensitive to anything it thinks is a misspelling. On top of that, Biopython's parsing functions can, sometimes, mess with 
the original organism names (way before the 'correct spelling' function is called). This usually happens with organism names 
containing special characters or things such as dates or researchers' names. This isn't something TaIGa handles very well as of 
now, so expect some perfectly normal organism names to, sometimes, appear on the 'missing' file.

## Handling Missing TaxID/Corrected Names:

TaIGa has a certain command line optional argument (flag), the '-c' option. This will disable TaIGa's name correcting function,
so TaIGa will search for an organism's TaxID without trying to correct its name. Why is this useful?

First of all, when you run TaIGa for a large list of organisms, there's a good chance you'll get a 'TaIGa_missing' file, probably
with a decent number of organisms missing information. Why is it so? Well, most of it has to do with Entrez and Taxonomy 
themselves. Sometimes, requests will return broken responses, and TaIGa isn't very good at handling those (yet), so it will 
assume that the particular information it was trying to get is missing. This has a lot to do with the (lack of) stability of 
NCBI's servers and the Entrez api itself, but also has to do with the sensitivity of Entrez's search functionality. One way of 
handling a big number of organisms with missing data is to go to the output 'TaIGa_missing' file, get all the names there and 
simply run TaIGa again (don't bother with duplicates, you may leave them there, as TaIGa is smart enough to handle them). 
This will, sometimes, at least reduce the number of organisms with missing data (because this is a countermeasure to broken 
responses). You may have to repeat this several times.

Some organisms, though, will persist with having missing data. Most of the time, those will fall into 'Missing TaxID'. The 
reason behind this has to do with Entrez's Correct Spelling functionality. Sometimes, when TaIGa runs its name correcting 
function, Entrez will change that name to a non-valid organism name, thus returning a valid corrected name, but no TaxID related 
to it. To work around this, gather all the names on your 'TaIGa_missing' file (after doing what was discussed in the previous 
paragraph) and run TaIGa again with the '-c' argument at the end. This will make TaIGa skip its spell-correcting function, thus 
solving this particular issue with Correct Spelling. In the end, your final 'TaIGa_missing' file shouldn't have many organisms 
anymore.

## Extra Info:

Don't bother with duplicate names. TaIGa will automatically evaluate if your input file contains duplicates and will remove them.

As TaIGa depends on Entrez to function, it means that it also depends on NCBI's servers and, for that sake, its runtime may vary 
highly depending on date and time. For working hours in the US, it is common that NCBI's servers get cluttered and responses 
form the server get very slow. Still on that, TaIGa requires an internet connection to work and its running time will vary 
according to the speed of your connection. Those variables aren't something I can easily go around, so have that in mind when 
using TaIGa. And, last but not least, TaIGa will obviously take longer to run if you have a lot of organisms.

Don't run TaIGa with Python2.7 or earlier versions. In fact, prefer running TaIGa with Python3.6x. The reason behind this is how 
differently Python3.6 and Python2 handle iterable objects like sets and dictionaries. In Python2.x, those objects are not 
ordered, meaning you cannot trust that the order by which you added elements to it will be maintained for further operations. In 
Python3.6, however, those iterables are ordered and indexable. TaIGa takes advantage of this particular behaviour to work  
properly, and you can actually test the difference yourself by running a list of 30~ organisms or something with Python2.7 and 
Python3.6x. Check the TaIGa_result file, particularly for rows that all of your organisms share (like 'genus' if your input is a 
list of species names). You'll see that the informations for certain rows are scrambled, because Python2.7 doesn't maintain the 
order for the items inside the objects TaIGa uses. So, to sum up, if you run TaIGa with Python2.7 it will not work as intended 
and your result file will be scrambled. Don't do it. TaIGa will know.

TaIGa is free to use, free to distribute and free to modify. TaIGa is a rather simple script and didn't take me much time to 
build it. It also depens heavily on pre-built modules to work properly. Being so, I won't ever complain if someone would take 
TaIGa and build, on top of it, a more sofisticated or efficient version of it and call it theirs. That said, I was very happy 
and accomplished when I finally had a stable version of TaIGa at hand, and I feel very proud of it. Because of that, I kindly 
ask that, if you want to use TaIGa for your work or want to build something on top of it, please, cite the original script and 
my name as its creator.

TaIGa's repository has an 'examples' folder on it. On that folder, you'll find an input file named 'names.txt' and both an 
output file and a missing information file. Those were run with the '--name' option on.

As said in the introduction, another inspiration for TaIGa's name is the cute romance anime character Taiga, from the japanese 
animation ToraDora. I highly recommend it.

Share love and knowledge and, on top of all, respect people.
