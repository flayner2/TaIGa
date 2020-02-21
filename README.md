# TaIGa - *Ta*xonomy *I*nformation *Ga*therer

This is a simple program that interacts with various utilities from NCBI's Entrez api in order to retrieve relevant
taxonomic information for a collection of organisms. As of now, TaIGa is able to handle multiple types of Genbank format
files, as well as a text file format list of organism names or Taxon IDs, separated by lines. TaIGa recieves a file as 
input, an output folder path, a valid user e-mail and a variety of optional arguments. TaIGa uses Entrez to retrieve the 
Taxon ID, Genome ID and all taxonomic information for the provided organisms. Then, it builds a DataFrame and outputs it to 
a .csv file, so the user can visualize it as a table. TaIGa's goal is to make easier for researchers to gather mass 
taxonomic metadata for their projects. Therefore, TaIGa is best used when you have a big list of organisms or a big 
collection of genomes in a file. TaIGa is also a very cute anime character from the Japanese romance animation ToraDora. You 
should watch it. 

TaIGa is developed and maintained by Maycon Douglas de Oliveira, 2020. 

## 1 Setup and running TaIGa

### 1.1 Dependencies

- Python 3.6x
- Biopython version 1.76 or newer
- Pandas version 0.25.3 or newer
- An internet connection

### 1.2 Installing dependencies

- **Python**: 
...You might want to use your own package manager and installation procedures to run TaIGa. Most Unix distributions come with 
Python 3 installed. You can try typing `$ python` into your shell to verify it (note that the `$` character should not be typed, 
it only represents the prompt). If you don't have Python 3 installed, check the [Python official website](https://www.python.org/) 
for how to install it on your system.

- **Biopython** and **Pandas**:
...Make sure you have Python and its package manager [pip](https://pypi.org/project/pip/) installed. After that, all you need to do 
is type `pip install --user biopython pandas` and it shall be done. If you're too scared of manually installing these dependencies, 
you could also run the command `make install-deps` from inside the TaIGa root directory. Note that this script calls pip in order 
to run, so you still need to install it. Additionally, you will also need a [C compiler](https://gcc.gnu.org/) and 
[CMake](https://cmake.org/) installed.

### 1.3 How to run

**Important**: All of the following examples (as well as the "Example run" section) take in account that you have multiple Python 
versions installed on your machine and that one of them is python3.6. It might be the case that your machine only has python3.6 or 
greater (eg. python3.7) installed or that the command "python" is an alias for any python version >= 3.6. In those cases, 
running TaIGa with only "python" instead of "python3.6" is totally fine.

To run TaIGa, on your command line, from the TaIGa root directory, type the following:
    
```shell
$ python3.6 -m taiga [input file] [output path] [valid e-mail] --[optional arguments]
```

or:

```shell    
$ path/to/python3.6 -m taiga [input file] [output path] [valid e-mail] --[optional arguments]
```

**Note**: The `-m` after the Python call is very important and TaIGa will not run if you forget it.

If you want further information on how to run TaIGa and what are the required and optional arguments, you may run:
    
```shell
$ python3.6 -m taiga -h
```
    
or:
    
```shell
$ path/to/python3.6 -m taiga -h
```

For further explanation on the required and optional arguments, refer to this same documentation in the **Arguments** section below.

### 1.4 Example run

If you want to test TaIGa's functionality before using your own data (even though TaIGa should never mess with your input files),
the software is shipped with an **examples** folder, which can be found (from the root TaIGa directory) in `/taiga/examples/`. 
There, you will find a file called `inputs.tar.gz`. It contains input files for all of TaIGa's input modes. To run each of them, 
unzip it and then issue the following command (again, from the TaIGa root directory):

```shell
$ python3.6 -m taiga taiga/examples/[name_of_input_file] [output_folder] [user_email] --[mode]
```

**The --[mode] argument**: this will be better explained further down, but this argument will change the expected input file type.
If you don't use this option at all, TaIGa will expect a text format file (eg. `input.txt`) with a list of organism names separated
by lines. If you use the argument `--tid`, TaIGa will expect a text format file with a list of Taxon IDs separated by lines. If you
use the argument `--gb-mode [0, 1, 2, 3]`, TaIGa will expect a Genbank format file (`.gb` or `.gbff`) with multiple records for
different organisms (`--gb-mode 1`), a single record (`--gb-mode 2`) or multiple records for the same organism (`--gb-mode 3`).
The option `--gb-mode 0`, which is the default value, is the same as running TaIGa without `--gb-mode`, so it will expect the 
default input type (text file with names).

## 2 Arguments

### 2.1 Positional (required) Arguments:

```shell
$ python3.6 -m taiga [input file] [output path] [user e-mail] --[optional arguments]
```

**[input file]**: This is the full path to the file you will use as an input for TaIGa. By default, TaIGa expects it to be a 
list of organism names separated by line in a text-like file (`.txt`). You can change this behaviour so TaIGa would expect: a line 
separated text file with a collection of Taxon IDs; a Genbank format file with multiple records, all from the same organism; 
a Genbank format file with only one record; or a Genbank format file with multiple records from multiple organisms. Organism names 
refer to any valid taxonomic level that is available on NCBI's Taxonomy database.

**[output path]**: This is the full path to the output folder. This is where TaIGa will automatically create the output file, 
discussed below, and the missing file (also discussed below) if there is need for one. This folder must be a valid path on your 
system, but it doesn't need to be pre-created. TaIGa will check if the folder exists and, if it doesn't, it creates it at the 
provided path.

**[user e-mail]**: This is just a valid e-mail of yours. Nothing will be sent to this e-mail, and neither TaIGa itself neither me 
will ever use it for anything other than running TaIGa (in fact, I will never have access to this information. You may check the 
code yourself to confirm this). TaIGa only requires this field because it is standard procedure to pass on this information when 
sending requests to Entrez. This is all TaIGa will use the e-mail for. You may pass on gibberish, if you so want, but I advise 
you not to. TaIGa will run fine anyways, as long as you provide something to this argument field.

### 2.2 Optional Arguments:

```shell
$ python3.6 -m taiga [input file] [output path] [user e-mail] --[optional arguments]
```

**--gb-mode [0, 1, 2, 3]**: *Default: 0*. This changes TaIGa's default input type to instead expect a Genbank format file. This 
argument exepects one numeric option from the available ones. Those are:

- *1*: A Genbank format file containing multiple records from multiple, differently named organisms (eg. *Escherichia coli*, 
*Bos taurus*, *Mus musculus*, all in the same `.gb` or `.gbff` file).
- *2*: A Genbank format file containing a single record (eg. an annotation for a *COX 1* gene for *Homo sapiens*).
- *3*: A Genbank format file containing multiple records for a single organism (eg. many annotations for *Apis mellifera* genes).

**--tid**: This changes TaIGa's behaviour to, instead of expecting any sort of name-based input, to expect a text file with a list 
of valid Taxon IDs for a collection of organisms (or taxon levels). This is incompatible with the '-c' option, as TaIGa skips the 
spelling correction when run with Taxon IDs. 

**-c**: This enables TaIGa's name correcting functionality. The usefulness of this is discussed below. This is incompatible with 
'--tid'. See '--tid' above.

**-t T**: *Default: 5*. This sets the maximum number of retries TaIGa will do when fetching for taxonomic information for an 
organism. This can be very useful as Entrez will many times return broken responses.

**-v**: This disables TaIGa's standard verbose mode, so TaIGa will automatically generate a log file called **TaIGa_run.log** inside
the current working directory. This log file will contain all information about that particular TaIGa run.

## 3 Output files

### 3.1 TaIGa_result.csv

After running successfuly, TaIGa will create the output files at the provided output path. If the output folder doesn't exist 
(or its parent folders), TaIGa will check it and create them for you. Do note that you still need to provide a valid path for 
TaIGa to run successfuly. Check it twice before running TaIGa. The created file will be named **TaIGa_result.csv**. This is 
default and can only be changed on the source code, which you can surely do if you know what you're doing. It will be a .csv 
format file. To better visualize the results, import it to any spreadsheet viewer of yours.

The file will contain a number of rows equal to the number of input organisms. Each row will be named for the corresponding taxon.
Each column will be a variable for a particular taxonomic information. The first two are always the organism's Taxon ID and Genome
ID (if it has one). The rest are the valid taxon rank names available on Taxonomy and their corresponding value for each organism.
If there's any missing value (eg., lack of a Genome ID or lack of a *tribe* for *Homo sapiens*), the value will be **N/A**.

### 3.2 TaIGa_missing.txt

TaIGa will also create a file named 'TaIGa_missing.txt'. This will be created regardless if TaIGa was able to run without issues or
any missing information. If any organism happens to be missing one of the core informations TaIGa needs to be able to run (those
being a valid `Name` or `Corrected Name`, and `Taxon ID`), that organism will be outputed to this file within the correct class of
missing information.

## 4 Potential issues
 
### 4.1 Missing Taxon ID/Corrected Names when using the '-c' option

TaIGa has a certain command line optional argument (flag), the `-c` option. This will enable TaIGa's name correcting function,
so TaIGa will search for an organism's TaxID after trying to correct its name. Why isn't this enabled by default?

So, Entrez has a utility called 'correct spelling', which is a usefull way of checking some misspelling on some organism names 
before submiting them to TaIGa's search functions. The problem is, Entrez's search utility is very sensitive to anything it *thinks*
is a misspelling. On top of that, Biopython's parsing functions can, sometimes, mess with the original organism names (way before 
the 'correct spelling' function is called). This usually happens with organism names containing special characters or things such as
dates or researchers' names. This isn't something TaIGa handles very well as of now, so expect some perfectly normal organism names 
to, sometimes, appear on the 'missing' file if you choose to use the `-c` option.

### 4.2 Missing information due to server/connection issues

When you run TaIGa for a large list of organisms, there's a good chance you'll get a 'TaIGa_missing' file, probably with a decent 
number of organisms missing information. Why is it so? Well, most of it has to do with Entrez and Taxonomy themselves. Sometimes, 
requests will return broken responses, and TaIGa isn't very good at handling those (yet), so it will assume that the particular 
information it was trying to get is missing. This has a lot to do with the (lack of) stability of NCBI's servers and the Entrez api
itself, but also has to do with the sensitivity of Entrez's search functionality. One way of handling a big number of organisms with
missing data is to go to the output 'TaIGa_missing' file, get all the names there and simply run TaIGa again. This will, sometimes, 
at least reduce the number of organisms with missing data (this is a weak countermeasure to broken responses). You may have to 
repeat this several times.

Also, TaIGa doesn't handle issues with your own internet connection very well yet. So, if your connection drop during one of the
fetching steps, TaIGa will probably throw an error and stop running.

### 4.3 TaIGa is taking too long to run!

As TaIGa depends on Entrez to function, it means that it also depends on NCBI's servers and, for that sake, its runtime may vary 
highly depending on date and time. For working hours in the US, it is common that NCBI's servers get cluttered and responses 
from the server get very slow. Still on that, TaIGa requires an internet connection to work and its running time will vary 
according to the speed of your connection. Those variables aren't something I can easily go around, so have that in mind when 
using TaIGa. And, last but not least, TaIGa will obviously take longer to run if you have a lot of organisms.

### 4.4 What do I do if I have many names/Taxon IDs and I think some of them are duplicates?

Nothing. TaIGa will handle duplicate names or Taxon IDs from your input files automatically.

### 4.5 TaIGa complains about my Python version!

Don't run TaIGa with Python versions earlier than Python 3. In fact, prefer running TaIGa with Python3.6x. The reason behind this is
now legacy (TaIGa's code was heavily changed), and it lied on how differently Python3.6 and Python2 handles orders in iterable 
objects like dictionaries. TaIGa doesn't rely on that functionality anymore (at least to my capacity), but it was coded taking in
account all syntax an functionalities of Python 3.6.x, as well as its dependencies. So, for the sake of safety, TaIGa will always
check the Python version used to run it before it actually tries to run, and will throw an error and quit if your Python version is
less than 3.6.x. If you don't want this functionality, you can simply go to the file `taiga/taiga` (from the root TaIGa directory) 
and comment the line that says `helpers.sanitize_version()`. This should disable the version checking and should not break anything
else in the code.

## 5 Licensing 

TaIGa is licensed under the MIT license. You can check the information inside the LICENSE file. To make it short, I wanted it to be
free and open, so that anyone can contribute to it.

## 6 Ending regards

As said in the introduction, the major inspiration for TaIGa's name is the cute romance anime character Taiga, from the japanese 
animation ToraDora. I highly recommend it.

Share love and knowledge and, on top of all, respect people.
