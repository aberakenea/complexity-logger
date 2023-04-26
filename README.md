## complexity-logger
'''

@author: Abera Kenea
'''
#Add a logger method to the complexity project
#Add a parameter using parseargs that sets nmax

from timeit import default_timer as timer 
from argparse import ArgumentParser 
from argparse import RawDescriptionHelpFormatter 
from datetime import datetime

import unittest 
import sys 
import os 
import hashlib
from pathlib import Path

import logging
logger = logging.getLogger(__name__) 


def initLogger(md5string):

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  
    log.setLevel(logging.DEBUG)
    for hdlr in log.handlers[:]:  
        log.removeHandler(hdlr)
    log.addHandler(fileh)      
    log.addHandler(handler)
    logging.info("+" + "*"*78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*"*78 + "+")
    logging.debug("debug mode is on")


def parseArgs():
    '''parse out Command line options.'''
    try:
        parser = ArgumentParser(description="a program to calculate Fibonacci's number", formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-m", "--max_number", dest="maxnumber", action="store", help="max value to calculate Fibonacci number [default: %(default)s]")
        parser.add_argument("-n", "--nmax", dest="nmax", action="store", help="max number of Fibonacci to be calculated [default: %(default)s]")
        
       
        args = parser.parse_args()
        global maxNumber, nmax
        maxNumber = args.maxnumber
        nmax = int(args.nmax) if args.nmax else 15
        
     
        if not maxNumber:
            maxNumber = 15
        
        if maxNumber:
            logger.info("max Fibonacci number to calculate is <" + str(maxNumber) + ">")
        else:
            logger.error("you must specify a max Fibonacci number")
            sys.exit(1)
    except KeyboardInterrupt:
        
        return 0
    except Exception as e:
        logger.error(e) 
def timeSimpleVersusDynamicFibo(nMax):
    '''
    compare times for simple and dynamic methods to calculate Fibonacci's number
    return fibonacciTimes
    '''
    n=0
    fibonacciTimes=[]
    while n<nMax:
        logger.debug("--- " + str(n))
        simpleStartTime = timer()
        fibSimple(n)
        simpleEndTime = timer()
        dynamicStartTime = timer()
        fibDynamic(n)
        dynamicEndTime = timer()
        fibonacciTimes.append({"n":n, "simple": (simpleEndTime-simpleStartTime), "dynamic": (dynamicEndTime-dynamicStartTime) })
        n+=1
    return fibonacciTimes
def generateTimingPlot(fibonacciTimes, nmax):
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    dfFibonacciTimes = pd.DataFrame(fibonacciTimes)
    dfMelt = dfFibonacciTimes.melt(id_vars=['n'], value_vars=['simple','dynamic'])
   
    plt.xlabel("Fibonacci number")
    plt.ylabel("runtime (s)",)
    sns.scatterplot(data=dfMelt, x="n", y="value", hue='variable').set(title='Title of Plot')
    timingPlotFile = os.path.join(os.getcwd(), "fibonacci_timing_0_to_" + str(nmax) + ".png")
    plt.savefig(timingPlotFile) 
def fibSimple(n):
    
    if n == 0:
        return(0)
    if n == 1:
        return(1)
   
    return(fibSimple(n-1) + fibSimple(n-2)) 
def writeTimingDataToFile(fiboTimes, nmax):  
    import csv

    timingDataFile = os.path.join(os.getcwd(), "fibonacci_timing_0_to_" + str(nmax) + ".tsv")
    print("timing data will be written to <" + timingDataFile + ">")
    with open(timingDataFile, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['n', 'recursive', 'dynamic'])
        for timing in fiboTimes:
            writer.writerow(timing.values())    
    
    print("---done")
    
def fibDynamic(n):

    
    
    if n == 0:
        return(0)
    if n == 1:
        savedFibNumbers[1] = 1
        return(1)
    
    if savedFibNumbers[n] != 0:
        return savedFibNumbers[n]
    
    savedFibNumbers[n] = fibDynamic(n -1) + fibDynamic(n - 2)

    return(savedFibNumbers[n])    
 
class CheckFibonacciDynamic(unittest.TestCase):
    
    
    def test_negative(self):
        
        nmax= 15
        global savedFibNumbers
        savedFibNumbers=[0]*nmax
        
        message = "fibonacci calculation is wrong !"
        
       

        self.assertEqual(fibDynamic(15), 55, message)
  
        

def main(argv=None): 

    if argv is None:
        argv = sys.argv
    nmax = 15
    
    parseArgs()
    
    global savedFibNumbers
    savedFibNumbers=[0]*nmax
    

    
    
    fibonacciTimes = timeSimpleVersusDynamicFibo(nmax)
    writeTimingDataToFile(fibonacciTimes, nmax)
    generateTimingPlot(fibonacciTimes, nmax)
           
    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    parseArgs()
    
    initLogger(md5String)
    
if __name__ == '__main__':

    sys.exit(main())

## the out is as following:
timing data will be written to <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\fibonacci_timing_0_to_15.tsv>
---done
+******************************************************************************+
project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\__20230425_183318__3697017c62772d71b4445895933026a4.log>
+******************************************************************************+
debug mode is on
#the code for the above out put is seen below when preveiwed:
![image](https://user-images.githubusercontent.com/130226484/234334708-3e862c2b-781e-4db0-b9d7-f7472ac0afeb.png)

#For the GC calc project: Add a logger method and calculate GC frequency for all miRNAs for by changing speciesCode for human, rat and mouse in the argument
# ''' That can be as the following in argument path for each species respectively:
#     -f "C:\Users\THIS-PC\AP\day4\data\mature.fa" -s hsa
#     -f "C:\Users\THIS-PC\AP\day4\data\mature.fa" -s rno
#     -f "C:\Users\THIS-PC\AP\day4\data\mature.fa" -s mmu
#  '''

import json
import sys
import os
from locale import atof, setlocale, LC_NUMERIC
from datetime import datetime
import hashlib
import logging


from pathlib import Path

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2023-04-11'
__updated__ = '2022-04-11'

import sequence

DEBUG = 1
TESTRUN = 0
PROFILE = 0


def initLogger(md5string):

    ''' setup log file based on project name'''
    projectBaseName = ""

    projectBaseName = Path(fastaFile).stem

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  
    log.setLevel(logging.DEBUG)
    for hdlr in log.handlers[:]: 
        log.removeHandler(hdlr)
    log.addHandler(fileh)      
    log.addHandler(handler)
    logging.info("+" + "*"*78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*"*78 + "+")
    logging.debug("debug mode is on")

def parseArgs(argv):
    '''parse out Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    
    program_license = '''%s
    i
      Created by Simon Rayner on %s.
      Copyright 2023 Oslo University Hospital. All rights reserved.
    
      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0
    
      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.
    
    USAGE
    ''' % (program_name, str(__date__))

    try:
       
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-f", "--fasta_file", dest="fastafile", action="store", help="fasta file for which you want to calc GC% [default: %(default)s]")
        parser.add_argument("-s", "--species_code", dest="speciescode", action="store", help="three character species code [default: %(default)s]")

       
        args = parser.parse_args()

        global fastaFile
        global speciesCode

        fastaFile = args.fastafile
        speciesCode = args.speciescode

        if fastaFile:
            print("fasta file is <" + fastaFile + ">")

        if speciesCode:
            print("speciesCode is <" + speciesCode + ">")

    except KeyboardInterrupt:
       
        return 0
    except Exception as e:
        print(e)
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def calcAverageGCPercent():
    '''
    calculate GC percent for each sequence and return the average value
    :return:
    '''
    totalGCPercent = 0
    sCount = 0
    for seqLine in sequenceLines:
        seq = sequence.Sequence(headerLines[sCount], seqLine)
        seq.calcGC()
        print("for sequence <" + seq.getHeaderLine() + "> GC% is <" + str(100.0*seq.getGCPercent()) + ">")
        totalGCPercent = totalGCPercent + seq.getGCPercent()
    return totalGCPercent/len(sequenceLines)

def readFastaFile(filename):
    '''
    load specified fasta file
    :param self:
    :return:
    '''
    global headerLines    
    global sequenceLines

    
    try:
        fFA = open(filename, 'r')
        fastaLines = fFA.readlines()
        fFA.close()
    except Exception as e:
        raise(e)

    headerLines = []
    headerLine = ""
    sequenceLines = []
    sequence = ""

    s = 0
    for fastaLine in fastaLines:
        if fastaLine[0] == '>':
            if s > 0 and headerLine.startswith(speciesCode):
                headerLines.append(headerLine)
                sequenceLines.append(sequence)
                sequence = ""
            headerLine = fastaLine[1:].strip()
        else:
            sequence = sequence + fastaLine.strip()
        s += 1
    if headerLine.startswith(speciesCode):
        headerLines.append(headerLine)
        sequenceLines.append(sequence)
    return len(headerLines)


def main(argv=None): 

    
    if argv is None:
        argv = sys.argv

    
    parseArgs(argv)
   
    n = readFastaFile(fastaFile)
    print("found <" + str(n) + "> sequences")

    avGCPercent = calcAverageGCPercent()
    print("average GC % = <" + str(100.0*avGCPercent) + ">")

    global filename

    filename = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"
    
    

    print ("The number of sequences are: ", str(n))
    
    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    parseArgs(argv)
    
    initLogger(md5String)    
    print(md5String) 

    
if __name__ == '__main__':

    sys.exit(main())
## for hsa the out put are as the following:
## the out put of logger method for GC calc of "hsa":
2023-04-25 20:17:42,519 - root - INFO - +******************************************************************************+
2023-04-25 20:17:42,520 - root - INFO - project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\mature__20230425_201742__3697017c62772d71b4445895933026a4.log>
2023-04-25 20:17:42,520 - root - INFO - +******************************************************************************+
2023-04-25 20:17:42,521 - root - DEBUG - debug mode is on

average GC % = <70.05512368618675>
The number of sequences are:  2656
fasta file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa>
speciesCode is <hsa>
+******************************************************************************+
project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\mature__20230425_192641__3697017c62772d71b4445895933026a4.log>
+******************************************************************************+
debug mode is on
## the number of mature_uniqueseeds for "hsa" are:2094 
write unique seed sequences to fasta file
output fasta file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature__uniqseeds.fa>
['GAGGUAG', 'UAUACAA', 'UGUACAG', 'UGUACAA', 'UAUACGA', 'UAUACGG', 'UAUACAG', 'AGCAGCA', 'AGGCCAU', 'CAGUAUU', 'AAAGUGC', 'CUGCAGU', 'AAGGUGC', 'CUGCCCU', 'GUUUUGC', 'GUGCAAA', 'CUGCAUU', 'AGCUUAU', 'AACACCA', 'GUUCUUC', 'AGCUGCC', 'GGGUUCC', 'UCACAUU', 'GCCUACU', 'GGCUCAG', 'GGCGGAG', 'AUUGCAC', 'UCAAGUA', 'CUAUUCU', 'CUGUUCU', 'GGGCUUA', 'UCACAGU', 'AGGAGCU', 'ACUAGAU', 'CUGAUUU', 'AGCACCA', 'GUAAACA', 'UUUCAGU', 'GGCAAGA', 'GCUAUGC', 'AAUUUAG', 'UGCAUUG', 'AAUGUUU', 'GGUUGGG', 'GGUGGGG', 'CUGCUGA', 'CAAUAAA', 'UCAACGG', 'UUGGCAC', 'AUCAUGU', 'ACCCGUA', 'AAGCUCG', 'AAGCUUG', 'AGUUAUC', 'ACAGUAC', 'CUGGUUU', 'UGGUUUC', 'GCUUCUU', 'GCAGCAU', 'CAAAUGC', 'CGGAUGU', 'UGCAAUG', 'CAAUAUU', 'UGACCUA', 'UGCCAAU', 'AGGUAGU', 'AACAACA', 'GGGUAGA', 'UCACCAC', 'GUCCAGA', 'CCAGUGU', 'CAGUAGU', 'AGCUUUU', 'UAAGACG', 'UUUUUGC', 'AGCCCUU', 'AAGUUCU', 'CAGUGCA', 'UGGGAGA', 'CUACAGU', 'GGAGACG', 'UGUGUGG', 'GGAAGAC', 'AACAAAU', 'ACCCUGU', 'AAAUUCG', 'CAGAUUC', 'GGCAGUG', 'AAUCAGC', 'ACAUUCA', 'CCACUGA', 'UCACUGA', 'ACCAUCG', 'UUGGCAA', 'GGUUCUA', 'AUGGCAC', 'UGAAUUA', 'GCUACAA', 'CGUGUCU', 'GGCAACA', 'GUGGUUC', 'UGAAAUG', 'UCCCUUU', 'CUGGGAA', 'CCUUCAU', 'AUUUCAG', 'GCCCCUG', 'UGUGCGU', 'CAGGGAC', 'CCUUGGC', 'AACAGUC', 'CCAUCGA', 'GCCUGUC', 'CAGCAGG', 'CUGUCAU', 'AAUCUCA', 'CACAGUG', 'ACUGCAU', 'AUCAGUU', 'UGUGCUU', 'UGGUUCC', 'AUGGUUC', 'GAUUGUC', 'GAGUUGA', 'CCUGGCA', 'GCUACAU', 'UCAGUAG', 'GUGUAUU', 'GUCAGUU', 'CAAGUCA', 'AAAUGGU', 'AUCUUAC', 'AAUACUG', 'UGCGCAA', 'GGAAUGU', 'GAAUCAU', 'GAGCUUA', 'UGGGAGG', 'GGAGUGU', 'ACGCCAU', 'GUGUUCA', 'AAGGCAC', 'CCCUGAG', 'CGGGUUA', 'GGGGCCG', 'CUCUUUU', 'AGUGCAA', 'CCGUGGC', 'GCUGGUA', 'UUGGUCC', 'AUGGCUU', 'AUAGGGA', 'UGUAGGG', 'CGGGUAU', 'UAUUGCU', 'GCUGGUG', 'CUAUUUC', 'AGUGGUU', 'ACCACAG', 'AUCUUCC', 'AACACUG', 'AUAAAGU', 'GUAGUGU', 'GUGCAGU', 'GAGAUGA', 'GAUAUCA', 'ACAGUAU', 'UCCAGUU', 'GAUUCCU', 'GGUUCUG', 'UGCAUAG', 'CAUUUUU', 'AACGGAA', 'CUGCGCU', 'CUUUGGU', 'UAAAGCU', 'CAGGUGA', 'CACAAGU', 'AUUAUUA', 'CGUACCG', 'UGAAGCU', 'CGGAUCC', 'GUGACUG', 'CUGUGGG', 'CUCCAUU', 'AUCAUCG', 'CUACUUC', 'GAGAACU', 'CUCUGAA', 'CUGGCUC', 'GGGAGGG', 'CUCCCAA', 'UGGUACA', 'AGGUUAU', 'AUCAUAC', 'GGACGGA', 'GGAGAGA', 'GGGGCUG', 'AAAGAAU', 'CCCAAAG', 'AUCCCUU', 'UCCCACA', 'GAUAUGU', 'UAUAUAU', 'GGGUCUU', 'ACUGGCC', 'GUAACAG', 'CCUUCUC', 'AAAGCUG', 'GUCUUAC', 'CAUACUU', 'UAAUGCU', 'UCCUACA', 'CAGUGGG', 'CGCACUG', 'GACCGAU', 'CUUAAAC', 'AAGUGCU', 'CGGUUAU', 'GAAUUGU', 'AGGCAGU', 'AAUCACU', 'AUCACUA', 'GGUUUAC', 'AUGUGGG', 'CUCUGAC', 'GGGCCCC', 'AGGGUUG', 'CUCUUUC', 'UAUCAGA', 'CCCCCAG', 'AUCCUUG', 'ACACACC', 'GGGUGGA', 'GGGACUU', 'AAUGCCC', 'CUUUAAC', 'UUAACAU', 'CUGUUGC', 'GUGGAUA', 'ACAUAGA', 'GAUCGAC', 'AUAAUAC', 'AGGUCAC', 'CCUGCUG', 'CUCAAAC', 'AGUGCCG', 'CUCAAAU', 'CUCAAAA', 'UAUAAUA', 'UUAUCAG', 'CGACGAG', 'UUGUUCG', 'UAGAUUC', 'UCAUAGA', 'GAGGUUG', 'UCACACA', 'UCCUGAC', 'CUGGACU', 'GGUAGAC', 'AUGUAAC', 'GGUUGAC', 'AUGUAAU', 'GCGAGGU', 'AUACAAG', 'AAGUUGU', 'AUCAUUC', 'GAUCAGA', 'CAGCACU', 'UAUAAAG', 'CCGUCUC', 'CUCUGGG', 'CAAAGCA', 'GGGGGGC', 'UGGCCCU', 'GGGGUGC', 'CUCACAC', 'AACGGCU', 'UCCUAUA', 'GGUGGUC', 'ACAUUAC', 'CGAGGAG', 'UAGACUG', 'AGUUCUG', 'UAGGUAU', 'CCCCUGG', 'GCAUCCC', 'CCACUGC', 'ACAAUAU', 'CCAGCAU', 'CCCUGUC', 'GAGCGCC', 'CAAGAGC', 'UUUUCAU', 'CUAGUAG', 'CUGACUC', 'CCCUGAA', 'GUCUGCC', 'UUCCUAG', 'CGACAGC', 'GAGGGGC', 'GCUCGGU', 'AAAACGU', 'AUGACAC', 'UCGGGAA', 'GCCCUAA', 'CUGUAGU', 'UGCAUAU', 'UUUGCGA', 'UUGGGAA', 'GUCUUGC', 'AGGUCGU', 'ACGGUGA', 'UCAUGAU', 'AGGUUUU', 'AACCGUU', 'ACUGUUU', 'UCAUCUG', 'GGUUACC', 'AAUGUUG', 'GGUCGAC', 'CUUCACC', 'GGUUGUC', 'AUAUAAC', 'AGACGGG', 'CACUCCU', 'CAGGCUC', 'GAGGCUG', 'UCAUACA', 'CCUGUAC', 'GGGGCAG', 'UGGUUAU', 'CCAGAUA', 'UGAAAGG', 'GUCGUAU', 'UGACAUC', 'CAUGGAU', 'AACCUGG', 'GUGGGGA', 'UUAUGCA', 'UGUCUUU', 'AUGUGUA', 'CCCUGUG', 'UCCUAUG', 'GAGGUAU', 'GGACCUG', 'UGUACAU', 'GAAGGUC', 'CUUGGAG', 'UGGAUGG', 'GAAACAU', 'AAGUUGC', 'AACAAAC', 'GAGUAUU', 'GGGGUUU', 'AAACCAC', 'CACCGGG', 'ACUCAGC', 'AGUGCUG', 'UUCAAGC', 'AAGCACC', 'UCAAGAU', 'UCUCCAA', 'AGUGCCU', 'CUCUAAA', 'AGUGCUU', 'UCUAGAG', 'AAGUGCA', 'UCCAGAG', 'UCUUGAG', 'AAGGCGC', 'AACGCGC', 'AAAGCGC', 'CUCUACA', 'CUCUGGA', 'UACAAAG', 'CUCUAGA', 'UCGUGCA', 'CUCCAAA', 'ACGCACU', 'CUAGAGG', 'CAAAGUG', 'UCUGGAG', 'GCUUCCU', 'AAGCGCU', 'UGCAAAG', 'UCUCGAG', 'UAAGACU', 'ACAUCAC', 'AAUCCUU', 'UGCACCU', 'AUCCUUU', 'AUGCACC', 'UCCUUGC', 'UUGGGGA', 'AGCAGCG', 'GGGUAUU', 'GACCCUG', 'GGAGUGC', 'GGAGCCA', 'GUCAACA', 'UCACAGG', 'AAAUUUC', 'AUUCAGG', 'UUUGCAC', 'ACUCCAG', 'GAUUGUA', 'ACUGCAG', 'GAUUGGU', 'ACUCAGG', 'UUGAAAC', 'ACUCUGG', 'UUGACAC', 'AUGCCUU', 'CUCCCAC', 'AUGUGCC', 'CAGUCCA', 'GAGAAAU', 'UUCUGCA', 'CAGUAAA', 'CAGCAAA', 'GUAGAUU', 'AUCGUAC', 'CGACCCA', 'UUUAACC', 'ACAGGUG', 'AAACGGU', 'CUAGUCC', 'GGGACGG', 'GGGUAAG', 'AUGAGCU', 'UAUUACC', 'GAGCUGC', 'AAAGUAA', 'UCAAGGA', 'AAAGUUU', 'AAGUAGC', 'GGCACGG', 'GUAUGUU', 'UGUAUAA', 'AAAUCAA', 'GUUAAUG', 'AAGGUAA', 'GAAAACA', 'GAGUUGG', 'UCCGCUC', 'UGAAGUG', 'GAGUGUG', 'ACGCUCA', 'AGCCAGU', 'UUCUAAU', 'AGAUGUG', 'AGAUAAA', 'UUCUUGU', 'CGCGGUU', 'UCAUUUG', 'AAUGAUU', 'UGAGAAU', 'CUUGUGU', 'UACAGUU', 'AACUGGU', 'AAAGAGG', 'UAUGGUU', 'CAGUUCC', 'UAGCACA', 'GGGCGUA', 'AAAACUG', 'AUGCAUU', 'UUCCAUA', 'AAGAACC', 'UGGCCAC', 'GAGAACC', 'CAGAACA', 'GUGCCUG', 'AAUUUUA', 'GACCAUG', 'UGUGUCA', 'GGCACCA', 'GUCUCUG', 'AAGUGUG', 'AGCCUGC', 'GUGUCAC', 'GGUUCUC', 'CGGUGAU', 'ACGUCAU', 'UUGUGUC', 'CUUACAG', 'GGUCUAG', 'ACACGGG', 'ACACACU', 'GGCUGCG', 'AAAUCCC', 'GAAGGCA', 'AACUACU', 'UUCAAAU', 'GGGGUGG', 'GGGUGUU', 'GAGCUAA', 'CGAGGAC', 'CUGGGCA', 'AACGCCU', 'GGGGUCC', 'CCGAGCC', 'GUCAUUG', 'AAAAAUC', 'GACUUCC', 'AACUCUA', 'CUGGGAU', 'ACCUGGA', 'UGGAGAU', 'GCUAGCA', 'CAGUCUG', 'UCCCUUG', 'AGUACCA', 'ACAAGGU', 'GGGGGAA', 'ACUAUAG', 'GCUGUCU', 'UGAGUCU', 'CUUUUCU', 'UGCUGAC', 'CUAGUAA', 'GGGUUUA', 'UUCUCCC', 'GUAUUCU', 'GACCUGG', 'UGUCUGC', 'UAAUAGU', 'ACCAGCA', 'CUUGGGC', 'GUGCUUG', 'CUGGGGG', 'GGGAUCG', 'UCGCUGC', 'UGAUCCA', 'AAGACAU', 'UCCCUCU', 'GACACAU', 'CUUGUAU', 'GUGUGGC', 'CUAGGCU', 'AGCAGCU', 'UGGCUGC', 'AGUGUGC', 'AACCUGU', 'GGAGGCA', 'UUAGGAU', 'AAGGAAA', 'AACCCUA', 'AUGGCGC', 'AAAAACC', 'GCCUGGG', 'CCCACGU', 'GGCGGGG', 'AGCCACA', 'UGUUGAA', 'UCACUGG', 'AGUAGAC', 'GGUGGGC', 'AUGUCUG', 'GAGGUUA', 'UAAUACA', 'GGUUGCC', 'AUAUUAU', 'GCUCAUC', 'GACAACU', 'GCAGGUU', 'GCGGAGG', 'GGACCUU', 'UUGGUUC', 'ACCCAUU', 'CCUCCUG', 'UCAACAG', 'CGGGGAU', 'GUGACAG', 'AUGGUUG', 'UUGUGAC', 'AAGUCUU', 'GGAAGCC', 'CCGGUUC', 'GCGCCUC', 'GUCACUC', 'GCACCAU', 'CUGCUCA', 'UGAGGAC', 'CCCACCU', 'UAGGGCC', 'AGUGGGG', 'UCCGUUU', 'GCAAAAU', 'CAAAACU', 'GCUCUAG', 'UGCAGCU', 'CCCUAUC', 'GAGGAUA', 'UGCUAGU', 'CUACAAA', 'GAGACCU', 'UGGGAUC', 'GGAGGAA', 'CUCCAGC', 'GGGCUGG', 'AGUAACA', 'UCCCUGA', 'UUCCUCA', 'UCAUUCG', 'AUCUGGG', 'UUUGUGC', 'CAGGUGC', 'CAGAGUG', 'GGAGGAG', 'CCAGUAC', 'GGUGCGG', 'UGUAUGC', 'GCAGAAG', 'GCAACGA', 'GUGGCAC', 'ACUGUGU', 'UUUGCAA', 'GGCCCCA', 'UGCCCUG', 'ACUUGGA', 'GCAACUU', 'ACUCAAA', 'ACUGACA', 'ACUGGCU', 'AAGGAUU', 'AUGGCUG', 'UAAUAUC', 'AUACCUC', 'CUGGAAA', 'GGAUUUC', 'GGUGGUU', 'AACUAGA', 'GGAAACA', 'UGUGCGG', 'CUAAAUG', 'GCGGGGC', 'UGUUGCC', 'CCAUUAC', 'GGCAGCG', 'UAGAGGA', 'CCUCUUC', 'UUGGGAG', 'UGAACGG', 'CCAGGAG', 'CAGGAAC', 'GAGACUG', 'AACAUUC', 'UUAGCAG', 'GGCUCUG', 'AAUCUCU', 'CACACUU', 'GGGAGCU', 'UAGUGAG', 'CAGCAGA', 'GAGUCUU', 'GUGCGCA', 'GUCUACU', 'CAGUUAC', 'CAGUAGA', 'UGAGUCA', 'UCCGCGC', 'GCCCUUA', 'GGGGAGC', 'CCUGGGC', 'AGGCAGG', 'ACCCGGC', 'CUUCUCU', 'ACAUGGC', 'UGACUGU', 'AAUUAUU', 'UGUAUGU', 'AGGGUCA', 'UGCUCAC', 'AGCAUUC', 'GACCCAC', 'UUCCGGC', 'CGUCGCC', 'AGGGUCU', 'ACUGUAG', 'CUGCAGC', 'UGGGUAC', 'GAGCCCC', 'UGAGGGC', 'CACCAGC', 'UGGGGCC', 'GUGCCAC', 'UGGGCGG', 'CACACCU', 'UGGGUAG', 'UCUCACC', 'UGUCUGG', 'GUGGGAG', 'GAGCCCU', 'CGGCCUG', 'GAGUGAC', 'CUCUUCC', 'GGGGGCG', 'CCUUCUG', 'UGAGUGG', 'UUCCUCG', 'UCCUGAG', 'UGCCAGC', 'CCGGAGC', 'GUGGCCC', 'CGUGGCC', 'CUGCAGG', 'GUUCAUG', 'GGCAGGG', 'CAGCUGG', 'CACUGUU', 'AAAAGCA', 'AAAAACU', 'AUCUCAC', 'GCAGGAC', 'GCUGGAU', 'UCUAGCC', 'GGAGUCC', 'GGAUUUU', 'GGCCCUG', 'AAAGUAC', 'GGGUGGU', 'GUGAGGU', 'UAGGCCG', 'UCUGGAA', 'AAAGUAU', 'UGGGACA', 'UUAGAGA', 'UUGAGGC', 'CUCACUG', 'UUUCAAC', 'ACUGGAU', 'GCAAAAG', 'AGUAGUU', 'AGUGAUC', 'AUGGAUU', 'CCCGUCC', 'CCCGGGA', 'CCUUCUU', 'GGAGGGA', 'CGCCCUU', 'CGGUGCU', 'CAUUUUC', 'CUCUAGC', 'GCUUUGC', 'GAGAAGA', 'GGAUGAG', 'GGCAUUG', 'GUGAAUG', 'GUUAGGA', 'UCCCACC', 'AAACUGU', 'UGGAUAA', 'UGGGUGA', 'UGGUACC', 'AAAAGUA', 'AAAGGUA', 'AGGAUGU', 'CUCAGGG', 'CCUGUUC', 'CUGUUGA', 'GGGCGUG', 'UGGACUG', 'AUGAUGA', 'UGGGGGA', 'AAAGAGC', 'AAUUGCU', 'AAUAUAU', 'ACGUAGA', 'AGCAAAA', 'AGUACUG', 'CAUAUUG', 'CGCCUCC', 'CGUUUGC', 'CUAUACA', 'CAGAUCA', 'GGACUGC', 'GGGAACG', 'CGCGCCC', 'GAAGGAA', 'AAAUGAG', 'ACCACUU', 'CUGGCUA', 'AUUCAUU', 'CACCUCC', 'CGUUGGC', 'CGACCGG', 'CUCGGCG', 'UCACAAG', 'AAUGUCA', 'UCUCAAG', 'AGGGAGG', 'AUGAUGC', 'AGGACAC', 'CAGACAG', 'UCGGCGC', 'CCCUCCG', 'CCCGCGU', 'GCUGUAA', 'AAACCGU', 'GGCCCGG', 'CCUGCGC', 'CAUAGCC', 'CCAGUGC', 'GAGGCAG', 'CGGCCGC', 'GAGUGCC', 'GCAGGGG', 'CAGUCCU', 'AGGCAGA', 'GAGUACC', 'ACCAGGC', 'UCAUUGC', 'ACCCAGA', 'CUGCCCC', 'CCUGUGC', 'GAGGGGU', 'CCUUGCC', 'CCCAGGG', 'CAGGCCA', 'CCGUGCA', 'CUCCUGC', 'GUUUUGA', 'UGUUAAU', 'UGUAAUA', 'UGGGGAA', 'AGUCCCU', 'GAGCCUC', 'GCUUCCA', 'AUCAGAA', 'GUUCUUA', 'CUCCCAU', 'GUUCUCU', 'CUGGUGC', 'CCCUCUG', 'CUGCAAG', 'GCGCGGG', 'GACAGCG', 'AGAGCAG', 'UUUUACC', 'AUCAUGG', 'GCCUCUU', 'GGACCCA', 'UUCCGCC', 'GGGCCUG', 'UUAGGGC', 'UAUGGGU', 'GCCUGGA', 'GACACUA', 'UAGGACU', 'GGCUUUU', 'CUGUCUG', 'ACAGCAA', 'CCUUUGC', 'AAAUAGA', 'UUGGGAC', 'AGAGAAU', 'UCGCGGG', 'CUUUCCU', 'UGGCCAA', 'AGAGGAA', 'GAGGGAC', 'AUCUGGC', 'UCAGGGC', 'CCCCUUC', 'CUGGCAA', 'AACUAAU', 'ACCCAGU', 'CUGCACC', 'AAAGAAC', 'GAUGGAU', 'GCCUAGG', 'UACACAU', 'UGACUGA', 'GGCCCAA', 'CCUGAGG', 'CUGUAGC', 'GUGGACA', 'CCUGAAU', 'GCUUUUG', 'AAAAGUG', 'AAAACCA', 'AGGGCGG', 'AGGCCUU', 'UAACAUU', 'AAAGACU', 'GGGGACC', 'UAUACCU', 'ACUCCAA', 'GAUAUUU', 'GCGACAA', 'AUGCUAG', 'GUUGGGC', 'GCUACAG', 'GGAAAAA', 'UUGUAUG', 'AACCUCG', 'UGGGGAG', 'GUGGGGC', 'CUGAUCC', 'UUGCCUC', 'GUGUUAG', 'GGGAAAG', 'UUCCUGC', 'AUAUCAG', 'AGAAGGG', 'AAGAUCU', 'UCCCACU', 'UCAGCCA', 'UGCCCUA', 'CUGCAGA', 'AGGGCUU', 'AGGAUUA', 'GCUUUCU', 'GAGCUGA', 'UGAUAAG', 'UAGGGAG', 'CCCUACC', 'AUAAAAU', 'GUGACUU', 'GGUGGAU', 'GCAGACA', 'AGUUCUA', 'AGGACUG', 'UGGGGUU', 'GAUGUAU', 'GCCCUGC', 'AAGGAGG', 'GGAUGGU', 'CCAAUAC', 'AGUGAGU', 'GGGGAGA', 'CUGGCCU', 'GUGUACA', 'GCACGGC', 'GGGCGCG', 'GAAGGGG', 'UUCCAGA', 'GGGGCGG', 'UCGGGCC', 'CUUCUGU', 'CCUCUCU', 'AAGUCUC', 'GAAGAAG', 'CAACAAA', 'CAGCACC', 'AGGCGUC', 'CACGCGG', 'UGGCCAU', 'GAGGCUU', 'GCCCCAU', 'CCUUGGG', 'AAGCUGG', 'CUGGCCA', 'GUGGAAG', 'UCUCUGG', 'GGGGACG', 'CUGGGAG', 'UCUGAUC', 'CCUGCGU', 'GCCAGCC', 'GCUCUGC', 'GCGCCGG', 'GAGGCGC', 'UGGAGUC', 'GGGACUG', 'AUCUGAG', 'ACCUUGC', 'GGAUAUG', 'GGAAGGG', 'UGUGGGC', 'GCCUUCC', 'AGCCUGA', 'GGAGUCU', 'CCCACUA', 'CUGGUGA', 'CGGCAUG', 'CAGUGUG', 'UCUGAGC', 'CUAGACA', 'AUGUUUU', 'CCCUGGA', 'CAGCAUU', 'AAAGAGA', 'GCCUUGU', 'GCCCCCU', 'CGCUUUC', 'GUGAGGC', 'UCUGGGA', 'ACUGUGG', 'GGAUUCU', 'CAUUGCC', 'UAGCGGU', 'AGCCCCA', 'CCUGAGA', 'UCUGACC', 'CAGAGGU', 'CCCGCCA', 'AGUUGGG', 'UUGGGGC', 'GGGCAUG', 'CUGAGAA', 'CCUGGAG', 'AGUGUUC', 'GCCACUG', 'UGCACUU', 'GUUCCUC', 'GCUUGCA', 'UAGGAGG', 'CCAGCUC', 'UCCCCAG', 'UUGUCCU', 'GCUCCUC', 'CAGGCAC', 'UUCUAAG', 'CUCAGUC', 'CAGGGAG', 'GGGGAAG', 'AUUCAAC', 'UGUUCUC', 'UCAGUGA', 'CAAUUAC', 'AGCAGUC', 'GGUCCCG', 'CAGUUCU', 'UCUCCUC', 'UAGGGGG', 'AGUGUAG', 'AAAAUUU', 'CGGCGAG', 'GGGGCUC', 'GGCUCAC', 'CCCCACU', 'CUCCCUU', 'CCCUGGG', 'CAUUGUG', 'GCCCUCC', 'UCAGCAG', 'CUGAGAC', 'CUCAGAU', 'CACCCAG', 'CAGUUUU', 'GAGGAUG', 'CUCCGUG', 'UAGUGAA', 'AAAGUGA', 'AAUCGGA', 'UGUGAAG', 'GUUGUAC', 'CAAAAAA', 'CACUUGG', 'AGCCUUC', 'CUCUCGG', 'GAGGGCA', 'AUCAGCA', 'GUCUACA', 'GGACCAU', 'UGGGCUG', 'CACCCUG', 'GCGGGUC', 'AGGCACG', 'CACCUGA', 'GGCAUGG', 'AAAUGAA', 'GCCGCGG', 'GGGACCU', 'GGUGUGU', 'AUAGCCC', 'GGCUGGA', 'ACUGGAC', 'CUUGUCG', 'GUGUCCC', 'UUAAGAA', 'GAGUGUU', 'CUGACAG', 'AAAAUGA', 'CUGGUCU', 'GAGCACC', 'ACUCUGU', 'CUCAGGA', 'GCAGGUG', 'AAGACCC', 'CCUUCCU', 'AUGUAGA', 'GAGCUCA', 'UCAAAUA', 'UGAGACU', 'UUGUAGA', 'AUGGGGC', 'AUCUCUA', 'AGUGGCC', 'UCGUGGG', 'CCGUACA', 'UGCAGAG', 'UUCCCCC', 'ACUCACU', 'UUUGCAU', 'AGUGGAU', 'UACUUCU', 'GAUGAUA', 'GCGACAU', 'UAGACCU', 'UUCCUAC', 'UCUGUAA', 'GUGGCAA', 'AUGGAAA', 'GUGAUAU', 'CCUGGAC', 'GUGGAUG', 'CCAAGUC', 'CUGCUGG', 'UUCCACA', 'GUAUCCG', 'AAGGCAG', 'GGUGCUC', 'AGCAAUG', 'GUCCUCU', 'AAGGCAU', 'GUGUGGA', 'UGUCCAU', 'AACGCAU', 'GACAUCA', 'AGGAACC', 'UGAGGAA', 'CUCGGAC', 'CAGGGCC', 'CAGAGAA', 'GAGGAGA', 'CUGAUUA', 'CUCUGAG', 'CAAGGCC', 'ACUAGUA', 'UAUGUAU', 'AGAGAAC', 'CUCCAGU', 'GGCCAAA', 'CCUAUCA', 'AGGUAGA', 'CUUCAAC', 'UGUCCUA', 'GAAGCUC', 'GAGGAAC', 'AGGCUGA', 'CAGGUGU', 'GCUCAGG', 'GUAGAUA', 'AAGGGGU', 'CAGGCGG', 'AUUCCCU', 'ACGCGCA', 'UGGGUUG', 'AGCCCGG', 'UACACAC', 'AGCAAUA', 'UUCAGAU', 'AGCCCCC', 'GUGCAGC', 'UCGGGCU', 'GGGCAUA', 'ACUUAGC', 'GUUCCCU', 'UGUGCCU', 'CUUACUC', 'CCUGUCU', 'AAAGCAU', 'GUUGCCU', 'UAGGCAC', 'CUGGGCU', 'GAGUUAA', 'GUUGGGA', 'AAGAUGG', 'CUGAAUA', 'AAGGAGA', 'AAAACGA', 'CGACUCU', 'AAGACUC', 'GUCCCAC', 'CAGGAGU', 'GGAGAAG', 'UGGCCAG', 'CAGGACA', 'GGGCUCA', 'ACAGGCU', 'GUCGUGG', 'CCGGACA', 'UGGAGGC', 'UCGAGUU', 'GAUUGUU', 'ACGGCAA', 'UUUCCCU', 'AGGGCUG', 'GUGGGGG', 'GCUCCUU', 'GUCCCGG', 'GGGGAUU', 'GGUAGAG', 'UGAAUUC', 'AGCUUGG', 'GAUCCGA', 'GGGUGUG', 'CUGGUGG', 'CACAAGG', 'UAGUGGU', 'GCUGGAG', 'GACACGG', 'AGACUGG', 'AGGUUUG', 'GGCGGCG', 'GAGCAGA', 'CUCCCUC', 'GGCAAAC', 'GGGAACU', 'UAGUGCU', 'UAGUCUC', 'UGUGGCU', 'AAGGGAC', 'AGGAAGG', 'UAUUAAG', 'UUAAGGA', 'GGAGGUG', 'GCGCGGC', 'GCCAAGU', 'GAGUGGG', 'UUCUAUU', 'AAAGGCG', 'CCGCCUG', 'AACGGCC', 'CUGGGCG', 'GAGCUGG', 'GGGGCUA', 'CUGGUAA', 'AUGUGGA', 'GAAGGCC', 'CAGACUG', 'AUGUAAA', 'AGGAAAC', 'UCCGGGA', 'AGACUGA', 'AUGUGAC', 'CUGAUGA', 'UUAAGCA', 'GUGACAA', 'GGCUGGG', 'AAUGGGU', 'AGGAUGG', 'CGGGGCU', 'CUAAAGG', 'GAGGGAG', 'AAGAACU', 'AGGGCCU', 'GACUGAC', 'CAGGCAG', 'GGACUGG', 'AAUAUGA', 'AGCAGUG', 'CUGCGUG', 'UGGACAG', 'CUAAGGA', 'GACUCUG', 'ACCGAGA', 'UAGCAGC', 'GAGACAG', 'GGGGGAU', 'GGUCUGC', 'CAUUAUA', 'GGCCAUC', 'UUGGACU', 'CCAGCAG', 'UGGAGAA', 'GGAAGGA', 'GAUGGAG', 'UGGACCU', 'GUGGUAG', 'CGUGCAU', 'AAAGGCA', 'GAGCCGA', 'CUGAACU', 'UAGUCCU', 'GCGGCGG', 'CAAAGUA', 'AAGGUCA', 'GAGGCUA', 'AUAGAGA', 'UGCUUCA', 'UGGAAAG', 'AGGGCAG', 'GCCGCCC', 'UAUGCCU', 'GGAGCUA', 'GGCGCGA', 'CUUGAAG', 'ACUCGUG', 'ACUAACU', 'CUCGGCU', 'CUGGACA', 'UGCUAAG', 'CACUCUC', 'GGGCCAG', 'ACCCCCU', 'GCCCAUG', 'UGGCAUC', 'ACACAUG', 'CCAGGCA', 'GACAGUA', 'UUGUCCC', 'GUGGGAC', 'GGGCGAG', 'CUGAGGC', 'CAGGCCU', 'GGUAGAA', 'GGGGACU', 'UUCUGUU', 'GGAGUUA', 'GUGGGAU', 'ACCGGGG', 'CCCUCGU', 'GGGCUGA', 'UGAGUGU', 'UGCCAUG', 'UUCUUCU', 'GCAGCUC', 'ACUAGCU', 'AGGAUCC', 'UAGCCAA', 'AAGAUAG', 'GCUGAGC', 'UUCCGGU', 'UCGGCCG', 'UACAUGU', 'AUACAAU', 'CUGGGGA', 'CCCUCCU', 'GGGAAAA', 'AAAAUCC', 'GAUGUCC', 'GAAUUGC', 'GUGUCCG', 'AGCGACC', 'GAAGUUA', 'CCGAAGA', 'UAGUGCA', 'UACACAG', 'UGGGCUC', 'GGGCUGU', 'UGUUCUU', 'CUGUGAG', 'AGGUAUU', 'CUGUGAU', 'GAACUCU', 'CUGAAUU', 'ACGGGAA', 'CUGAGUU', 'GGAGAUC', 'UCUCUAC', 'GUUGCAA', 'CCAGGGC', 'AUCUGCU', 'AGCCCUC', 'GGCUGUU', 'UCCUGGG', 'AGGGGCA', 'UGAGGAG', 'AGCAGGC', 'CAGCCCA', 'UCCUCCA', 'CAGCCAC', 'UACUGUG', 'GAGAGUG', 'GGUGUUA', 'AAAUGGA', 'AGGAGGC', 'GAUCUCA', 'GCAAGAC', 'GUCAGUG', 'CAAAAUG', 'GAAGAUU', 'AUUUACU', 'ACAGGAC', 'GUAGUUG', 'ACACUAG', 'CAGUCAC', 'CAAUCAC', 'GCGGGGA', 'CCCCGGC', 'GCCCGCC', 'GAGAUGC', 'GCAAGGC', 'CAACAGU', 'UGAAGAG', 'AGUGGUC', 'UGAACUG', 'GGUGAGG', 'GCAUCAG', 'GUGUCUU', 'AUGAGAG', 'UCUCCCA', 'GGGAUCC', 'ACUCUGA', 'CAACCUA', 'AGUUGGC', 'UGCCACC', 'CCAUGUU', 'AGGGGGA', 'ACAACAA', 'AGGCCAC', 'CACAUGG', 'GCUGUAC', 'CACAAAU', 'CUGGCAU', 'GCUUAAG', 'GAGGGCU', 'GCAGGAG', 'CCUGCCA', 'UUGGACA', 'GGGGGAG', 'AGCAAGA', 'ACUGAAC', 'UACCUUC', 'CCCAGGU', 'UCUGCCA', 'UAGUGGG', 'AUGCUGA', 'CAUUUAU', 'UGGCGGA', 'GCUGGGG', 'ACACAAG', 'GUAGAGC', 'CCCUGAC', 'AUCCCAA', 'CACCAGG', 'CUGCGGG', 'CUAAUUU', 'GGCAGGU', 'UGCGAGG', 'CUGGCUG', 'UGCCACA', 'CCAGCGC', 'GAAACUG', 'GGACUGA', 'CCCGAGA', 'CAGGCAA', 'CUGUAUU', 'GGCCGGA', 'UUCUGUC', 'CUAAAGA', 'UUAGUGU', 'CGGUCCC', 'GCGGUGC', 'GGGAAGG', 'AGGCCCG', 'GCCCCUC', 'UCGGGCG', 'CUGACCC', 'GAGGACC', 'UGUGGAU', 'AAGGCCA', 'UCUCUUU', 'AGUGCCC', 'UGCGGAC', 'UUCCCUU', 'GCCAGGC', 'CAGACUU', 'CAGAGAU', 'GGCCUCU', 'AUGACGU', 'GCCCCAC', 'AGGACUA', 'UUAGAUU', 'AAUUCAU', 'CAAGGUG', 'AGGGCAU', 'CAAAUCU', 'UUCUGAU', 'GCCUGCC', 'GGAUGUG', 'UAACUCC', 'GAGUGAU', 'CUGAAAG', 'UAGCAAU', 'GCGGGCG', 'UUCUCUC', 'CAGGAGA', 'CUGCCAU', 'GCAGACU', 'GAUCAGG', 'CUGCAAC', 'AGAACAG', 'CUGGUAU', 'UUGCCUA', 'UAAUUUU', 'UGGACCA', 'UUGCCAU', 'UCUAGAU', 'UACCUCA', 'AUUCUGU', 'CUUCUUC', 'AGGAGGG', 'CCCUUGA', 'UCCACUU', 'AGCGGGG', 'AUGUUGG', 'UCUGGAU', 'GCGCGCC', 'CCCGGUG', 'GAGUCGG', 'CAGAUGA', 'GAGGCUC', 'GCAGAGG', 'GAGACCA', 'GAAGCCA', 'CGGGGGU', 'AUGCGCC', 'UACGGAC', 'UAUACAC', 'ACACAUA', 'UCGCUUU', 'GAAUGGU', 'CAUCCUG', 'CUGCACU', 'GAAGUGG', 'UAUUAUU', 'GUCUAUA', 'AAAGUGG', 'ACAGAGU', 'CUCAGUA', 'UCGGUAU', 'ACUCACG', 'UCUAAAU', 'GUGGACC', 'AUCCGUC', 'AUGGAGG', 'ACAUGGA', 'AACAUAA', 'UGGACGG', 'GCUUAAC', 'GCUGUAU', 'CACUACC', 'AGUUCAG', 'UCUGCCU', 'AUUUGGU', 'GACUGCC', 'CACAACA', 'ACUUUUC', 'GAGGACA', 'UUGGAUU', 'GAAGUAA', 'AAGACCG', 'UGCCAGG', 'AGAGUCU', 'UCAUAUG', 'GAGGCCC', 'CUGUGCU', 'UGGACUU', 'CCUAAAU', 'UUUGUGU', 'AUAUAUA', 'UGCAUGG', 'GGUUUGU', 'AGGGCUC', 'CCCUUCU', 'UGGGAUU', 'UGCUACU', 'CGGGGCA', 'CGGAGAC', 'AUCCACG', 'GGAAAUG', 'AUCAGUG', 'GAGAUUG', 'GGGAUGA', 'CUGAAUC', 'AUCGGAC', 'GCCAACC', 'CAGUGAC', 'GGAUAGG', 'GGAGAGU', 'CCUCCUC', 'ACCCCUA', 'CAUCCUC', 'AAUGGCA', 'AGAAGAG', 'AGACAGG', 'AAUUCUC', 'UCCUAGG', 'UCAGAUC', 'UUGGGGU', 'AAAACCC', 'AUGGUAC', 'UAGCUUA', 'GGGCUAA', 'GCUGGCU', 'ACAUAUG', 'UCCAUGC', 'AAAGUUA', 'AAAACCG', 'AGGCACA', 'AAAACUU', 'AAUAUGG', 'AGGGAAA', 'AGUUCUU', 'GAAGUAC', 'UGAAUAG', 'AUCCAGC', 'AGAGUGA', 'UGGUCAC', 'CCCCGGG', 'GGCAGUU', 'AUAGGCC', 'AGUCCCA', 'GCUGGGU', 'GCACAUG', 'UGCCAUA', 'UACCCAU', 'AAACUGC', 'AGAAAUG', 'GAAAGGG', 'UAGCACC', 'UGCAAAA', 'ACAGAUG', 'ACUCUAG', 'AUAAUAU', 'UAGAACG', 'GGUAUUG', 'GCAUACA', 'CAGCUAC', 'UGCUCUG', 'AAAUAAU', 'UGCAUGU', 'CAGUGGC', 'AGAUCAU', 'CUCCAAG', 'UCAUUUA', 'CAAGUAG', 'GCCCCAA', 'AAUGCAU', 'UAUUGUC', 'GAGUCAG', 'UAGGCCA', 'GUUUCGG', 'CGUUUGA', 'UGAGCGA', 'CGGAGAG', 'CUGAGCC', 'GCGGCCG', 'CUGCGAG', 'GGCUAGG', 'CGGUUCC', 'UCUGCUG', 'CCUCAUC', 'GUAGUGA', 'AUAUUCA', 'CGGCCCA', 'GCAUGAC', 'GGAAGAG', 'CGCCUGA', 'UGGAAGC', 'CUAGUGC', 'GGAGCAG', 'AAUACGU', 'UUAUAUC', 'UCCGCCA', 'AGGGGCU', 'GAGGCCG', 'GGGAGCG', 'GGAAAAG', 'CGGAAGG', 'UGAAGGC', 'CUGGAAU', 'CUGGAUU', 'GCUGGUC', 'GCAGGGC', 'GAGGUGG', 'AGCAGGA', 'CGGGCGC', 'GCAGUGU', 'CACUUGU', 'GUUGCCA', 'CAGAGCA', 'AGACCAU', 'GGACUAG', 'AUUACAG', 'UGGAAUA', 'GACUUCU', 'CGUAUCA', 'AAGAAUA', 'AAAGUCC', 'CUAGAAA', 'GGGCCAU', 'UUAGGUA', 'UCCACUG', 'AGCAGGG', 'ACCGACU', 'CUCACCA', 'ACCAUUA', 'UCCAGCC', 'UUGGGAU', 'CAAGUGU', 'AUGGAGU', 'UGCCUGU', 'UGGAGGG', 'CUCUUCA', 'CAAACCA', 'UCAAACC', 'GGGAAUG', 'CCGAACU', 'GGCGAUG', 'UGCAGGC', 'CUGACAU', 'GCGCCUG', 'GGGCAGG', 'GGCGCAC', 'UGGGCCC', 'AUUCAGA', 'UCGCCCU', 'UCGGGGC', 'UGGGAUG', 'CUCUGCU', 'CAUCCCC', 'GAAAGGU', 'CUGACAC', 'GGGAGAG', 'CUAUUCC', 'AGGGGGU', 'AACCCUG', 'CAGUGUC', 'UGAGGGG', 'CCUUCCC', 'GGCCUGU', 'CAGCUCC', 'UGGGGUG', 'UUCUGCC', 'UUGUUCU', 'GUUUGGG', 'GUCUUCU', 'UGGGUGC', 'CGGCUCU', 'GUGGGGU', 'CCUGGGU', 'GCCGCUC', 'GGAUGAC', 'CGGGAGA', 'AGCCGCC', 'GGGGUGU', 'CCUGCCU', 'GUGGGUG', 'CCUGUCC', 'CGGGCCU', 'UCCUCCC', 'AGGGAAC', 'AACUCAC', 'UGGGGGU', 'GGGGGUG', 'CCCUGCC', 'ACCAGGG', 'GGUCUGU', 'AGGGUAG', 'GUUGUCA', 'GGGUGGG', 'AGGGAUG', 'ACACUGG', 'AGAGAGG', 'CUCAUUC', 'UGUGGGU', 'GACCUUU', 'AGGGAGA', 'CACUGUC', 'CUGAGAG', 'CCUACGC', 'GGGGCCA', 'GGCUGCU', 'UCCCCGG', 'CUCUGGU', 'UGAGGCG', 'CACCUGG', 'CGCAGAC', 'CACGUGC', 'ACACAGG', 'AAAGGCC', 'GGUGGGU', 'AGCCCCU', 'GAGAAGG', 'UGGCGGC', 'UCGGGAG', 'AAACCCC', 'GGGUGUA', 'UGCUCCU', 'CUGUCAC', 'CGUGUCC', 'CUGGGUG', 'AACCACC', 'CGGGGAG', 'CCACUCU', 'GCCUCCC', 'AGCCCUG', 'UCCUCUG', 'GGGCCGG', 'AGGGGUG', 'ACCUUUG', 'AGGGGAA', 'CCGGGGC', 'CUCACCC', 'CAUCGCC', 'CGGUGGG', 'GACGCCC', 'CUCAGCU', 'UCGCCAC', 'UAGGGGC', 'GGCGCCC', 'GACCUCG', 'GCCUCCU', 'UAAGCAG', 'GUGGGUU', 'CCCCAAC', 'UCACUCU', 'GGGGGGA', 'CCCCUCG', 'UGUGGGG', 'AAGCUCU', 'CAGGGGG', 'UACCCCC', 'GGGAGGU', 'UAGGUGA', 'ACCUCUC', 'GGUCAGA', 'CCCCUGC', 'UAGGUGG', 'UCACCCC', 'CCCUCGC', 'GAGGGUG', 'GCACCUG', 'AGGGGGC', 'GUAGGCA', 'UGAGCCA', 'UGUGACC', 'GGCAAGG', 'UUCUCUU', 'UGGGGAC', 'GCCUGUG', 'CGCUCUU', 'ACCUUGG', 'CCCAAGG', 'CUCGCAU', 'AGGUGGC', 'GGCUUCU', 'AAGGACC', 'CUCUCUG', 'UGUGUGA', 'UGUCUCU', 'AGCCUCU', 'GCGGCAG', 'UGCGUGG', 'GACCUCU', 'GGCUCUA', 'CAGGGUU', 'UAGGGGA', 'GGGGAGG', 'CGCUGAC', 'CAAUAGG', 'UCCCCUC', 'GGGAGCC', 'GGAAGCA', 'UCUGCUC', 'GGGCUGC', 'CAAGGAA', 'GUCUUUC', 'GACUAAC', 'GUAGAGA', 'CCCUUUU', 'UGUGGAA', 'UGAGGGA', 'AUGUCCC', 'GGGGGUA', 'AAAGCAC', 'CCCUUGU', 'UGCCUCC', 'CUUCACU', 'AGUCCUG', 'UGGGUUU', 'CCCCCGG', 'CCCAGGA', 'AGGGUAC', 'CCUUGCA', 'UGGCUGG', 'UGGUCUC', 'UCUUUGU', 'CUCUCCU', 'GGGGGCU', 'GACCCCU', 'CAGAGGA', 'GCUCAUG', 'CCAGCCU', 'UGCGGAA', 'CCGGCCG', 'GGCCCUU', 'CCUGGGG', 'GCGUGGG', 'GUUCAUU', 'AGCUCAG', 'GCGUUUC', 'ACAGCCC', 'UGGGGAU', 'GACUGAG', 'AGCCAGC', 'GACCCCC', 'CUGGGUA', 'GGACCUC', 'AGACGUG', 'UGAAGGG', 'CACCCUC', 'UAGAGGC', 'GUGUGUG', 'UCUCCCU', 'CUGGCAG', 'UCCUUCU', 'UGAGUAG', 'GCCGCGC', 'CUCAUCC', 'AUGGGAG', 'AGCACCC', 'CCAUGCC', 'AGAGGGA', 'UCUCUCU', 'UGGAGCU', 'UUCUUCC', 'GGGCCGA', 'GGGAGAA', 'UGGCCUC', 'GUCACCC', 'GGUGGAG', 'CGCCUUC', 'GGGGUAA', 'UCCUCUU', 'ACAAGUC', 'GCUGCCU', 'UCCCUAU', 'CCAUCAC', 'UUUGCUU', 'CCGCAGG', 'GCCCUUC', 'CCCCUCC', 'AGGAGAU', 'UCUGUCU', 'CUGUGCC', 'AUGGGGU', 'CACUGCC', 'UAAGGGA', 'CCCUCUC', 'GGAGGAU', 'UGCCUGC', 'AGGGCCA', 'GUCUCUC', 'GGGAGGA', 'GCUCCCU', 'UGUGGCC', 'CCCGCCC', 'UGGGGGG', 'AAGCCUC', 'CUCUCUC', 'UCCUCUC', 'CGGGCAG', 'GCAUCAC', 'CCAGGGA', 'CUCCCUG', 'CUGUGGA', 'UGGCAGG', 'AUCCAUC', 'UACAGGC', 'UUCCUGU', 'CUGGUCC', 'ACCAUGG', 'UCAUGAA', 'GGAGGAC', 'CUGGGGU', 'UGCAGCC', 'GCUCAAU', 'UGAACUA', 'AGAUCUU', 'UCAACAA', 'UUCUAUG', 'GCUGAGG', 'AGGGCCC', 'CUGAGGU', 'GAAGGGA', 'UUAGACU', 'UGCACUC', 'GGGGUCG', 'AUAGCUC', 'GAAGCGC', 'UGAAGCC', 'UGUCCCA', 'AGGAGUG', 'UGCAGAC', 'UUGCAGU', 'UCAUGUA', 'AAACUAG', 'AGGGACA', 'AGCGGAG', 'GUGGAGG', 'UACCCUC', 'ACAAUUG', 'UUUGGAC', 'ACCUGGG', 'AUGUAGU', 'GAGGUGA', 'UGGUGAG', 'UUUAAGG', 'GUGACCC', 'GGCUGUG', 'UCCUAGU', 'GCCCUGA', 'UCCCAGC', 'CUGGUGU', 'GGCGAUU', 'UUUGAGC', 'GUGGAUU', 'UGGCUCU', 'CAUGAAG', 'AGUGAUU', 'CAAAAUC', 'GCACACU', 'GUAGGAA', 'AAUGUGA', 'GUUUGUU', 'GAUGGUU', 'UGUGAUU', 'GGUGGAC', 'UAUGGCG', 'GCUGAUG', 'AUAUGGA', 'GCUGAGU', 'GUCUAGG', 'AAGGACA', 'UUGAGUC', 'AGGACUU', 'AAUACUA', 'GCUAGUC', 'AAGACUU', 'CUCGGUA', 'ACACACA', 'AGGGAAG', 'ACUUACC', 'GGGCGCC', 'CAGGUCC', 'GGUCGCC', 'CCAGAAA', 'UAUCCUC', 'CAUGUGU', 'CGCUUCG', 'CGGCCCC', 'GAAUUCU', 'UGAUGGA', 'UGUAUUC', 'GCCCCGG', 'CCUUGAC', 'AUAGAUC', 'GGCUCCC', 'CCCGGAG', 'UCUCGGA', 'GUGUGGG', 'CCUCGCC', 'GAAGAAU', 'ACAAUGA', 'AGGAUGC', 'GACUAUG', 'AAAGGGG', 'AAGCAAA', 'UCAGGUC', 'CGGCUGU', 'CACGUCU', 'UUCCUCU', 'CUCUAGG', 'GGUGCCA', 'AGAGAUC', 'GCUGGCA', 'AGGUGGA', 'AGUGGAG', 'UAGGCUU', 'AAGUGGA', 'UCUGAGG', 'AAGGAAC', 'UGCCACG', 'GUCACAG', 'UAUUCAU', 'AGGAAAU', 'CUCUCAC', 'ACUUGGG', 'AGAGGUU', 'UCAGGGA', 'AUGUACU', 'GGGUUGU', 'UUGGAGA', 'AUUACUG', 'AAAGGUU', 'AAAAAGU']
         
## mature_uniquseeds_logoplt for hsa is the following:
![image](https://user-images.githubusercontent.com/130226484/234348871-5bdc9de7-aa3b-428b-b810-ae4943b1d58b.png)
## the calculation of GC frequency for "hsa":
          A         C         G         T
0  0.249284  0.236867  0.282235  0.231614
1  0.244986  0.227794  0.295129  0.232092
2  0.205826  0.247373  0.346705  0.200096
3  0.218243  0.225883  0.323782  0.232092
4  0.225883  0.257402  0.284623  0.232092
5  0.239733  0.244031  0.278415  0.237822
6  0.217287  0.267431  0.256447  0.258835

## for "rno" the out put are as the following:
## the out put of logger method for GC calc of "rno":
2023-04-26 12:05:31,673 - root - INFO - +******************************************************************************+
2023-04-26 12:05:31,674 - root - INFO - project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\mature__20230426_120531__3697017c62772d71b4445895933026a4.log>
2023-04-26 12:05:31,674 - root - INFO - +******************************************************************************+
2023-04-26 12:05:31,675 - root - DEBUG - debug mode is on
average GC % = <67.4683526671524>
## The number of sequences are:  764
fasta file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa>
speciesCode is <rno>
+******************************************************************************+
project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\mature__20230425_195834__3697017c62772d71b4445895933026a4.log>
+******************************************************************************+
debug mode is on
## the number of mature_uniqueseeds for "rno" are:635
write unique seed sequences to fasta file
output fasta file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature__uniqseeds.fa>
['AGCAGCA', 'AACAUGA', 'GGUGGUC', 'ACAUUAC', 'CUCUGAC', 'AGUGCAA', 'GCAUCCC', 'CACUGCC', 'CUAGUAG', 'UUAUUGA', 'GGGGCAG', 'CUCUGGG', 'CUUGAGG', 'GAGGUAG', 'UAUACGA', 'GGGGGCA', 'UGGCCCU', 'GAGGUUU', 'ACACACC', 'CAAAGCA', 'GUCUUGU', 'CCCCUGG', 'AGUGGUU', 'ACCACAG', 'CAAGAGC', 'CACCCUU', 'CUGGAUU', 'GGCGUCA', 'UCAGCUC', 'AAGUUCU', 'CAGUGCA', 'ACAAUAU', 'CCAGCAU', 'CCCUGUC', 'GAGCGCC', 'UAUAAAG', 'CCGUCUC', 'CGGUCGA', 'GGGGUGC', 'CUCACAC', 'CUCCCUC', 'CAGGCUC', 'GAUCUAG', 'GCUGACC', 'CCUGAAC', 'GUCUGCC', 'GUCCCUC', 'AGCCCUG', 'UUUUUGC', 'AGCCCUU', 'AAAGUGC', 'CUGCAUU', 'UCACAAA', 'GGAAGAC', 'CAACAAA', 'CCCUGAG', 'GUCAAGA', 'GAGUAGU', 'AUGGCUU', 'UGUAGGG', 'CGAGGAG', 'UAGACUG', 'CGGUUAU', 'ACAGUAC', 'UAUACAA', 'UGUACAG', 'UGUACAA', 'UAUACGG', 'UAUACAG', 'UGCGCAA', 'AACAAGU', 'CUUUGGU', 'UAAAGCU', 'ACCCUGU', 'AAAUUCG', 'CCUGUAG', 'CAGAUUC', 'GAAUCAU', 'CCAAUAU', 'CUGCAGU', 'AAGGUGC', 'CUGCCCU', 'GUUUUGC', 'GUGCAAA', 'CGUUUUG', 'AGCUUAU', 'AACAGCA', 'GUUCUUC', 'AGCUGCC', 'GGGUUCC', 'UCACAUU', 'GGUUCCU', 'UGCCUAC', 'GGCUCAG', 'GGCGGAG', 'AUUGCAC', 'UCAAGUA', 'CUAUUCU', 'CUGUUCU', 'GAGCUUA', 'UCACAGU', 'GGGCUUA', 'AGGAGCU', 'ACUAGAU', 'UGGUUUC', 'AGCACCA', 'CUGAUUU', 'UUCAUAU', 'GACCGAU', 'GUAAACA', 'UGGGAGA', 'UUUCAGU', 'UGGGAUG', 'GGCAAGA', 'GCUAUGC', 'CAAUUUA', 'UGCAUUG', 'AAUGUUU', 'GGCAGUG', 'AUCACUA', 'AUCAGCA', 'GGUUGGG', 'GGUGGGG', 'CUGCUGA', 'UUGGCAC', 'AAUCAUG', 'ACCCGUA', 'AAGCUCG', 'AAGCUUG', 'CAGUUAU', 'GCUUCUU', 'GCAGCAU', 'CGCACUG', 'GGAGUGU', 'ACGCCAU', 'GUGUUCA', 'AAGGCAC', 'CAGGUGA', 'CGGGUUA', 'CAAGUCA', 'AUUAUUA', 'CGUACCG', 'UGAAGCU', 'CGGAUCC', 'GGGGCCG', 'CACAGUG', 'CUCUUUU', 'CUCUUUC', 'CCGUGGC', 'AACAGUC', 'GCUGGUA', 'UUGGUCC', 'GUGACUG', 'UGUGGGC', 'GUAGGGA', 'CUCCAUU', 'AUCAUCG', 'CGGGUAU', 'UAUUGCU', 'GCUGGUG', 'CUAUUUC', 'GGCUACU', 'CUACAGU', 'GGAGACG', 'CCAUCUU', 'AACACUG', 'AUAAAGU', 'GUAGUGU', 'GUGCAGU', 'GAGAUGA', 'GAUAUCA', 'ACAGUAU', 'UCCAGUU', 'GAUUCCU', 'GAGAACU', 'CCUGUGA', 'CUCCCAA', 'UGGUACA', 'GGUUCUG', 'UCAUUUU', 'UGCAUAG', 'AGGUUAU', 'AUCAUAC', 'ACAUUCA', 'CCAUCGA', 'CCACCAA', 'UCACUGA', 'UCACUGG', 'AUGGCAC', 'GAAUUAC', 'GGACGGA', 'GGAGAGA', 'UUCCUCU', 'AAAGAAU', 'CCCAAAG', 'GGCUACA', 'CGUGUCU', 'GAUAUGU', 'CUAUAUA', 'AACGGAA', 'CUGCACU', 'UGACCUA', 'UGCCAGU', 'GGGUCUU', 'ACUGGCC', 'GUAACAG', 'CAGUGGG', 'CAAUAUU', 'AGGUAGU', 'CGGCAAC', 'CCAGUGU', 'CAGUAGU', 'GUCUUAC', 'AAUACUG', 'AUCUUAC', 'GUGGUUC', 'UGAAAUG', 'UCCCUUU', 'CUGGGAA', 'CCUUCAU', 'CAUGCUU', 'GGAAUGU', 'AGCUUUU', 'UAAGACG', 'GCCACUG', 'UGUGCGU', 'GCAAGGA', 'CCUUGGC', 'GAGUUGU', 'CAGCAGG', 'AAUCUCA', 'ACAGUGG', 'ACUGCAU', 'UCAGUUC', 'UGUGCUU', 'AUGGUUC', 'AACAUGG', 'GAUUGUC', 'GAGUUGC', 'GAAUUGU', 'CCUGGCA', 'GCUACAU', 'GCUCAGU', 'GUGUAUU', 'GUCAGUU', 'CUCAAAC', 'AUCAAAG', 'AAGUGCU', 'AGUGCCG', 'GGGCCCC', 'AGGGUUG', 'UGUAUGU', 'GCAGAGG', 'GGAACUA', 'GGUUUAC', 'AUGUGGG', 'UGAAGAG', 'AUGCAAG', 'CCUUCUC', 'AAAGCUG', 'CGACAGC', 'GCCUCAU', 'UCAACAG', 'ACAUCCU', 'UGCAUAU', 'AGCUAAC', 'UUUGCGA', 'AUUGGGA', 'AGGGACU', 'AAUGCCC', 'GUCUUGC', 'ACGGUGA', 'UCAUGAU', 'AACCGUU', 'UGGUAAU', 'GUCGUAU', 'UGACAUC', 'AGAUCAG', 'CACAGCA', 'UUCUCCU', 'AUCCUUU', 'AUGCACC', 'UAUCAGA', 'CCCCAGG', 'AAGUCAC', 'AAUGGUG', 'AGACGGG', 'ACUCCUC', 'AGGUCAC', 'CCUGCUG', 'GAGGUUG', 'GAAUCAC', 'GGUCGAC', 'UUCACCU', 'CACAUAC', 'CUGGUCA', 'AAGGGUC', 'GGUCAGA', 'AGUUGCC', 'AACAUUC', 'GAGAAAU', 'AUACAAG', 'AGGGAUU', 'GUGGCGA', 'UCGGGGA', 'GUGACAG', 'UGUACAU', 'GAAGGUC', 'GGUAGAC', 'UAUGUAA', 'GGUUGUC', 'GAAACAU', 'UGGAUAU', 'ACAUAGA', 'UCAUAGA', 'GUAGAUU', 'UCGUAGA', 'GCGAGGU', 'GUGGUUA', 'AUCGUAC', 'AAGUUGU', 'AUCAUUC', 'GAGGCUG', 'AUACACG', 'GGUUACC', 'AUGUUGC', 'GAUCGAC', 'AUAAUAC', 'UAUAAUA', 'UUAGCAC', 'GGGUGGA', 'AGCAGCG', 'GAGUAUU', 'UCCUGAC', 'CUGGACU', 'GGAGCCA', 'UCAACAC', 'UAAGACU', 'ACAUCAC', 'UGGCUGG', 'AUUCAUU', 'CAAACCA', 'AUGUGUG', 'UACAUAC', 'UAUACAU', 'GUGAUGU', 'AUACAUG', 'GUUCAGA', 'AAAGACA', 'AUUCAGA', 'GACUGGC', 'AGGUUAC', 'GAACUAU', 'GGCCCCA', 'UGCCCUG', 'UAGAGGA', 'CAUGACA', 'GAGGCUU', 'UUAUGGC', 'ACUCAGA', 'ACUCCAU', 'AACUGUG', 'GCUGAGA', 'AACUGCA', 'UGUGCGG', 'CACCGGG', 'UUGGCAA', 'AUCCCUU', 'UCCCACA', 'ACUAAAU', 'CAACAAC', 'CGACGAG', 'UUGUUCG', 'UGGUUGA', 'AUGUAGU', 'UUCCUAG', 'AUAUAAC', 'AGUAGAC', 'AUGUAAC', 'GAGGGGC', 'GCUCGGU', 'AUGACAC', 'UCGGGAA', 'GCUCGAC', 'UUGAACC', 'AUGUGCC', 'CAGUCCA', 'ACCUAAU', 'GAUAGAC', 'ACGUAGU', 'GAAAGGU', 'AACAAAC', 'AUCCUUG', 'AUGCCUU', 'CUCCCAC', 'CGGUGAU', 'ACGUCAU', 'CCGGUUC', 'CACACAG', 'UCACAGC', 'CCGGGAC', 'CACUGAG', 'ACAGCUC', 'AACUAGA', 'ACUCACA', 'AAAGCCA', 'AAAGACG', 'GGUUGAC', 'UUGUGAC', 'CCCUCAG', 'GGCUCUG', 'GCACCAC', 'UGGGCCU', 'CAGGAAC', 'AGACUGA', 'GGGACGG', 'CCAGAUA', 'UGAAAGG', 'CAACCCU', 'AUGGCGC', 'CUAGGGA', 'AAAUCAA', 'GCGACCC', 'UCCUAUG', 'GAGGUAU', 'CAUGGAU', 'AACCUGG', 'UGAGGAC', 'AUGUGUA', 'GACCCUG', 'UUGUUAA', 'UUCUGCA', 'UGUAUAA', 'ACAGUUG', 'ACCUGUU', 'UUGUGUC', 'GGGGUCC', 'UGCUGAC', 'UGUCUGU', 'UGAAACA', 'GUCACUC', 'CAGUAAC', 'CGGAGAG', 'GGUGCGG', 'GUAUGCC', 'AUACCUC', 'GGAUUUC', 'CCAGGAG', 'CAGUUAC', 'ACUCAGU', 'GAACAGC', 'GUGCCGC', 'UCAAAAU', 'CUCAAAU', 'AGUGCUA', 'AUUUAGA', 'GAUCAGU', 'CACUUCA', 'UUGGUAC', 'GGUGCUG', 'GACACCU', 'GUGCUCA', 'AGGAGGC', 'GCGGGCA', 'GCUGCAG', 'CAGAGUG', 'UCUCGGG', 'GGACCCU', 'GUAUUAC', 'CCCUCCC', 'AGGCUCU', 'AGGAGUC', 'ACUCCUG', 'GUGGUCC', 'UGAACUG', 'UAACCCG', 'GAGCACC', 'AGAACUU', 'CUGGCUG', 'AGCACUG', 'UCUUGCC', 'CGAGCCC', 'CUUCACC', 'GAAAUUC', 'GUCAGGC', 'GUAUAAC', 'GGCUGCA', 'AGUGAAU', 'AUAUCAG', 'UAUGCCA', 'AACGAAA', 'CCGAUUU', 'CCUUCUG', 'CACAGGA', 'GAUUUCA', 'AAGUGCA', 'UGUGUGU', 'CAUAGAA', 'CUGUGGA', 'GACAGAC', 'UGUAGUA', 'AAAUCCU', 'CUGUGUC', 'AUCCAGG', 'UGGGGCA', 'UCGGUUA', 'UACAGCU', 'GGGAAAU', 'AUGGUUA', 'CAGCACC', 'UGCCUAA', 'ACCAAAA', 'AAUCUCU', 'ACACUUA', 'GAGGACA', 'GUACAAU', 'ACACACU', 'AUACCAG', 'UACACUU', 'GGUGUGA', 'GAGGCUC', 'UCCUAAC', 'UGGUAAA', 'CAGCCGC', 'CUGCUGG', 'CACACCA', 'CUGUCCC', 'GCUCAUU', 'ACAGAAC', 'CGAAUAU', 'CGCAUGU', 'GCUUGUA', 'GACUAGG', 'GGGGUUC', 'ACUAUGC', 'AGCUGGU', 'GCAUGAA', 'CUGACUU', 'GGAGGAG', 'CCCGCAU', 'CUAUACA', 'UCACAAG', 'GAACGGC', 'UACUAGA', 'GAUGUCC', 'ACACUGA', 'UGAGUGG', 'GCCCACC', 'ACUAUAC', 'CACAAGU', 'AGGAAAC', 'GCCCUUU', 'AGCACAA', 'GUUAUAG', 'CGCCCUU', 'ACACCAU', 'AGUGUUG', 'GGCCUCC', 'UUCAACC', 'CCAGGGC', 'ACACCGC', 'ACAGCAA', 'AAGAUGC', 'UGGAUUU', 'CUGAGUC', 'GCGGCCG', 'AGGUAGA', 'GACCUAC', 'UGGAGUU', 'GUGAAUU', 'AUGGAGG', 'CUGAUCG', 'ACUCCAG', 'GAUUGAC', 'CACCUCC', 'ACGUUGG', 'CUAGGGC', 'UUAGGGU', 'GUGGACU', 'AUACACA', 'GGCACUG', 'AGGCUAG', 'UCUAGCC', 'CUGGACA', 'UCUACUC', 'AUGUCUG', 'UGAACCU', 'UGCCUGG', 'CAUUCUC', 'GUCUGGG', 'GAACAUG', 'ACUGCAG', 'UCAUUCG', 'GUCUGUA', 'GGAACUC', 'CAGUAGG', 'AAGUCAG', 'UGGUCUG', 'GGACUGU', 'GGCCUGC', 'AUGUGAC', 'AUAAUGA', 'UUUGGUG', 'GCCCCCU', 'AAGGUCU', 'GCAGGGA', 'AGGCUGG', 'GUGGCUC', 'UAAUGCU', 'UCCUACC', 'CCCGUCC', 'GGGAACG', 'UUCAUGU', 'CUGGCUC', 'GGACGGG', 'GGCCAUA', 'UGGAGGU', 'GGUGGGU', 'GGGUUUU', 'GUCCAGG', 'GGCUGGC', 'GGCUGUU', 'UUGCCAC', 'UGAAGGU', 'UCUGCCA', 'UUGAAGG', 'UUGGGGA', 'ACUGUUU', 'CAGUCUC', 'GUGGGGC', 'CUGAUGU', 'CGGGGCU', 'CGUCCUG', 'CCUGUAC', 'GAGGCAG']
## mature_uniquseeds_logoplt for "rno" is the following:
![image](https://user-images.githubusercontent.com/130226484/234350294-8b683253-f7e9-4361-97a8-6df26e1612bd.png)
## the calculation of GC frequency for "rno":
 A         C         G         T
0  0.277165  0.223622  0.288189  0.211024
1  0.251969  0.220472  0.270866  0.256693
2  0.242520  0.237795  0.307087  0.212598
3  0.244094  0.223622  0.264567  0.267717
4  0.226772  0.262992  0.231496  0.278740
5  0.277165  0.215748  0.255118  0.251969
6  0.222047  0.280315  0.222047  0.275591

## for "mmu" the out put are as the following:
## the out put of logger method for GC calc of "mmu":
2023-04-26 12:09:42,152 - root - INFO - +******************************************************************************+
2023-04-26 12:09:42,153 - root - INFO - project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\mature__20230426_120942__3697017c62772d71b4445895933026a4.log>
2023-04-26 12:09:42,153 - root - INFO - +******************************************************************************+
2023-04-26 12:09:42,153 - root - DEBUG - debug mode is on
average GC % = <70.50310821391368>
The number of sequences are:  1978
fasta file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa>
speciesCode is <mmu>
+******************************************************************************+
project log file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\logfiles\mature__20230425_201742__3697017c62772d71b4445895933026a4.log>
+******************************************************************************+
debug mode is on
## the number of mature_uniqueseeds for "mmu" is: 
found <1588> unique seed sequences
write unique seed sequences to fasta file
output fasta file is <C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature__uniqseeds.fa>
['GAGGUAG', 'CUGUACA', 'UGCGCAA', 'CAUACUU', 'GGAAUGU', 'AGCAGCA', 'GAAUCAU', 'GGUUCCU', 'UCACAUU', 'GAGCUUA', 'UCACAGU', 'CUGGUUU', 'AGCACCA', 'GUAAACA', 'UUUCAGU', 'UGGGAUG', 'ACCCGUA', 'AAGCUCG', 'CAGUUAU', 'ACAGUAC', 'GUGUUCA', 'AAGGCAC', 'CCCUGAG', 'CAGGUGA', 'CAAGUCA', 'AUUAUUA', 'CGUACCG', 'UGAAGCU', 'CGGAUCC', 'GGGGCCG', 'CACAGUG', 'CUCUUUU', 'AGUGCAA', 'CUUUGGU', 'UAAAGCU', 'ACCGUGG', 'AACAGUC', 'CUGGUAA', 'UUGGUCC', 'GUGACUG', 'UGUGGGC', 'AUGGCUU', 'AUAGGGA', 'CUCCAUU', 'UCAUCGU', 'CGGGUAU', 'UAUUGCU', 'GCUGGUG', 'CUAUUUC', 'AGUGGUU', 'ACCACAG', 'AUCUUCC', 'AACACUG', 'AUAAAGU', 'GUAGUGU', 'GAUAUCA', 'ACAGUAU', 'UCCAGUU', 'UUCCUGG', 'GAGAACU', 'CUGUGAA', 'CUGGCUC', 'AGGGAGG', 'CUCCCAA', 'UGGUACA', 'CGAGGAG', 'UAGACUG', 'AGGUUCU', 'CAGUGCA', 'UCAUUUU', 'UGCAUAG', 'AGGUUAU', 'AUCAUAC', 'UAAUGCU', 'UCCUACC', 'ACCCUGU', 'AGAUUCG', 'UUUUUGC', 'AGCCCUU', 'ACAUUCA', 'CCACCGA', 'UUGGCAA', 'UGGUUCU', 'AUGGCAC', 'UGAAUUA', 'CUUAUCA', 'GGACGGA', 'GGAGAGA', 'GGGGCUG', 'AAAGAAU', 'CCCUAAG', 'GGCUACA', 'CGUGUCU', 'AUCCCUU', 'UCCCACA', 'UGCCUAC', 'GGCUCAG', 'GAUAUGU', 'CUAUAUA', 'AACGGAA', 'CUGCACU', 'GGGUCUU', 'ACUGGCC', 'GUAACAG', 'CAGUGGA', 'CAAUAUU', 'CCAGUGU', 'CAGUAGU', 'AUCUUAC', 'AAUACUG', 'ACUCAGU', 'GAACAGU', 'UCCUAUG', 'GAGGUAU', 'GUGGUUC', 'UGAAAUG', 'UCCCUUU', 'CUGGGAA', 'CCUUCAU', 'AUUUCAG', 'CAUGCUU', 'CUUCUCC', 'GGAGUGU', 'AACGCCA', 'GUGCAGU', 'GAGAUGA', 'CUCAAAC', 'AAGUGCC', 'AUCAAAG', 'AAGUGCU', 'GUGCCGC', 'CUCAAAA', 'CUCAAAU', 'GGGCCCC', 'AGGGUUG', 'UGUAUGU', 'GCAGAGG', 'AGGAACU', 'GGUUUAC', 'AUGUGGG', 'UGAAGAG', 'AUGCAAG', 'CUCUGAC', 'CUUAAAC', 'GGCAGUG', 'AUCACUA', 'UAUACGA', 'AAAGUGC', 'CUGCAGU', 'CGCACUG', 'CUCUUUC', 'GUUUUGC', 'GUGCAAA', 'UGGGAGA', 'AAGUUCU', 'UGACCUA', 'UGCCAAU', 'AGGUAGU', 'AACGACA', 'CGGCAAC', 'AGCUUUU', 'UAAGACG', 'UAUACAA', 'UGUACAG', 'UGUACAA', 'UAUACGG', 'UAUACAG', 'AGGCCAU', 'CAGUAUU', 'CCAAUAU', 'AAGGUGC', 'CUGCCCU', 'CUGCAUU', 'AGCUUAU', 'AACAGCA', 'GUUCUUC', 'AGCUGCC', 'GGGUUCC', 'UCAAGUA', 'CUAUUCU', 'CUGUUCU', 'CUGAUUU', 'GACCGAU', 'GGGCUUA', 'GGCAAGA', 'GCUAUGC', 'GGUGGGG', 'AUUGCAC', 'CUGCUGA', 'UUGGCAC', 'AAUCAUG', 'AUCAGCA', 'GCUUCUU', 'GCAGCAU', 'AACAUGA', 'GGUGGUC', 'ACAUUAC', 'GCAUCCC', 'CACUGCC', 'CUAGUAG', 'UUAUUGA', 'GGGGCAG', 'CUCUGGG', 'GGGGGCA', 'UGGCCCU', 'GAGGUUU', 'ACACACC', 'CAAAGCA', 'UAGGUAU', 'CCCCUGG', 'GGCGUCA', 'CAGCUCC', 'ACAAUAU', 'CCAGCAU', 'CCCUGUC', 'GAGCGCC', 'UAUAAAG', 'CCGUCUC', 'GGUCGGC', 'CGGUCGA', 'GGGGUGC', 'CUCACAC', 'GUCAGGC', 'GAUCUAG', 'CUGACCC', 'CUGAACU', 'GUCUGCC', 'GGCAGGG', 'AAGUGCA', 'UCACAAA', 'GUCAAGA', 'UGUAGGG', 'CGGUUAU', 'UACAGUA', 'AAAUUCG', 'AGUUUUG', 'GGCGGAG', 'AGGAGCU', 'ACUAGAU', 'AAUUUAG', 'CAAGCUU', 'CUACAGU', 'GGAGACG', 'GUCUUAC', 'GCCACUG', 'UGUGCGU', 'CCUUGGC', 'CCAUCGA', 'GCCUGUC', 'CAGCAGG', 'AAUCUCA', 'ACAGUGG', 'UGUGCUU', 'AACAUGG', 'AUGGUUC', 'GAUUGUC', 'GAGUUGC', 'GUGUAUU', 'GUCAGUU', 'CCUUCUC', 'AAAGCUG', 'UGCAUUG', 'AAUGUUU', 'CAAGGAC', 'CCUGGCA', 'GCUACAU', 'UCAGUAG', 'AAGUCAC', 'AAUGGUG', 'UGGUUUC', 'GUAGGGA', 'GGUUGGG', 'GGCUACU', 'UCACUGA', 'CGGGUUA', 'GGAAGAC', 'AACAAAU', 'AACAAGU', 'ACUGCAU', 'AUCAGUU', 'CAGUGGG', 'GAAUUGU', 'UAUCAGA', 'CCCCCAG', 'AUCCUUG', 'AGGUGGA', 'GGGACUU', 'AAUGCCC', 'CGACGAG', 'UUGUUCG', 'GUAGAUU', 'UCGUAGA', 'GAGGUUG', 'UCACACA', 'UCCUGAC', 'CUGGACU', 'GGUAGAC', 'AUGUAAC', 'UGGUUGA', 'AUGUAGU', 'GCGAGGU', 'AUACAAG', 'AAGUUGU', 'CAUUCAC', 'GAUCAGA', 'CACAGCA', 'CAAGAGC', 'UUUUCAU', 'CUGGUCA', 'CUGUCAU', 'UUCCUAG', 'CGACAGC', 'GGUUACC', 'AAUGUUG', 'GGUUGUC', 'AUAUAAC', 'UGGAUAU', 'UCAUAGA', 'AGUAGAC', 'GGUCGAC', 'UCACCUG', 'AGGUCAC', 'CCUGCUG', 'AUGACAC', 'UCGGGAA', 'GUCUUGC', 'AGGUCGU', 'ACGGUGA', 'UCAUGAU', 'CUCGACU', 'UUGAACC', 'AACAUCC', 'UGCAUAU', 'UCUUACC', 'AGCUAAC', 'UUUGCGA', 'UUGGGAA', 'AACCGUU', 'GUUUGCA', 'CAGUCUC', 'ACCUAAU', 'GAUAGAC', 'AUUUAGA', 'AUCAGGG', 'AUGUGUG', 'AUACAUA', 'AUAUACA', 'ACUGAUG', 'AUGACUG', 'UCUUGGA', 'ACCAGUA', 'ACGUAGU', 'GAAAGGU', 'AUGCCUU', 'CUCCCAC', 'GUCAUAU', 'AGACGGG', 'CACUCCU', 'CAGGCUC', 'GAGGCUG', 'GUCAUAC', 'CCUGUAC', 'UGGUGGC', 'AAGGGUC', 'GGUCAGA', 'AGUUGCC', 'AACAUUC', 'GAGAAAU', 'AGGGAUU', 'GGCGAAC', 'UCGGGGA', 'GUGACAG', 'ACUUGAG', 'UUGGUAC', 'CUGUUGC', 'GAAACAU', 'ACAUAGA', 'GGUUAUC', 'AUCGUAC', 'GAUCGAC', 'AUAAUAC', 'UUGGGGA', 'AGCAGCG', 'AGUAUUG', 'CUUUAAC', 'AGUGCUU', 'UGAGGAC', 'CCCACCU', 'CCCGUCC', 'GGGAACG', 'GGUGCGG', 'UGUAUGC', 'GCGGGGC', 'UGUUGCC', 'UAUAAUA', 'GUUGUAU', 'AAUCUCU', 'CACUUAC', 'UUGUGUC', 'CAUCACG', 'GGUUGAC', 'UUGUGAC', 'GGUCCUC', 'AAAUCUU', 'AAAUCAA', 'CGACCCA', 'GGAGGGA', 'CGCCCUU', 'GGAAGCC', 'CCGGUUC', 'UAAGUGU', 'GUCACUC', 'AUGGAGG', 'CUGAUCG', 'GGGGCCU', 'CCAGGAG', 'GGUGCUG', 'GACACCU', 'GCACCAC', 'GUGGGCC', 'GGGCUGG', 'CUGGAGG', 'CAGUAAC', 'CGGAGAG', 'CACACAG', 'UCCCUGA', 'UUCCUCA', 'UCAUUCG', 'AUCUGGG', 'GUGCUCA', 'GGAGGCC', 'UUCCUCU', 'CUCUAGG', 'CAGCCAC', 'UCCUGUC', 'CAGCUUC', 'AGGCUAG', 'GGCUCAU', 'CUGCCUA', 'CAGCUGG', 'GUUGUGU', 'CAUAACA', 'GCGGGCA', 'GCUGCAG', 'CAGAGUG', 'GGUUGCC', 'GAGUAUU', 'UCACAGC', 'CCGGGGC', 'CCCUCAG', 'GGCUCUG', 'CACUGAG', 'ACAGCUC', 'GUCUGGU', 'CCAGAUA', 'UGAAAGG', 'UCAGUGA', 'AAGCCAG', 'UCUCGGU', 'AAACCAC', 'GAGGGGC', 'GCUCGGU', 'GACUGUG', 'GCAAGGU', 'AAGUUGC', 'AACAAAC', 'GGCAUCU', 'AGCCUCG', 'UGCAGUC', 'AGUUGCU', 'CUGCUGU', 'GUUUUCC', 'UUGCUUC', 'UCUCGGC', 'UAUCCUG', 'CGCAGGC', 'AAGGCUA', 'UUCCUGA', 'UCUCUUU', 'AGCCACA', 'CAGCUUU', 'UGAAAAU', 'CCCUAGG', 'GUUUUGU', 'UAAGUGC', 'UAUACAU', 'UAGUUGU', 'ACACACA', 'GAUUGGG', 'UAAGACU', 'AACAUCA', 'CGUGUGC', 'AUGUGCC', 'CAGUCCA', 'GUGGGGA', 'UUAUGCA', 'ACAUCCU', 'GUGGGUG', 'AUUCUCG', 'AAGGCUC', 'ACGCGGG', 'UAGCCGC', 'AUCUAUU', 'UGAGUGG', 'GCCCACC', 'AAACCUU', 'GACAUGU', 'GUGGGAG', 'GAGAAAC', 'AGUCAUG', 'AACUAGA', 'GAGGCAG', 'CAAGUCU', 'GGACCCG', 'UCCUUCA', 'GCGAGUC', 'GCACUGA', 'GACGAGG', 'AUGCACC', 'AUCCUUU', 'UCAGACA', 'UUUGCAG', 'GGAGCCA', 'GUCAACA', 'UUCCGCC', 'AACCCUA', 'AUGGCGC', 'CAUGGAU', 'AACCUGG', 'CUCUACA', 'CGUCCUG', 'GGGGUCC', 'CCGAGCC', 'GUGAGUU', 'ACGUAGA', 'GAGAGAU', 'ACUCACA', 'AAAGCCA', 'AUUCAGA', 'AAAGACA', 'CCACCGG', 'GUUCAGA', 'GACUGGC', 'GAGGCUU', 'CUUAUGG', 'ACUCAGA', 'ACUCCAU', 'AGAGAGA', 'ACUGUGU', 'GCUGAGA', 'AACUGCA', 'ACUGAGA', 'CUGAAUG', 'GGCCCCA', 'UGCCCUG', 'GGAUUUC', 'CAAGUGC', 'GGAAACA', 'UGUGCGG', 'ACUGCCC', 'GGGGUUU', 'CUUGAGG', 'CUCCCUU', 'UCAUUAA', 'UCAACAG', 'GAUGUGU', 'UACAUAC', 'ACGUGUG', 'AUACACA', 'UACAGAC', 'GUGUGCA', 'ACGCACG', 'AAGUGCG', 'UGUACAU', 'GAAGGUC', 'GACCCUG', 'GGGAGAG', 'ACUCCAG', 'GAUUGAC', 'UGUAUAA', 'GAGUGUG', 'ACGCUCA', 'AAUUUUA', 'GGUAAGC', 'AUGUCUG', 'GGGACGG', 'GUGUGUG', 'AUCUAGU', 'CAUGACA', 'AGGUUAC', 'GAACUAU', 'CAGGAAC', 'AGACUGA', 'AUACCUC', 'CUGAAAA', 'UAGAGGA', 'GUCCUCU', 'UGCCUUU', 'AUGUGUA', 'CUUGUUA', 'UUCUGCA', 'CGGUGAU', 'ACGUCAU', 'UGUUGAA', 'UCACUGG', 'UACAGUU', 'AACCUGU', 'UGUGUGU', 'AUAAAUA', 'GUGCAUG', 'AUGCAUA', 'GCAUUGU', 'CUUGUGU', 'UACACAC', 'ACAUACU', 'GGGUACA', 'GCAUAUA', 'UGCAUGG', 'GGUGUGA', 'CCGAGGC', 'CGUGUGU', 'UAUACAC', 'CAGCUGA', 'AGUCUUA', 'AAUGAGU', 'GUCUUGU', 'GAAUAUA', 'GAGUUCG', 'AGGACAC', 'AUGUGUU', 'AGCUAGC', 'CUGAGUC', 'GCGGCCG', 'GAGGUGC', 'UUUGGAU', 'CAACUCG', 'ACCAGUC', 'CCCCGAG', 'CUUCUUC', 'GCGAUGG', 'GCCGCCC', 'UUUGGGG', 'GCAGCAG', 'UCUCUGA', 'UUCUGCU', 'GGUCAAG', 'UCUCCCC', 'CAAGGGA', 'GCGCGGG', 'CGCUCGU', 'ACCUCUG', 'UCUAGGA', 'AGCUCAU', 'CCUCCAU', 'UGCAAGG', 'UUGCGGA', 'GUCAUGG', 'CAGGACC', 'CUGGUCC', 'GGAUGAC', 'AACUGAC', 'GGUGGGA', 'GGGAGAU', 'AUCUUAG', 'CAGAUGU', 'AGGUGCC', 'CUUCGCG', 'ACCACCU', 'CGUUGGC', 'GCCGGGC', 'GGACGAG', 'UAUGAGU', 'UUAGGCA', 'UAUACCA', 'CUGCAUC', 'AGUUGUG', 'UAGUGGA', 'CUCCACC', 'GGGAAAG', 'CUGCAGA', 'GUCCCAG', 'AGCAUUG', 'UUUGUGU', 'GUCCAGG', 'AGGAAAG', 'CAGUGCU', 'GGGACGA', 'GCUGGAG', 'CGACUUC', 'AGGGAGC', 'UUCUGAC', 'GAGGAUC', 'GCAGCUG', 'CCACCUC', 'AGAUGGA', 'CCGGGCA', 'GUGUCAC', 'UAAAGGC', 'AGGUAGA', 'GACCUAC', 'AUCUAAC', 'UGGGAGG', 'CUCACCC', 'UGGGUGU', 'CCGGCGG', 'GCUGCGC', 'UGAACCC', 'CUUGGAG', 'UGUGUCA', 'GCACCAU', 'UGGCUGG', 'AUUCAUU', 'GGGCCUG', 'UUGGAGC', 'CCCACAG', 'CUCUAGC', 'GCUUUGC', 'GGAGUGC', 'CAUAGCA', 'AGUGGGC', 'UACCUUU', 'GAGAAUG', 'GGGGACA', 'GGAGAUC', 'GAGGAAU', 'CUGGCUG', 'GCCACAC', 'CAACAAA', 'CAGCACC', 'UGGUUGC', 'ACUUAAU', 'GUUCUCA', 'CAAGCGG', 'UGGAGUU', 'GUGAAUU', 'UGGCAGU', 'UGGACAC', 'GCCCCUG', 'GGUGCUA', 'CUCAUUU', 'UCAUCAA', 'GGGACCC', 'GCCCCCU', 'UGGUCAC', 'UGAUGUC', 'UUCCUGC', 'AUAUCAG', 'GUCUGGG', 'UCCUUGG', 'ACAGGGG', 'GCACUCU', 'UGCGGAC', 'UGACUCC', 'AAAGCCU', 'UGCUGGG', 'UUGAUCU', 'AGGCUCA', 'GAAGCGC', 'CCUCGGG', 'ACUGGAG', 'UGCGCUC', 'ACAGAGU', 'ACAUGGC', 'GGCUGGG', 'CCGAAAC', 'UUGAAGG', 'UCUGCCA', 'GGUGCCA', 'AGAUUGU', 'CCAAUGA', 'ACAUACA', 'AUAACAU', 'UGUGUGC', 'AUACAUG', 'AGGGCAG', 'AACUCAC', 'CUAUGCA', 'UCAUGAG', 'GAGUUCA', 'GCAUCUG', 'UCUGGGU', 'CCCAGGU', 'AUGGGUC', 'GGGGAAA', 'AAUGGGG', 'GCACCCC', 'GUGGACA', 'GUUGGGG', 'CUUUAAA', 'AGCUUUC', 'GGACACU', 'ACAGGUG', 'UCAGACC', 'CCUAACA', 'UCUGCUG', 'UGGGAAC', 'UGUGACA', 'GUACCAU', 'AGCUUUG', 'AGGCUUC', 'AUUUAGC', 'AGGCUCU', 'GUCAGUC', 'UCUAGCC', 'GUGGUGC', 'UCUACUC', 'AGCACCC', 'GAGGGAG', 'AACCUCU', 'AGGGGGC', 'CGCUCUG', 'GAGCAAG', 'CUGCUUA', 'GUUUCAG', 'UCUCUAA', 'GUCUAGG', 'AUGGAUG', 'AGGGCCA', 'UCUGCCU', 'CACUCCA', 'UAAUACA', 'CUUAGCA', 'CAUAGAA', 'AGUUUGU', 'GCUCAUU', 'ACAGAAC', 'UCCUGGC', 'CCUGGCC', 'CACUUUG', 'CACUCUG', 'GAGAUCC', 'AAUAGCC', 'CCUGGGA', 'GCGGCGG', 'GCCCUCA', 'GGAGCUC', 'GGUAGUA', 'GUAUCCC', 'UAAGGUA', 'UGGAGGU', 'CGAUCGU', 'GCUUAUC', 'UGGACUU', 'CAGUACU', 'GCUGCCA', 'GCUUGUC', 'GAAUCCC', 'CCUAAAG', 'AGGUUGU', 'UCCCCAC', 'GGCUGGA', 'GCUCCCG', 'UUACAUG', 'UAGAUCG', 'CGAAUCC', 'UUGUUUG', 'CAUCCGG', 'UGUGCUA', 'GGUCUGU', 'CUCUCCA', 'GGGCAGA', 'UAGAGCA', 'GAGGAGG', 'AGCUCAG', 'CAGAGGA', 'CUGGAGA', 'UUGAUAG', 'AGGUUAG', 'AUCUCAU', 'UUGGGGC', 'GCUUGUG', 'CGCGGGA', 'UGCACAC', 'GUAGAUC', 'GUCCAGU', 'CUGCCUG', 'GGGGAAC', 'ACACUUG', 'CGGGCGG', 'AAUUGGG', 'UGUGGGG', 'UGGAGCG', 'GGCGCCC', 'CGUGGGG', 'UGAUGUU', 'CUGGAGC', 'UGGCAGA', 'CGGGUGG', 'GGUCUAG', 'UAUGCGA', 'GAAAGGC', 'CUCCUGC', 'GUUUACC', 'UUGGUUG', 'GUCUCAG', 'CUGGGCA', 'GCCCUGU', 'ACUUGUG', 'UAAGUGA', 'AGGCGGC', 'ACCCUUU', 'GUCACUA', 'GAAAAGC', 'CGAGGCA', 'CAGUCAU', 'GGAGGUC', 'GGGCCCU', 'UCACCAC', 'UUAGCUG', 'GGCAGCU', 'UUCUGGC', 'UGGCGCC', 'AUUGUAA', 'UAAGGCA', 'CCGGAAG', 'UGAUUCA', 'CCCCAUC', 'AGCAGUU', 'GUAUAAC', 'GAGGGUG', 'CAGGGCU', 'AGUCAGG', 'UACGCAC', 'AAGUCUU', 'CUUGGGA', 'CCGUCGC', 'UAAGGAU', 'UAUUAGU', 'AGCGUUG', 'ACGGCGG', 'GGGGCGU', 'UGCAUAC', 'CUCGGAU', 'GAGGGCU', 'UCAGCAG', 'GUAAAGU', 'UUUUCCU', 'AGCACUG', 'AGCAGGU', 'AGUGCGG', 'CAGCCUU', 'GCCCUGG', 'ACAAUGC', 'CCCCAGA', 'AGUGUGG', 'GAUGUUG', 'AGUGUUG', 'UGGGCAU', 'GUAGACU', 'CGAUCCC', 'UAGCCUU', 'GCUAAGG', 'UGGGAAG', 'UCAAGUG', 'GAACUGG', 'AUGGGAA', 'AGUUAUA', 'AAUUGGC', 'AGUUUAU', 'CAUUUGU', 'GGUUCAC', 'CCCUGGA', 'GUAAGUG', 'AGGAUGC', 'CUUUCCU', 'CAAGGAU', 'GGAAUGG', 'GCAGAUG', 'AAUUGAU', 'GAGGCAC', 'CUGGCGC', 'ACGCAGG', 'UGCCCAC', 'GGGCCUU', 'GGAGACU', 'UCAAUGA', 'UGCAAUG', 'UCUUGCU', 'UACACUC', 'CAAAGGC', 'GCGACCU', 'GGCUGUC', 'ACAUGGG', 'GCGGGUG', 'AGCAUGU', 'CGAAACC', 'CUGGCCU', 'UAAGUUC', 'UCCGUAU', 'CAAAGAG', 'AGCAAAC', 'CAGGGGA', 'CUGCAAC', 'CUAAUCC', 'CAAGUUC', 'GGGAGCA', 'UAAACAU', 'ACUGAAA', 'UUGCAGU', 'CAUGUAU', 'UGAGUUU', 'UAUCAGG', 'GCGGGCU', 'CUGAAGC', 'ACAGGCA', 'UUGUAGG', 'UUCUCUC', 'GUAAGGG', 'AGUCUCU', 'GCAGAAC', 'UGUCCUU', 'GCCAGGA', 'GGUGAUG', 'GCAGGAA', 'ACCUUCU', 'UCCUCAC', 'CAUGUGA', 'UGGGAAA', 'GUGUGUC', 'CCUGGGG', 'CUGAUCU', 'AAUAACC', 'UGGAUGG', 'AAUUCCG', 'UCUGGGA', 'ACGUGCG', 'UGCCUGU', 'CACUCUC', 'AGGGGGG', 'GCCUUCC', 'GGGGGUA', 'CAUGGUU', 'GGGCUGA', 'CAAUCAC', 'ACAGGGA', 'UACCUGA', 'GGGAGAC', 'CUCUACU', 'CUGUCCA', 'AACUUCC', 'AGGAGGA', 'UUGGUCU', 'GUGGAAG', 'UCACUUC', 'GCUGAGG', 'GAGCCUG', 'AGGCCAC', 'CUCCUCC', 'CACAAUG', 'GACUACU', 'GUUCUGC', 'UGAGGGC', 'CACUCCC', 'GAGGAUG', 'CCUUAAU', 'GAGGAAG', 'CCCCUUU', 'CACCAUC', 'UGAGGGG', 'CUGAGCU', 'GGGGAUU', 'AGCUGUU', 'GGGAGGA', 'AGGUGCU', 'AGGGGGU', 'UGGCCUC', 'GGGGUGG', 'CCAACCU', 'UGUAGGU', 'GGCGGCA', 'GGUGUUU', 'UGGUGGA', 'CCUCUGC', 'GUUCUUU', 'UGAGCCA', 'GGCGCGC', 'AGCUGUA', 'CACAGGU', 'UGGAGCU', 'CAUCUGG', 'AGGCAGC', 'CUCACUU', 'UACCUUC', 'UGGGGGA', 'UAACUCG', 'UAUGGUU', 'GGGGUGA', 'UGCUUCC', 'UGAGAGC', 'AACUCUU', 'AAAGGGA', 'UUCUUCU', 'AUCCAGU', 'CAGCCUC', 'GUCCUUC', 'AAGUGGG', 'UCCCCUU', 'CCUGGAG', 'UCUGUCU', 'UGUAUUU', 'UUUUUUC', 'UGGUGGG', 'AGGGGCA', 'UGACCCA', 'GCAGCCA', 'GGGCCUC', 'CACCUGU', 'GAGUGGG', 'GACCGGC', 'AGAAGAG', 'GGAUCUC', 'CAGCUGC', 'GGGAAAC', 'UGCACCG', 'AGGAUGA', 'GCCAACC', 'GACUGGG', 'CCUCUUU', 'AGGACCU', 'AUCUGCC', 'GUGGGAU', 'GCCUCUU', 'UGGGGUG', 'UUCUUGU', 'CUGUGCU', 'GAGCACC', 'GCUGUGU', 'CAUCUUU', 'GGGGGCG', 'CAAGCGC', 'AGCAGAG', 'GUUGUUC', 'CACGCCA', 'GGGGGAG', 'UACUGUA', 'ACGGUGG', 'ACUCUAA', 'UGAGGCA', 'CUCCACU', 'CUGGGGA', 'CUCUCCU', 'AGGGAAG', 'UCUCCCU', 'UGAGGCG', 'AGGCGUU', 'AGGGGUG', 'CGGCUUC', 'GGGAGGC', 'UGUGUCU', 'UGGGGGG', 'GACUUGC', 'UGAGGAG', 'AACAGCC', 'UGGAGGA', 'GGCCCCU', 'UGGAAGG', 'GAAUCAA', 'CUGAAAG', 'ACUUUCU', 'ACUGAGG', 'GGGGAGC', 'UUUUCAC', 'GGGACAG', 'ACCCAUU', 'GGGUGGA', 'GACCUCU', 'UGGGGAC', 'GUCCACU', 'CCAGGGU', 'GCCCUGC', 'GGAGAAC', 'CAGCCCU', 'GCCUGUG', 'CUCACUG', 'UUGCGUU', 'CAAAGGU', 'ACGAUCU', 'UGGGAGU', 'GUGUCCC', 'GCACAGG', 'GGUGUCU', 'AACAGGC', 'CAAACCU', 'UGGGCAG', 'GAGCUGC', 'AGGAAGG', 'UUCAGCU', 'ACCCACC', 'GCUCACA', 'UGGCUUC', 'GCUCUCU', 'CCCGGGU', 'UCCGUGG', 'CCACUCC', 'UGUGCCU', 'GGGGGUG', 'CAGAAGA', 'CCAUCCA', 'GUGCUUC', 'UGGGGUC', 'CUUUUCC', 'GGAAAGG', 'GGUUCUC', 'GGAGGAU', 'CUGCUUC', 'AGGAGAG', 'GACCUGU', 'AUGAAGA', 'CACACUU', 'UGGGUGC', 'CUGUGCA', 'CCUGCCU', 'GAGGGUU', 'CCCUGCU', 'UGAGCAG', 'CACCCUG', 'CACCUUG', 'GGGAUGG', 'ACCCCUC', 'UAGGAGU', 'CUCUGUG', 'CUGGCAA', 'CAAGCCU', 'CAGCAGU', 'GUGAGCU', 'AGCACAA', 'UCUGAGA', 'GUGCUUU', 'CGGGUUU', 'CUUCUCU', 'CUGGACG', 'GUGGCUC', 'GGGGAGG', 'UUGCUCU', 'CUGAGAG', 'ACCCUCU', 'UCCUCUC', 'CUCCAGG', 'GCCUCCA', 'CCGGGAG', 'UGGCCUU', 'GCUGGGA', 'CUGAGCC', 'GCGGGGU', 'CGUCCUC', 'AGGUGGC', 'GUCUCCU', 'GUAGAAG', 'ACUGCUC', 'GAGGGAU', 'AUACGGA', 'GCCACCC', 'GAGGGAA', 'AGAGACA', 'GUCCCUU', 'GUGAAAG', 'CUGUGCC', 'UGUGGUG', 'UGCAGCC', 'CUCCCCC', 'GUAGGGU', 'GAACCAC', 'CACCCUC', 'CUUCCAU', 'CUGAAGG', 'GCAGCCC', 'CAGGAGA', 'CUCAGCC', 'CACCAGG', 'GAGCCAC', 'GUGGGAA', 'CUCUGCC', 'UCCUGUG', 'CCAACUU', 'UGCUACU', 'GUGGAGG', 'CUCACCU', 'AAUUUGA', 'AGAGGGG', 'UCGUUCC', 'CCGGGGA', 'GUUUGCU', 'UGAGUGC', 'CUACUCU', 'AUGGAAA', 'CUAACUU', 'UGAGAGG', 'CAGAGGC', 'AGGGCCC', 'GUGUAGC', 'UGGUUUU', 'GGGUUGG', 'CUACCCA', 'GCAGGGG', 'GAUUGCU', 'CUGGUGG', 'UGAGGCU', 'UGGGGGC', 'UGACGUC', 'AGGGGAG', 'UAGCUUC', 'GGGGGAA', 'GGGAAGA', 'AGGCCUU', 'AACCAUG', 'GACUACC', 'GCCAAGG', 'ACUUUUU', 'GACUUCU', 'UAGGAGC', 'UCCCCUG', 'ACGGGCA', 'CCCCCCU', 'CGGGGCU', 'CUGCUCU', 'AGAGGAU', 'CCUCUGA', 'AGUGGGG', 'AGCUGGC', 'AGAGGAG', 'CCGUUCU', 'GGCAGGU', 'GCUGGCU', 'CAAGUGG', 'UGACCUU', 'CGUGGGA', 'CAGUGUC', 'CUCAGAG', 'UUGAUGU', 'UAGGGGU', 'GUGGCUU', 'AAGGCAA', 'UUCCAUC', 'CUGAAGA', 'ACUCUGC', 'CGGGGGU', 'CUUGGUC', 'GAAGACA', 'CUGGGGG', 'GGGGAAG', 'AACUCUA', 'AAAAACA', 'AACAUUG', 'UAUAUGU', 'CUUUCCC', 'GGAGUGA', 'AACACAC', 'CUGGGUC', 'ACCUCAA', 'GUUUUCU', 'UCAGUCC', 'CUCUUUA', 'GCUCUGA', 'GGAGAGC', 'GGAAUUU', 'ACUUGUA', 'GAGAACC', 'GCAGGGU', 'UAACCUG', 'GUGUUAG', 'AUCCGGG', 'GUCUCUU', 'GCACUUG', 'CCAGGAC', 'UCUCGAA', 'GGUAGGC', 'CCACUGA', 'CGUAGAC', 'UUCUAGA', 'CCAGGGA', 'GACACAG', 'CUGAUGA', 'GGCGACC', 'GUGUUCG', 'ACACAGA', 'GACUGCU', 'UGAGUUG', 'UGGGGAA', 'ACAUGGA', 'GGUUGAA', 'GUUAGGG', 'AUUGUCU', 'UGUUUUC', 'AACGUCU', 'CUGACUU', 'CUAAAAC', 'UGGACCG', 'UGGUAGC', 'AUCCUGU', 'UGCGGAU', 'CCAUCCU', 'UGGAGAG', 'GCCGUUG', 'CUUGCAU', 'AGUAGAU', 'UUGGGCU', 'AUGGCUC', 'GGCGCCU', 'AGCGAAG', 'CAAGAGA', 'CACUGGU', 'ACUACCA', 'AGUCCAG', 'CGCGUUC', 'GACAGCA', 'AUCCUCU', 'UUUUGAU', 'AGUGGAA', 'CAGUUUG', 'GAUGUCC', 'GAAUUGC', 'CGCGGGG', 'UCCUGGU', 'AAGGGCA', 'CUGGCUU', 'AAGGGGC', 'CCCUACU', 'GAGUCGC', 'GAGCGGG', 'GGCCCAC', 'UUAUGAC', 'CAGGCUG', 'AUAGUUG', 'UCCUUGC', 'GUGGGCG', 'CCACCUU', 'GGAGCCG', 'GGAGUCA', 'UUGCCCU', 'UCUGUAA', 'AGAAAGA', 'CUGCUUG', 'GUGUGGA', 'AGUUUAC', 'CUCUUAG', 'GCCAGGU', 'UAGCCCG', 'AUGCAGC', 'AGCCAUC', 'GGGACUG', 'ACUCUGA', 'GGUGCAG', 'GUACCAC', 'AGCUGGG', 'UUCAGAU', 'GGCUCCG', 'AGCGGGC', 'GGUGGAC', 'CGGCUGC', 'UUGACUG', 'UCCUUAC', 'CGGAGGA', 'CAGUCUG', 'GAGGUGU', 'CCCAUUC', 'UUACGAA', 'GUUUCAU', 'GGGCAGU', 'CCGGUGC', 'CUUUUUG', 'AAGCCCA', 'AGAGCCC', 'CAGGCUA', 'GACUCAC', 'GUGGUGA', 'CGGCUGU', 'UACAAGG', 'UUCCUUA', 'GGCUGGC', 'AGGAUCC', 'UUCCAGU', 'CAGGCAC', 'CAAACCA', 'UGCUCAG', 'GAUCAGG', 'GAAAGGG', 'GGCACGA', 'CUGUGGG', 'UCCGUGU', 'AAGAAUC', 'UUGCCUG', 'GCUGACU', 'GACGUCC', 'GAGAAUU', 'UUAUUAC', 'GCGUACC', 'CCUUCCG', 'GACUUGG', 'CUUCCAC', 'UGCUCGG', 'CUUAAAA', 'AGUGCCC', 'GGCGGGG', 'CACUCAG', 'GCUGAGC', 'UAGGGAA', 'UAGAGCC', 'GCUGUCA', 'CAGAACA', 'AAAUAGG', 'AAGCGCA', 'CCAUCUG', 'UGGGGAG', 'ACUGAAG', 'AAGGAUU', 'UGGCUCU', 'GGCACGG', 'ACAAGGA', 'CUAUGGU', 'CAAGGCG', 'CCAUAAA', 'GGAGGAA', 'GCGAGCG', 'CACGCGG', 'CUCCUGU', 'GGCGACU', 'GGAUAUG', 'GACUCUG', 'GGAAACU', 'AAAGGCC', 'CGCCGCG', 'AGCGUGG', 'CUCCGCC', 'CACCCAU', 'CGCAAUA', 'CUCGUGU', 'ACAAACA', 'GUGCCUU', 'AAAAGUG', 'AGAAGGG', 'GGCUGCA', 'GUGAUUU', 'UACAUUA', 'CCAGUUA', 'UGCUGAC', 'GGCACUC', 'CUGCCUU', 'GGCACAG', 'AGGGUCU', 'GUGCCUA', 'GAGUACA', 'ACCAGGC', 'GGGCCAU', 'AUGGGUG', 'CUCUGUC', 'GGAGCUG', 'UGAGGGU', 'CCAGGCC', 'CAGCGCC', 'AAUUAGG', 'GUGGUCC', 'UGAACUG', 'UUGCAUG', 'CAGACGG', 'UCGGUUA', 'UACAGCU', 'UUGGGGG', 'GGGACUC', 'AAGGGGA', 'UCAUUCA', 'AUGUCAC', 'GGACAGC', 'AGGCUGG', 'CUCUGAA', 'UGGCCUG', 'CCGUGGG', 'CCAUGGA', 'UUAGUGU', 'AACACCA', 'GUGGGGU', 'GGGUGAC', 'CCGGUCC', 'CAAGGUC', 'GGGUAUA', 'AGACAAG', 'AGGGGCU', 'CCUUUCU', 'ACUCAGG', 'GUGAAUC', 'GGAAGGA', 'GGUGCCU', 'GCCGGAG', 'UGGAAUC', 'UGAGUCU', 'AGCCGGG', 'GUGUUUC', 'AACCAAC', 'ACUAUGC']
done
## the mature_uniquseeds_logoplt for "mmu" is the following: 
![image](https://user-images.githubusercontent.com/130226484/234354535-a2454698-34a0-4486-bd4b-d64b6bf10473.png)
## the calculation of GC frequency for "mmu":
 A         C         G         T
0  0.241814  0.243703  0.291562  0.222922
1  0.236146  0.211587  0.301008  0.251259
2  0.205919  0.250630  0.335013  0.208438
3  0.217884  0.220403  0.307305  0.254408
4  0.211587  0.241814  0.286524  0.260076
5  0.217254  0.246851  0.301008  0.234887
6  0.212217  0.270151  0.253778  0.263854
