
#this is where it begins
from multi_key_dict import multi_key_dict
import module_manager

class FASTATranslator(object):
    def __init__(self, path):
        self.path = path
        self.aminoacidsequence = []
        self.title = None
        self.codonDict = self.makeCodonDict()
        sequence = self.FASTAtranslate(path)
    
    def getSequence(self):
        return self.aminoacidsequence

    def getTitle(self):
        return self.title

    def FASTAtranslate(self,path):
        dnaSequence = self.readToString(path)
        rnaSequence = self.transcribe(dnaSequence)

        self.aminoacidsequence = self.translate(rnaSequence)
    
    #from http://www.cs.cmu.edu/~112/notes/notes-strings.html#basicFileIO
    def readFile(self, path):
        with open(path, "rt") as f:
            return f.read()

    def readToString(self,path):
        raw = self.readFile(path)
        dnaSequence = ""
        for line in raw.splitlines():
            if (line[0]==">"):
                self.title = line[0][1:]
            else:
                dnaSequence += (line)
        return dnaSequence
        
    def transcribe(self,dnaSequence):
        rnaSequence = ""
        rnaSequence = dnaSequence.replace("T", "U")
        return rnaSequence

    def translate(self,rnaSequence):
        aminoacidsequence = []
        for i in range(0,len(rnaSequence),3):
            codon = rnaSequence[i:i+3]
            if (codon == "stop" or len(codon)!=3):
                break
            else:
                aminoacidsequence.append(self.codonDict[codon])
        return aminoacidsequence


    def makeCodonDict(self):
        codonDict = multi_key_dict() #TODO: replace with amino acid objects
        codonDict["UUU","UUC",] = "Phe"
        codonDict["UUA","UUG","CUU","CUC","CUA","CUG"] = "Leu"
        codonDict["AUU","AUC","AUA"] = "Ile"
        codonDict["AUG"] = "Met"
        codonDict["GUU","GUC","GUA","GUG"] = "Val"
        codonDict["UCU","UCC","UCA","UCG","AGU","AGC"] = "Ser"
        codonDict["CCU","CCC","CCA","CCG"] = "Pro"
        codonDict["ACU","ACC","ACA","ACG"] = "Thr"
        codonDict["GCU","GCC","GCA","GCG"] = "Ala"
        codonDict["UAU","UAC"] = "Tyr"
        codonDict["UAA","UAG","UGA"] = "stop"
        codonDict["CAU","CAC"] = "His"
        codonDict["CAA","CAG"] = "Gln"
        codonDict["AAU","AAC"] = "Asn"
        codonDict["AAA","AAG"] = "Lys"
        codonDict["GAU","GAC"] = "Asp"
        codonDict["GAA","GAG"] = "Glu"
        codonDict["UGU","UGC"] = "Cys"
        codonDict["UGG"] = "Trp"
        codonDict["CGU","CGC","CGA","CGG","AGA","AGG"] = "Arg"
        codonDict["GGU","GGC","GGA","GGG"] = "Gly"
        return codonDict

    
testTranslator = FASTATranslator("test.fasta")
final = testTranslator.getSequence()
module_manager.review()
#print(final) 


class AminoAcid(object):
    def __init__(self):
        self.alphaCarbon = AlphaCarbon()
        self.amine = Amine()
        self.carboxyl = Carboxyl()

class FunctionalGroup(object):
    def __init__(self, charge,sulfide,hbond):
        self.charge = charge
        self.sulfide = sulfide
        self.hbond = hbond

class AlphaCarbon(FunctionalGroup):
    def __init__(self):
        super().__init__(self, 0,False,False)

class Amine(FunctionalGroup):
    def __init__(self):
        super().__init__(self,1,False,False)

class Carboxyl(FunctionalGroup):
    def __init__(self):
        super().__init__(self,-1,False,False)

class Phenyalanine(AminoAcid):
    def __init__(self):
        super.__init__(self)
        self.sidechain = "idk"