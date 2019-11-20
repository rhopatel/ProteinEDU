
from multi_key_dict import multi_key_dict
import module_manager
import copy

#biological inheritance hierarchy:

#AminoAcid superclass (will contain applicable method)
class AminoAcid(object):
    def __init__(self, sidechain):
        self.alphaCarbon = AlphaCarbon()
        self.amine = Amine()
        self.carboxyl = Carboxyl()
        self.sidechain = sidechain
  

#FunctionalGroup superclass (will contain applicable method)
class FunctionalGroup(object):
    def __init__(self, charge,sulfide,hbond):
        self.charge = charge
        self.sulfide = sulfide
        self.hbond = hbond


#three common functional groups
class AlphaCarbon(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,False)

class Amine(FunctionalGroup):
    def __init__(self):
        super().__init__(1,False,False)

class Carboxyl(FunctionalGroup):
    def __init__(self):
        super().__init__(-1,False,False)


#all the unique amino acids
class Alanine(AminoAcid):
    def __init__(self):
        sidechain = Methyl()
        super().__init__(sidechain)
    

class Arginine(AminoAcid): #TODO: add other sidechains
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Asparagine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)

class AsparticAcid(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)

class Cysteine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Glutamine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class GlutamicAcid(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Glycine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Histidine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Isoleucine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)

class Leucine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Lysine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Methionine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)

class Phenylalanine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Proline(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Serine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Threonine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Tryptophan(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Tyrosine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


class Valine(AminoAcid):
    def __init__(self):    
        sidechain = None 
        super().__init__(sidechain)


#unique sidechain functional groups (add the rest)
class Methyl(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,False)



#fasta file to AminoAcid sequence translator
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
                aminoacidsequence.append(copy.copy(self.codonDict[codon]))
        return aminoacidsequence


    def makeCodonDict(self):

        codonDict = multi_key_dict() 
        codonDict["UUU","UUC"] = Phenylalanine()
        codonDict["UUA","UUG","CUU","CUC","CUA","CUG"] = Leucine()
        codonDict["AUU","AUC","AUA"] = Isoleucine()
        codonDict["AUG"] = Methionine()
        codonDict["GUU","GUC","GUA","GUG"] = Valine()
        codonDict["UCU","UCC","UCA","UCG","AGU","AGC"] = Serine()
        codonDict["CCU","CCC","CCA","CCG"] = Proline()
        codonDict["ACU","ACC","ACA","ACG"] = Threonine()
        codonDict["GCU","GCC","GCA","GCG"] = Alanine()
        codonDict["UAU","UAC"] = Tyrosine()
        codonDict["UAA","UAG","UGA"] = "stop"
        codonDict["CAU","CAC"] = Histidine()
        codonDict["CAA","CAG"] = Glutamine()
        codonDict["AAU","AAC"] = Asparagine()
        codonDict["AAA","AAG"] = Lysine()
        codonDict["GAU","GAC"] = AsparticAcid()
        codonDict["GAA","GAG"] = GlutamicAcid()
        codonDict["UGU","UGC"] = Cysteine()
        codonDict["UGG"] = Tryptophan()
        codonDict["CGU","CGC","CGA","CGG","AGA","AGG"] = Arginine()
        codonDict["GGU","GGC","GGA","GGG"] =  Glycine()
        return codonDict

    
testTranslator = FASTATranslator("test.fasta")
final = testTranslator.getSequence()
module_manager.review()
print(final) 


