#Project: ProteinEDU
#Author: Rohan Patel


#imports
from multi_key_dict import multi_key_dict
import module_manager
import copy
import pygame
import sys
import os

#biological inheritance hierarchy:

#AminoAcid superclass (will contain applicable method)
class AminoAcid(object):
    def __init__(self, sidechain):
        self.alphaCarbon = AlphaCarbon()
        self.amine = Amine()
        self.carboxyl = Carboxyl()
        self.sidechain = sidechain

    def __repr__(self):
        return self.__class__.__name__
    
    def addToChain(self):
        self.amine = None
        self.carboxyl = None

#FunctionalGroup superclass (will contain applicable method)
class FunctionalGroup(object): #TODO: Update to be more biologically significant
    def __init__(self,charge,sulfide,hbond,polar):
        self.charge = charge
        self.sulfide = sulfide
        self.hbond = hbond
        self.polar = polar

#three common functional groups
class AlphaCarbon(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,False,False)
class Amine(FunctionalGroup):
    def __init__(self):
        super().__init__(1,False,False,True)
class Carboxyl(FunctionalGroup):
    def __init__(self):
        super().__init__(-1,False,False,True)

#all the unique amino acids
class Alanine(AminoAcid):
    def __init__(self):
        sidechain = AlanineSideChain()
        super().__init__(sidechain)   
class Arginine(AminoAcid): #TODO: add other sidechains
    def __init__(self):    
        sidechain = ArginineSideChain() 
        super().__init__(sidechain)
class Asparagine(AminoAcid):
    def __init__(self):    
        sidechain = AsparagineSideChain() 
        super().__init__(sidechain)
class AsparticAcid(AminoAcid):
    def __init__(self):    
        sidechain = AsparticAcidSideChain() 
        super().__init__(sidechain)
class Cysteine(AminoAcid):
    def __init__(self):    
        sidechain = CysteineSideChain() 
        super().__init__(sidechain)
class Glutamine(AminoAcid):
    def __init__(self):    
        sidechain = GlutamineSideChain() 
        super().__init__(sidechain)
class GlutamicAcid(AminoAcid):
    def __init__(self):    
        sidechain = GlutamicAcidSideChain() 
        super().__init__(sidechain)
class Glycine(AminoAcid):
    def __init__(self):    
        sidechain = GlycineSideChain() 
        super().__init__(sidechain)
class Histidine(AminoAcid):
    def __init__(self):    
        sidechain = HistidineSideChain() 
        super().__init__(sidechain)
class Isoleucine(AminoAcid):
    def __init__(self):    
        sidechain = IsoleucineSideChain() 
        super().__init__(sidechain)
class Leucine(AminoAcid):
    def __init__(self):    
        sidechain = LeucineSideChain() 
        super().__init__(sidechain)
class Lysine(AminoAcid):
    def __init__(self):    
        sidechain = LysineSideChain()
        super().__init__(sidechain)
class Methionine(AminoAcid):
    def __init__(self):    
        sidechain = MethionineSideChain() 
        super().__init__(sidechain)
class Phenylalanine(AminoAcid):
    def __init__(self):    
        sidechain = PhenylalanineSideChain() 
        super().__init__(sidechain)
class Proline(AminoAcid):
    def __init__(self):    
        sidechain = ProlineSideChain() 
        super().__init__(sidechain)
class Serine(AminoAcid):
    def __init__(self):    
        sidechain = SerineSideChain() 
        super().__init__(sidechain)
class Threonine(AminoAcid):
    def __init__(self):    
        sidechain = ThreonineSideChain() 
        super().__init__(sidechain)
class Tryptophan(AminoAcid):
    def __init__(self):    
        sidechain = TryptophanSideChain() 
        super().__init__(sidechain)
class Tyrosine(AminoAcid):
    def __init__(self):    
        sidechain = TyrosineSideChain() 
        super().__init__(sidechain)
class Valine(AminoAcid):
    def __init__(self):    
        sidechain = ValineSideChain() 
        super().__init__(sidechain)

#unique sidechain functional groups (add the rest)
class AlanineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,False,False)
class ArginineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(1,False,True,True)
class AsparagineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,True,True)
class AsparticAcidSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(-1,False,False,True)
class CysteineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,True,False,False)
class GlutamineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,True,True)
class GlutamineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,True,True)
class GlutamicAcidSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(-1,False,False,True)
class GlycineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
class HistidineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(1,False,True,True)
class IsoleucineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
class LeucineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
class LysineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
class MethionineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,True,False,False)
class PhenylalanineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
class ProlineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
class SerineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
class ThreonineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
class TryptophanSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
class TyrosineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
class ValineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)

#fasta file to AminoAcid sequence translator
class FASTA(object):
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
        self.bond()
        
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

    def bond(self):
        for i in range(1, len(self.aminoacidsequence)-2):
            currentAminoAcid = self.aminoacidsequence[i]
            if (isinstance(currentAminoAcid,AminoAcid)):
                currentAminoAcid.addToChain()

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
    
def FASTAtest():
    testTranslator = FASTA("FADS.fasta")
    final = testTranslator.getSequence()
    module_manager.review()
    print(final) 

def main():
    FASTAtest()
    
    pygame.init()
    screen = pygame.display.set_mode((800,800))
    screen.fill((255,255,255))
    pygame.display.set_caption('first game')
    done = False
    
    screenSize,screenSize2 = pygame.display.get_surface().get_size()
    loadBarHeight = screenSize*(0.10)
    loadBarWidth = screenSize*(0.80)
    gameScreenWidth = loadBarWidth
    
    gameScreen = pygame.Rect(0,loadBarHeight,gameScreenWidth,screenSize)
    loadBar = pygame.Rect(0,0,gameScreenWidth,loadBarHeight)
    infoBar = pygame.Rect(gameScreenWidth,0,screenSize,screenSize)
    green = (0,255,0)
    black = (0,0,0)
    gamefile = FASTA("FADS.fasta")
    while not done:
        pygame.time.delay(100)
        screen.fill((255,255,255))
        pygame.draw.rect(screen,black,gameScreen, 10)
        pygame.draw.rect(screen,black,loadBar,10)
        pygame.draw.rect(screen,black,infoBar,10)
        pygame.draw.rect(screen,green,(loadBarWidth*0.1,loadBarHeight*0.3,
                                        loadBarWidth*0.125,loadBarHeight*0.45))

        for event in pygame.event.get():
            if (event.type == pygame.QUIT):
                done = True
           
        pygame.display.update()
    pygame.quit()

if __name__ == '__main__':
    main()

