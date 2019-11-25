#Project: ProteinEDU
#Author: Rohan Patel


#imports
from multi_key_dict import multi_key_dict
import module_manager
import copy
import pygame
import sys
import os
from tkinter import filedialog
from tkinter import *
import random


#biological inheritance hierarchy:

#AminoAcid superclass (will contain applicable method)
class AminoAcid(object):
    def __init__(self, sidechain,color):
        self.alphaCarbon = AlphaCarbon()
        self.amine = Amine()
        self.carboxyl = Carboxyl()
        self.sidechain = sidechain
        self.x = None
        self.y = None
        self.r = 30
        self.color = pygame.Color(color)
        self.particle = None

    def __repr__(self):
        return self.__class__.__name__
    
    def addToChain(self):
        self.amine = None
        self.carboxyl = None
    
    def draw(self,screen):
        self.particle = Particle(self.color,self.x,self.y,self.r)
        self.particle.draw_particle(screen)
        
class Particle(object):
    def __init__(self,color, x, y ,r):
        self.x = x
        self.y = y
        self.color = color
        self.r = r

    def draw_particle(self,screen):
        pygame.draw.circle(screen,self.color,(int(self.x),int(self.y)),self.r)


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
        super().__init__(sidechain,sidechain.color)   
class Arginine(AminoAcid):
    def __init__(self):    
        sidechain = ArginineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Asparagine(AminoAcid):
    def __init__(self):    
        sidechain = AsparagineSideChain() 
        super().__init__(sidechain,sidechain.color)
class AsparticAcid(AminoAcid):
    def __init__(self):    
        sidechain = AsparticAcidSideChain() 
        super().__init__(sidechain,sidechain.color)
class Cysteine(AminoAcid):
    def __init__(self):    
        sidechain = CysteineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Glutamine(AminoAcid):
    def __init__(self):    
        sidechain = GlutamineSideChain() 
        super().__init__(sidechain,sidechain.color)
class GlutamicAcid(AminoAcid):
    def __init__(self):    
        sidechain = GlutamicAcidSideChain() 
        super().__init__(sidechain,sidechain.color)
class Glycine(AminoAcid):
    def __init__(self):    
        sidechain = GlycineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Histidine(AminoAcid):
    def __init__(self):    
        sidechain = HistidineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Isoleucine(AminoAcid):
    def __init__(self):    
        sidechain = IsoleucineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Leucine(AminoAcid):
    def __init__(self):    
        sidechain = LeucineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Lysine(AminoAcid):
    def __init__(self):    
        sidechain = LysineSideChain()
        super().__init__(sidechain,sidechain.color)
class Methionine(AminoAcid):
    def __init__(self):    
        sidechain = MethionineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Phenylalanine(AminoAcid):
    def __init__(self):    
        sidechain = PhenylalanineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Proline(AminoAcid):
    def __init__(self):    
        sidechain = ProlineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Serine(AminoAcid):
    def __init__(self):    
        sidechain = SerineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Threonine(AminoAcid):
    def __init__(self):    
        sidechain = ThreonineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Tryptophan(AminoAcid):
    def __init__(self):    
        sidechain = TryptophanSideChain() 
        super().__init__(sidechain,sidechain.color)
class Tyrosine(AminoAcid):
    def __init__(self):    
        sidechain = TyrosineSideChain() 
        super().__init__(sidechain,sidechain.color)
class Valine(AminoAcid):
    def __init__(self):    
        sidechain = ValineSideChain() 
        super().__init__(sidechain,sidechain.color)

#unique sidechain functional groups (add the rest)
class AlanineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,False,False)
        self.color = "green"
class ArginineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(1,False,True,True)
        self.color = "red"
class AsparagineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,True,True)
        self.color = "purple"
class AsparticAcidSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(-1,False,False,True)
        self.color = "blue"
class CysteineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,True,False,False)
        self.color = "yellow"
class GlutamineSideChain(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,True,True)
        self.color = "purple"
class GlutamicAcidSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(-1,False,False,True)
        self.color = "blue"
class GlycineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"
class HistidineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(1,False,True,True)
        self.color = "red"
class IsoleucineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"
class LeucineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"
class LysineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
        self.color = "red"
class MethionineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,True,False,False)
        self.color = "yellow"
class PhenylalanineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"
class ProlineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"
class SerineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
        self.color = "purple"
class ThreonineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
        self.color = "purple"
class TryptophanSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"
class TyrosineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,True,True)
        self.color = "purple"
class ValineSideChain(FunctionalGroup):
    def __init__(self):    
        super().__init__(0,False,False,False)
        self.color = "green"

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
                dnaSequence += (line.upper())
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
    
    def __repr__(self):
        return str(self.getSequence())

class Game(object):
    def __init__(self):
        self.gameFASTA = None
        self.gameSequence = None
        self.fileLoaded = False
        pygame.init()
        self.screen = pygame.display.set_mode((800,800))
        self.screen.fill((255,255,255))
        pygame.display.set_caption("ProteinEDU")

        #visual assets
        self.gamefont = pygame.font.Font('freesansbold.ttf', 20) 
        self.green = (0,255,0)
        self.black = (0,0,0)
        self.white = (255,255,255)

        self.done = False
        #initialize heights/widths based on window size
        self.screenSize,self.screenSize2=pygame.display.get_surface().get_size()
        self.loadBarHeight = self.screenSize*(0.10)
        self.loadBarWidth = self.screenSize*(0.80)
        self.gameScreenWidth = self.loadBarWidth
        self.gameScreenHeight = self.screenSize*(0.9)
        self.loadButtonX1 = self.loadBarWidth*0.1
        self.loadButtonY1 = self.loadBarHeight*0.3
        self.loadButtonX2 = self.loadBarWidth*0.125
        self.loadButtonY2 = self.loadBarHeight*0.45
        
        #load UI objects
        self.gameScreen =pygame.Rect(0,self.loadBarHeight,self.gameScreenWidth,
                                            self.screenSize)
        self.loadBar = pygame.Rect(0,0,self.gameScreenWidth,self.loadBarHeight)
        self.infoBar = pygame.Rect(self.gameScreenWidth,0,self.screenSize,
                                    self.screenSize)
        self.loadButton =pygame.Rect(self.loadButtonX1,self.loadButtonY1,
                                        self.loadButtonX2,self.loadButtonY2)
        self.loadText = self.gamefont.render("LOAD",True,self.black)
    
    

    def loadFASTA(self): 
        fileName = filedialog.askopenfilename() #prompts user to specify file
        print(fileName)
        if ("fasta" in fileName.lower()):
            self.gameFASTA = FASTA(fileName)
            sequence = self.gameFASTA.getSequence()
            start = random.randint(0,len(sequence)-6) #takes a fragment
            self.gameSequence = sequence[start:start+5]
            self.fileLoaded=True
        else:
            print("WRONG FILE TYPE")
        

    def FASTAtest(self,fileName):
        testTranslator = FASTA(fileName)
        final = testTranslator.getSequence()
        module_manager.review()
        print(final) 

    def getGameSequence(self):
        return self.gameSequence
    
    def initializeGameScreen(self):
        dx = self.gameScreenWidth/(len(self.gameSequence)+1)
        currentX = dx
        currentY = self.gameScreenHeight/2
        for aminoAcid in self.gameSequence:
            aminoAcid.x = currentX
            aminoAcid.y = currentY
            currentX+=dx
            
    
def main():
    game = Game()
    #game.FASTAtest("assets//FADS.fasta")
    
  
    while not game.done:
        pygame.time.delay(100)
        game.screen.fill((255,255,255))
        pygame.draw.rect(game.screen,game.black,game.gameScreen, 10)
        pygame.draw.rect(game.screen,game.black,game.loadBar,10)
        pygame.draw.rect(game.screen,game.black,game.infoBar,10)
        pygame.draw.rect(game.screen,game.green,game.loadButton)
        game.screen.blit(game.loadText,game.loadButton)
        if (game.fileLoaded):
            prevX = None
            prevY = None
            for aminoAcid in game.getGameSequence():
                
                if (prevX!=None):
                    pygame.draw.line(game.screen,game.black,(prevX,prevY),
                                        (aminoAcid.x,aminoAcid.y),10)
                prevX = aminoAcid.x
                prevY = aminoAcid.y
                aminoAcid.draw(game.screen)

        for event in pygame.event.get():
            if (event.type == pygame.QUIT):
                game.done = True
            elif (event.type == pygame.MOUSEBUTTONDOWN):
                mousePos = pygame.mouse.get_pos()
                mouseX,mouseY = mousePos
                left,right,middle=pygame.mouse.get_pressed()
                if (game.loadButton.collidepoint(mousePos) and left):
                    game.loadFASTA()
                    game.initializeGameScreen()
                elif(game.screen.get_at(mousePos)!=game.white \
                    and game.screen.get_at(mousePos)!=game.black \
                    and mouseX<game.gameScreenWidth \
                    and mouseY<game.gameScreenHeight):
                    print("clicked")
            elif (event.type == pygame.MOUSEBUTTONUP):
                pass
                
            elif (event.type == pygame.MOUSEMOTION):
                pass
        
        
        pygame.display.update()
    pygame.quit()


if __name__ == '__main__':
    main()

