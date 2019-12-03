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
import math


#biological inheritance hierarchy:

#AminoAcid superclass (will contain applicable method)
class AminoAcid(object):
    def __init__(self, sidechain,color,abbr):
        self.alphaCarbon = AlphaCarbon()
        self.amine = Amine()
        self.carboxyl = Carboxyl()
        self.sidechain = sidechain
        self.x = None
        self.y = None
        self.r = 30
        self.color = pygame.Color(color)
        self.particle=None
        self.name = self.__class__.__name__
        self.abbr = abbr
        self.satisfied = False
        self.associates = set()

    def __repr__(self):
        return self.name
    
    def addToChain(self):
        self.amine = None
        self.carboxyl = None
    
    def draw(self,game):
        self.particle.draw_particle(game)
        
class Particle(object):
    def __init__(self,name,abbr,color, x, y ,r):
        self.x = x
        self.y = y
        self.color = color
        self.r = r
        self.clicked = False
        self.name = name
        self.abbr = abbr
   
    def draw_particle(self,game):
        screen = game.screen
        self.boundingBox=pygame.draw.circle(screen,self.color,(int(self.x),int(self.y)),self.r)
        self.textRect = pygame.Rect(self.x,self.y,
                                        self.r*2,self.r*2)
        font = pygame.font.Font('freesansbold.ttf', 30) 
        self.nameText = font.render(self.abbr,True,pygame.Color("black"))
        game.screen.blit(self.nameText,self.textRect)


#FunctionalGroup superclass (will contain applicable method)
class FunctionalGroup(object): 
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
        abbr = "Ala"
        super().__init__(sidechain,sidechain.color,abbr) 
class Arginine(AminoAcid):
    def __init__(self):    
        sidechain = ArginineSideChain() 
        abbr = "Arg"
        super().__init__(sidechain,sidechain.color,abbr)
class Asparagine(AminoAcid):
    def __init__(self):    
        sidechain = AsparagineSideChain() 
        abbr = "Asn"
        super().__init__(sidechain,sidechain.color,abbr)
class AsparticAcid(AminoAcid):
    def __init__(self):    
        sidechain = AsparticAcidSideChain()
        abbr = "Asp"
        super().__init__(sidechain,sidechain.color,abbr)
class Cysteine(AminoAcid):
    def __init__(self):    
        sidechain = CysteineSideChain()
        abbr = "Cys"
        super().__init__(sidechain,sidechain.color,abbr)
class Glutamine(AminoAcid):
    def __init__(self):    
        sidechain = GlutamineSideChain()
        abbr = "Gln" 
        super().__init__(sidechain,sidechain.color,abbr)
class GlutamicAcid(AminoAcid):
    def __init__(self):    
        sidechain = GlutamicAcidSideChain() 
        abbr = "Glu"
        super().__init__(sidechain,sidechain.color,abbr)
class Glycine(AminoAcid):
    def __init__(self):    
        sidechain = GlycineSideChain() 
        abbr = "Gly"
        super().__init__(sidechain,sidechain.color,abbr)
class Histidine(AminoAcid):
    def __init__(self):    
        sidechain = HistidineSideChain()
        abbr = "His" 
        super().__init__(sidechain,sidechain.color,abbr)
class Isoleucine(AminoAcid):
    def __init__(self):    
        sidechain = IsoleucineSideChain()
        abbr = "Ile" 
        super().__init__(sidechain,sidechain.color,abbr)
class Leucine(AminoAcid):
    def __init__(self):    
        sidechain = LeucineSideChain()
        abbr = "Leu" 
        super().__init__(sidechain,sidechain.color,abbr)
class Lysine(AminoAcid):
    def __init__(self):    
        sidechain = LysineSideChain()
        abbr = "Lys"
        super().__init__(sidechain,sidechain.color,abbr)
class Methionine(AminoAcid):
    def __init__(self):    
        sidechain = MethionineSideChain() 
        abbr = 'Met'
        super().__init__(sidechain,sidechain.color,abbr)
class Phenylalanine(AminoAcid):
    def __init__(self):    
        sidechain = PhenylalanineSideChain()
        abbr = "Phe" 
        super().__init__(sidechain,sidechain.color,abbr)
class Proline(AminoAcid):
    def __init__(self):    
        sidechain = ProlineSideChain()
        abbr = "Pro" 
        super().__init__(sidechain,sidechain.color,abbr)
class Serine(AminoAcid):
    def __init__(self):    
        sidechain = SerineSideChain()
        abbr = "Ser" 
        super().__init__(sidechain,sidechain.color,abbr)
class Threonine(AminoAcid):
    def __init__(self):    
        sidechain = ThreonineSideChain()
        abbr = "Thr" 
        super().__init__(sidechain,sidechain.color,abbr)
class Tryptophan(AminoAcid):
    def __init__(self):    
        sidechain = TryptophanSideChain()
        abbr = "Trp" 
        super().__init__(sidechain,sidechain.color,abbr)
class Tyrosine(AminoAcid):
    def __init__(self):    
        sidechain = TyrosineSideChain()
        abbr = "Tyr" 
        super().__init__(sidechain,sidechain.color,abbr)
class Valine(AminoAcid):
    def __init__(self):    
        sidechain = ValineSideChain()
        abbr = "Val" 
        super().__init__(sidechain,sidechain.color,abbr)

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
        end = len(rnaSequence)
        i = 0
        while (i<=end):
            codon = rnaSequence[i:i+3]

            if (codon == "stop" or len(codon)!=3):
                pass
            
            else:
                aminoacidsequence.append(copy.copy(self.codonDict[codon]))
            i+=3
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
    def __init__(self,delay):
        self.gameFASTA = None
        self.gameSequence = None
        self.fileLoaded = False

        self.delay = delay
        self.done = False

        pygame.init()
        self.screen = pygame.display.set_mode((1000,1000))
        self.screen.fill((255,255,255))
        pygame.display.set_caption("ProteinEDU")

        #visual assets
        self.gamefont = pygame.font.Font('freesansbold.ttf', 20) 
        self.green = (0,255,0)
        self.black = (0,0,0)
        self.white = (255,255,255)
        self.red = (255,0,0)
        self.blue = (0,0,255)
        self.yellow = (255,255,0)
        self.purple = (255,0,255)
        #initialize heights/widths based on window size
        self.screenSize,self.screenSize2=pygame.display.get_surface().get_size()
        self.loadBarHeight = self.screenSize*(0.10)
        self.loadBarWidth = self.screenSize*(0.70)
        self.gameScreenWidth = self.loadBarWidth
        self.gameScreenHeight = self.screenSize*(0.9)

        self.buttonWidth= self.loadBarWidth*0.1
        self.buttonHeight = self.loadBarHeight*0.5

        self.loadButtonX1 = self.loadBarWidth*0.1
        self.loadButtonY1 = self.loadBarHeight*0.3
        

        self.clearButtonX1 = self.loadBarWidth*0.3
        self.clearButtonY1 = self.loadButtonY1
   
        self.helpButtonX1 = self.loadBarWidth*0.5
        self.helpButtonY1 = self.loadButtonY1

        self.infoImageX= self.gameScreenWidth+10
        self.infoImageY = self.gameScreenHeight * 0.15
        
        self.infoImageWidth = self.buttonWidth
        self.infoImageHeight = self.buttonWidth*3

        
        #load UI objects
        self.gameScreen =pygame.Rect(0,self.loadBarHeight,self.gameScreenWidth,
                                            self.screenSize)
        self.loadBar = pygame.Rect(0,0,self.gameScreenWidth,self.loadBarHeight)
        self.infoBar = pygame.Rect(self.gameScreenWidth,0,self.screenSize,
                                    self.screenSize)
        self.loadButton =pygame.Rect(self.loadButtonX1,self.loadButtonY1,
                                        self.buttonWidth,self.buttonHeight)

        self.loadText = self.gamefont.render("LOAD",True,self.black)

        self.clearButton =pygame.Rect(self.clearButtonX1,self.clearButtonY1,
                                        self.buttonWidth,self.buttonHeight)

        self.clearText = self.gamefont.render("CLEAR",True,self.black)

        self.helpButton =pygame.Rect(self.helpButtonX1,self.helpButtonY1,
                                        self.buttonWidth,self.buttonHeight)

        self.helpText = self.gamefont.render("HELP",True,self.black)

        self.infoText = self.gamefont.render("No selection",True,self.black)
        self.infoImage = pygame.image.load("assets/images/empty.png")
        self.infoImageBox = pygame.Rect(self.infoImageX,self.infoImageY,self.infoImageWidth,self.infoImageHeight)
        
        self.gameOverText = self.gamefont.render("GAMEOVER",True,self.white)

    def loadFASTA(self): 
        fileName = filedialog.askopenfilename() #prompts user to specify file

        if ("fasta" in fileName.lower()):
            self.gameFASTA = FASTA(fileName)
            sequence = self.gameFASTA.getSequence()
            start = random.randint(0,len(sequence)-6) #takes a fragment
            self.gameSequence = sequence[start:start+5]
            while("stop" in self.gameSequence):
                self.gameSequence = sequence[start:start+5]
            self.fileLoaded=True
            self.unsolvedSequence = copy.copy(self.gameSequence)
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
        self.bondLength = dx * 2.5
        currentY = self.gameScreenHeight/2
        for aminoAcid in self.gameSequence:
            aminoAcid.x = currentX
            aminoAcid.y = currentY
            aminoAcid.particle = Particle(str(aminoAcid),aminoAcid.abbr,
                                            aminoAcid.color,aminoAcid.x,
                                            aminoAcid.y,aminoAcid.r)
            currentX+=dx
    
    def helpScreen(self):
        self.helpText = self.gamefont.render("HELP INFORMATION",True,self.black)
        helpString = """NOTE: At the moment, the help function is still under development. Until eventually incorporated into the user interface, this feature will remain inthe print console. 
        Objective: Solve the puzzle by moving each amino acid into a satisfactory conformation.
        GREEN- nonpolar amino acids, prefer to be near the center of the protein due to hydrophobic interaction with water
        RED/BLUE - positive/negatively charged amino acids, prefer to be far from like charges
        PURPLE - polar amino acids, prefer to be away from the center of the protein due to hydrophilic interaction with water
        YELLOW - special case containing a sulfur atom, prefers to be near others with sulfur (if exists). Also nonpolar.

        The amino acids that have not yet reached their optimal conformation are listed in the console, along with feedback. 
        Load an example file, then click and drag each amino acid to move it. Keep in mind that the peptide bonds are of a fixed length."""
        print(helpString)        
        #self.screen.blit(self.helpText,self.gameScreen)
    
    def checkLegal(self,aminoAcid):
        print("unsolved:", end="")  #add to the real UI
        print(self.unsolvedSequence)
        otherAminoAcids = copy.copy(self.getGameSequence())
        otherAminoAcids.remove(aminoAcid)
        totalX = 0
        totalY = 0
        
        for otherAminoAcid in otherAminoAcids:
            totalX+=otherAminoAcid.x
            totalY+=otherAminoAcid.y
            if (aminoAcid.particle.boundingBox.colliderect(otherAminoAcid.particle.boundingBox)):
                print("Collision with other amino acid!")
                return False
        centerX = totalX/len(otherAminoAcids)
        centerY = totalY/len(otherAminoAcids)

        sidechain = aminoAcid.sidechain
        threshold = (self.gameScreenWidth/7)

        if (sidechain.sulfide):
            for other in otherAminoAcids:
                if(other.sidechain.sulfide):
                    distance = abs(other.x-aminoAcid.particle.x)
                    if(distance > threshold/2):
                        print("sulfide containing amino acids should be near other sulfide containing amino acids")
                        return self.unsatisfiedInteraction(aminoAcid, other)
                    else:
                        return self.satisfiedInteraction(aminoAcid, other)

        if (sidechain.hbond):
            for other in otherAminoAcids:
                if(other.sidechain.hbond):
                    distance = abs(other.x-aminoAcid.particle.x)
                    if(distance > threshold/2):
                        print("hbonding amino acids should be near other hbonding amino acids")
                        return self.unsatisfiedInteraction(aminoAcid, other)
                    else:
                        return self.satisfiedInteraction(aminoAcid, other)


        elif (not sidechain.polar):
            if (abs(aminoAcid.particle.x-centerX)>threshold or \
                abs(aminoAcid.particle.y-centerY)>threshold):
                print("nonpolar amino acid should be near the center")
                return self.unsatisfiedLocation(aminoAcid)
            else:
                return self.satisfiedLocation(aminoAcid)

        elif (sidechain.polar or sidechain.charge!=0):
            if (abs(aminoAcid.particle.x-centerX)<threshold/4 or \
                abs(aminoAcid.particle.y-centerY)<threshold/4):
                print("polar or charged amino acid should be near the edges")
                aminoAcid.satisfied = False
                return self.unsatisfiedLocation(aminoAcid)
            else:
                return self.satisfiedLocation(aminoAcid)

        elif(sidechain.charge==1):
            for other in otherAminoAcids:
                distance = abs(otherAminoAcid.x-aminoAcid.particle.x)
                if (other.sidechain.charge==-1):
                    if (distance < threshold*1.5):
                        print("positive charged amino acid cannot be near other positive charges")
                        return self.unsatisfiedInteraction(aminoAcid,other)
                    else:
                        #possibly animate
                        pass
                elif (other.sidechain.charge==1):
                    if (distance < threshold*1.5):
                        return self.satisfiedInteraction(aminoAcid, other)
     
        elif(sidechain.charge==-1):
            for other in otherAminoAcids:
                distance = abs(other.x-aminoAcid.particle.x)
                if (other.sidechain.charge==1):
                    if (distance < threshold*1.5):
                        print("negative charged amino acid cannot be near other negative charges")
                        return self.unsatisfiedInteraction(aminoAcid,other)
                    else:
                        #possibly animate
                        pass
                elif (other.sidechain.charge==1):
                    if (distance < threshold*1.5):
                        return self.satisfiedInteraction(aminoAcid, other)
    
        return True

    def unsatisfiedInteraction(self,aminoAcid,other):
        if (other in aminoAcid.associates):
            aminoAcid.associates.remove(other)
        if (aminoAcid in other.associates):
            other.associates.remove(aminoAcid)
        aminoAcid.satisfied = False
        if (aminoAcid not in self.unsolvedSequence):
            self.unsolvedSequence.append(aminoAcid)
        print("interaction bad") #include other

        return False 
    
    def satisfiedInteraction(self,aminoAcid, other):
        aminoAcid.associates.add(other)
        other.associates.add(aminoAcid)
        aminoAcid.satisfied = True
        other.satisfied = True
        if (aminoAcid in self.unsolvedSequence):
            self.unsolvedSequence.remove(aminoAcid)
        if (other in self.unsolvedSequence):
            self.unsolvedSequence.remove(other)
        print("interaction gud")
        return True
    
    def unsatisfiedLocation(self,aminoAcid):
        aminoAcid.satisfied = False
        if (aminoAcid not in self.unsolvedSequence):
            self.unsolvedSequence.append(aminoAcid)
        print("location bad")
        return False 
    
    def satisfiedLocation(self,aminoAcid):
        aminoAcid.satisfied = True
        if (aminoAcid in self.unsolvedSequence):
            self.unsolvedSequence.remove(aminoAcid)
        print("location gud")
        return True

    def checkAssociates(self,aminoAcid):
        otherAminoAcids = copy.copy(self.getGameSequence())
        otherAminoAcids.remove(aminoAcid)
        sidechain = aminoAcid.sidechain
        threshold = (self.gameScreenWidth/7)
                    

        for associate in aminoAcid.associates:     
            distance = abs(associate.x-aminoAcid.particle.x)
            if(distance > threshold/3):
                print("interaction broken")
                #aminoAcid.satisfied = False
                self.unsolvedSequence.append(aminoAcid)
                aminoAcid.associates.remove(associate)
                associate.associates.remove(aminoAcid)

    def checkForWin(self):
        for aminoAcid in self.gameSequence:
            if not aminoAcid.satisfied:
                return False
        print("SOLVED!!!")
        return True

    def undoMove(self,aminoAcid,iterations=10):
        
        dx = (aminoAcid.particle.x - aminoAcid.x) / 10
        dy = (aminoAcid.particle.y - aminoAcid.y) / 10

        for i in range(iterations):
            pygame.draw.circle(self.screen,self.white,(int(aminoAcid.particle.x),int(aminoAcid.particle.y)),aminoAcid.r)
            index = self.gameSequence.index(aminoAcid)
            if (index!=0):
                prevElement = self.gameSequence[index-1]
                prevX = prevElement.x
                prevY = prevElement.y
                pygame.draw.line(self.screen,self.white,(prevX,prevY),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)
            if (index!=len(self.gameSequence)-1):
                nextElement = self.gameSequence[index+1]
                nextX = nextElement.x
                nextY = nextElement.y
                pygame.draw.line(self.screen,self.white,(nextX,nextY),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)
            aminoAcid.particle.x -= dx
            aminoAcid.particle.y -= dy
            self.drawSequence()

            pygame.time.delay(15)
    
    def drawSequence(self):
        if (self.fileLoaded):
            prevX = None
            prevY = None
            for aminoAcid in self.getGameSequence():
                if (prevX!=None):
                    pygame.draw.line(self.screen,self.black,(prevX,prevY),
                                            (aminoAcid.particle.x,
                                            aminoAcid.particle.y),10)
                prevX = aminoAcid.particle.x
                prevY = aminoAcid.particle.y
                aminoAcid.draw(self)
                #print(aminoAcid.associates)
                for associate in aminoAcid.associates: #draw interactions
                    if (aminoAcid.sidechain.sulfide and associate.sidechain.sulfide):
                        pygame.draw.line(self.screen,self.yellow,(associate.particle.x,associate.particle.y),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)
                    elif (aminoAcid.sidechain.hbond and associate.sidechain.hbond):
                        pygame.draw.line(self.screen,self.blue,(associate.particle.x,associate.particle.y),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)
                    elif (abs(aminoAcid.sidechain.charge - other.sidechain.charge)==1):
                        pygame.draw.line(self.screen,self.purple,(associate.particle.x,associate.particle.y),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)                            

            pygame.display.update()
    
    def displayInformation(self,aminoAcid):
        self.infoText = self.gamefont.render(aminoAcid.name,True,self.black)  
        path = "assets/images/" +aminoAcid.name.lower()+".png"
        self.infoImage = pygame.image.load(path)
        self.infoImage = pygame.transform.scale(self.infoImage,[200,200])

    def resetSelection(self):
        self.infoText = self.gamefont.render("No selection",True,self.black) 
        self.infoImage= pygame.image.load("assets/images/empty.png")

    def gameOver(self):
        self.done = True
        self.screen.fill((0,255,0))
        self.screen.blit(self.gameOverText,self.gameScreen.center)
        
def main():
    game = Game(10)
    #game.FASTAtest("assets/data/lab.fasta")

    while not game.done:
        pygame.time.delay(game.delay)
        game.screen.fill((255,255,255))
        pygame.draw.rect(game.screen,game.black,game.gameScreen, 10)

        pygame.draw.rect(game.screen,game.black,game.loadBar,10)
        pygame.draw.rect(game.screen,game.black,game.infoBar,10)

        pygame.draw.rect(game.screen,game.green,game.loadButton)
        game.screen.blit(game.loadText,game.loadButton)

        pygame.draw.rect(game.screen,game.red,game.clearButton)
        game.screen.blit(game.clearText,game.clearButton)

        pygame.draw.rect(game.screen,game.blue,game.helpButton)
        game.screen.blit(game.helpText,game.helpButton)

        game.screen.blit(game.infoText,game.infoBar)
            
        game.screen.blit(game.infoImage, game.infoImageBox.topleft)
       
        game.drawSequence()
        
        

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
                elif (game.clearButton.collidepoint(mousePos) and left):
                    game.gameSequence = None
                    game.gameFASTA = None
                    game.fileLoaded = False
      
                elif (game.helpButton.collidepoint(mousePos) and left):
                    game.helpScreen()
    
                elif(game.screen.get_at(mousePos)!=game.white \
                    and mouseX<game.gameScreenWidth \
                    and mouseY>game.loadBarHeight \
                    and mouseY<game.gameScreenHeight):
        
                    for aminoAcid in game.getGameSequence():
                        sqx = (mouseX - aminoAcid.particle.x)**2
                        sqy = (mouseY - aminoAcid.particle.y)**2

                        if (math.sqrt(sqx + sqy) < aminoAcid.r):

                            aminoAcid.particle.clicked=True
                            game.displayInformation(aminoAcid)
                            
                else:
                    game.resetSelection()

            elif (event.type == pygame.MOUSEBUTTONUP):
                mouseX, mouseY = pygame.mouse.get_pos()
                if (game.fileLoaded):
                    for aminoAcid in game.getGameSequence():
                        if (aminoAcid.particle.clicked):                      
                            if(game.checkLegal(aminoAcid) \
                                and mouseX>0 and mouseY>game.loadBarHeight \
                                and mouseX<game.gameScreenWidth and mouseY < game.gameScreenHeight):
                                aminoAcid.x = aminoAcid.particle.x
                                aminoAcid.y = aminoAcid.particle.y
                                game.checkAssociates(aminoAcid)
                                if(game.checkForWin()):
                                    game.gameOver()
                            else:
                               
                                game.undoMove(aminoAcid)
    
                            aminoAcid.particle.clicked=False
                            
                
            elif (event.type == pygame.MOUSEMOTION):
                if (game.fileLoaded):
                    mousePos = event.pos
                    mouseX,mouseY = mousePos
                    seq = game.getGameSequence()
                    brokenLink = None
                    for i in range(len(seq)):
                        aminoAcid = seq[i]
                        if (aminoAcid.particle.clicked):
                            if (i!=0):
                                otherElement = game.gameSequence[i-1]
                                otherX = otherElement.x
                                otherY = otherElement.y
                                
                        
                            elif (i!=len(game.gameSequence)-1):
                                otherElement = game.gameSequence[i+1]
                                otherX = otherElement.x
                                otherY = otherElement.y

                            sqx = (otherX - aminoAcid.particle.x)**2
                            sqy = (otherY - aminoAcid.particle.y)**2

                            distance = math.sqrt(sqx+sqy)
                            if (distance>game.bondLength):
                                
                                aminoAcid.particle.clicked=False
                                brokenLink=aminoAcid
                            else:
                                
                                aminoAcid.particle.x = mouseX
                                aminoAcid.particle.y = mouseY
                    if (brokenLink!=None):
                        #tell the user what they did
                        game.undoMove(brokenLink,iterations=7)
                        brokenLink=None
                       
        pygame.display.update()
    pygame.quit()


if __name__ == '__main__':
    main()

