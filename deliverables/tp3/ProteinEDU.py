#Project: ProteinEDU
#Author: Rohan Patel
#this file contains the biological inheritance hierarchy, the FASTA translator, the autosolver algorithm, and the game itself
#translates DNA to amino acids, creates interactable amino acid objects, and includes a puzzle autosolve feature
#A comprehensive list of every external resource used in this project can be found at tinyurl.com/tp3citations

#import statements
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
import tkinter


#biological inheritance hierarchy: conserves code

#AminoAcid superclass
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
    
    def addToChain(self): #performs peptide bonding by removing amine and carboxyl groups
        self.amine = None
        self.carboxyl = None
    
    def draw(self,game):
        self.particle.draw_particle(game)
        
class Particle(object): #interface between the amino acid object and the pygame user interaction
    def __init__(self,name,abbr,color, x, y ,r):
        self.x = x
        self.y = y
        self.color = color
        self.r = r
        self.clicked = False
        self.name = name
        self.abbr = abbr
   
    def draw_particle(self,game): #draws the amino acid
        screen = game.screen
        self.boundingBox=pygame.draw.circle(screen,self.color,(int(self.x),int(self.y)),self.r)
        self.textRect = pygame.Rect(self.x,self.y,
                                        self.r*2,self.r*2)
        font = pygame.font.Font('freesansbold.ttf', 30) 
        self.nameText = font.render(self.abbr,True,pygame.Color("black"))
        game.screen.blit(self.nameText,self.textRect)


#FunctionalGroup superclass
class FunctionalGroup(object): 
    def __init__(self,charge,sulfide,hbond,polar):
        self.charge = charge
        self.sulfide = sulfide
        self.hbond = hbond
        self.polar = polar

#three common functional groups inherit from Functional Group superclass
class AlphaCarbon(FunctionalGroup):
    def __init__(self):
        super().__init__(0,False,False,False)
class Amine(FunctionalGroup):
    def __init__(self):
        super().__init__(1,False,False,True)
class Carboxyl(FunctionalGroup):
    def __init__(self):
        super().__init__(-1,False,False,True)

#all the unique amino acids inherit from AminoAcid superclass
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

#unique sidechain functional groups inherit from FunctionalGroup superclass
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

#FASTA file (dna) to AminoAcid sequence translator
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

    def FASTAtranslate(self,path): #performs all steps of the central dogma, produces finished bonded sequence
        dnaSequence = self.readToString(path)
        rnaSequence = self.transcribe(dnaSequence)
        self.aminoacidsequence = self.translate(rnaSequence)
        self.bond()
        
    #from http://www.cs.cmu.edu/~112/notes/notes-strings.html#basicFileIO
    def readFile(self, path): 
        with open(path, "rt") as f:
            return f.read()

    def readToString(self,path):  #reads a FASTA file by lines
        raw = self.readFile(path)
        dnaSequence = ""
        for line in raw.splitlines():
            if (line[0]==">"):
                self.title = line[0][1:]
            else:
                dnaSequence += (line.upper())
        return dnaSequence
        
    def transcribe(self,dnaSequence): #performs biological transcription
        rnaSequence = ""
        rnaSequence = dnaSequence.replace("T", "U")
        return rnaSequence

    def translate(self,rnaSequence): #performs biological translation
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

    def bond(self): #bonds all amino acids except for the terminus
        for i in range(1, len(self.aminoacidsequence)-2):
            currentAminoAcid = self.aminoacidsequence[i]
            if (isinstance(currentAminoAcid,AminoAcid)):
                currentAminoAcid.addToChain()

    def makeCodonDict(self):  #since multiple codons code for the same amino acid, a multi key dictionary was used
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

#autosolver using backtracking:

# from http://www.cs.cmu.edu/~112/notes/notes-recursion-part2.html#genericBacktrackingSolver

##############################################
# Generic backtracking-based puzzle solver
#
# Subclass this class to solve your puzzle
# using backtracking.
#
# To see how useful backtracking is, run with checkConstraints=True
# and again with checkConstraints=False
# You will see the number of total states go up (probably by a lot).
##############################################

class BacktrackingPuzzleSolver(object):
    def solve(self):
        self.moves = [ ]
        self.states = set()
        # If checkConstraints is False, then do not check the backtracking
        # constraints as we go (so instead do an exhaustive search)

        # Be sure to set self.startArgs and self.startState in __init__
        self.solutionState = self.solveFromState(self.startState)

        return (self.moves, self.solutionState)

    def solveFromState(self, state):
        if state in self.states:
            # we have already seen this state, so skip it
            return None
        self.states.add(state)
        if self.isSolutionState(state):
            # we found a solution, so return it!
            return state
        else:
            for move in self.getLegalMoves(state):
                # 1. Apply the move
                childState = self.doMove(state, move)
                # 2. Verify the move satisfies the backtracking constraints
                #    (only proceed if so)
                if ((self.stateSatisfiesConstraints(childState))):
                    # 3. Add the move to our solution path (self.moves)
                    self.moves.append(move)
                    # 4. Try to recursively solve from this new state
                    result = self.solveFromState(childState)
                    # 5. If we solved it, then return the solution!
                    if result != None:
                        return result
                    # 6. Else we did not solve it, so backtrack and
                    #    remove the move from the solution path (self.moves)
                    self.moves.pop()
            return None

    # You have to implement these:

    def __init__(self):
        # Be sure to set self.startArgs and self.startState here
        pass

    def stateSatisfiesConstraints(self, state):
        # return True if the state satisfies the solution constraints so far
        raise NotImplementedError

    def isSolutionState(self, state):
        # return True if the state is a solution
        raise NotImplementedError

    def getLegalMoves(self, state):
        # return a list of the legal moves from this state (but not
        # taking the solution constraints into account)
        raise NotImplementedError

    def doMove(self, state, move):
        # return a new state that results from applying the given
        # move to the given state
        raise NotImplementedError

##############################################
# Generic State Class
#
# Subclass this with the state required by your problem.
# Note that this is a bit hacky with __eq__, __hash__, and __repr__
# (it's fine for 112, but after 112, you should take the time to
# write better class-specific versions of these)
##############################################

class State(object):
    def __eq__(self, other): return (other != None) and self.__dict__ == other.__dict__
    def __hash__(self): return hash(str(self.__dict__)) # hack but works even with lists
    def __repr__(self): return str(self.__dict__)


class ProteinState(State): #contains the placed and unplaced amino acids
    def __init__(self, placed, unplaced):
        self.placed = placed
        self.unplaced = unplaced

    def __repr__(self):
        return "Placed: "+str(self.placed)+", Unplaced: "+str(self.unplaced)
 

class ProteinSolver(BacktrackingPuzzleSolver):
    def __init__(self, game): 
        self.game = game
        self.L = self.game.getGameSequence()
      
        self.L[0].particle.x = game.gameScreenWidth//2 
        self.L[0].particle.y = game.gameScreenHeight//2
        self.L[0].satisfied = True

        placed = [self.L[0]]
        unplaced = self.L[1:]
        self.startState = ProteinState(placed,unplaced)


    def stateSatisfiesConstraints(self, state):
        lastPlaced = state.placed[len(state.placed)-1]
        otherPlaced = copy.copy(state.placed)
        otherPlaced.remove(lastPlaced)
        if (otherPlaced):
            return self.game.checkLegal(lastPlaced,otherPlaced)
        return True

    def isSolutionState(self, state):
        if (state.placed and (not state.unplaced)):
            for aminoAcid in state.placed:
                if (not aminoAcid.satisfied):
                    return False
            return True
        return False
  
    def getLegalMoves(self, state): 
        if (not state.unplaced):
            return []
        legalMoves = []
        rangeOfMotion = [0,-1,+1,-2,+2,-3,+3]
        random.shuffle(rangeOfMotion)
        for dx in rangeOfMotion:
            for dy in rangeOfMotion:
                newAminoAcid = state.unplaced[0]
                lastPlaced = state.placed[len(state.placed)-1]
                newX = lastPlaced.particle.x + (dx*(self.game.bondLength/6))
                newY = newAminoAcid.particle.y = lastPlaced.particle.y + (dy*(self.game.bondLength/6))
                if (not (dx==0 and dy==0) \
                    and newX-newAminoAcid.r>0 \
                    and newY-newAminoAcid.r>self.game.loadBarHeight \
                    and newX-newAminoAcid.r<self.game.gameScreenWidth \
                    and newY+newAminoAcid.r<self.game.screenSize2):
                    legalMoves.append((dx,dy))
       
        return legalMoves


    def doMove(self, state, move):
        newAminoAcid = state.unplaced[0]
        dx,dy = move
        lastPlaced = state.placed[len(state.placed)-1]
        newAminoAcid.particle.x = lastPlaced.particle.x + (dx*(self.game.bondLength/6))
        newAminoAcid.particle.y = lastPlaced.particle.y + (dy*(self.game.bondLength/6))
        placed = copy.copy(state.placed)
        placed.append(newAminoAcid)
        unplaced = copy.copy(state.unplaced)
        unplaced.remove(newAminoAcid)
        newState = ProteinState(placed,unplaced)
        self.game.drawSequence()
        pygame.display.update()
        return newState


class Game(object): #main Game object for interaction with pygame
    def __init__(self,delay):
        #initializes game details
        self.gameFASTA = None
        self.gameSequence = None
        self.solvedSequence = []
        self.unsolvedSequence = []
        self.fileLoaded = False
        self.delay = delay
        self.done = False
        self.solutionDisplay = False
        self.aminoAcidInfo = ""
        self.advice = "Click and drag the amino acids to move them."
        self.issue = ""

        #sets up pygame environment
        os.environ['SDL_VIDEO_WINDOW_POS'] = "15,40"
        pygame.init()
        self.root = tkinter.Tk()
        self.width = self.root.winfo_screenwidth()
        self.height = self.root.winfo_screenheight()
        self.screen = pygame.display.set_mode((int(self.width),int(self.height)),pygame.RESIZABLE)
        self.screen.fill((255,255,255))
        pygame.display.set_caption("ProteinEDU")

        #visual assets
        self.smallfont = pygame.font.Font('freesansbold.ttf', 15)
        self.gamefont = pygame.font.Font('freesansbold.ttf', 20) 
        self.namefont = pygame.font.Font('freesansbold.ttf', 30) 
        self.titlefont = pygame.font.Font('freesansbold.ttf', 50) 
        self.green = (0,255,0)
        self.black = (0,0,0)
        self.white = (255,255,255)
        self.red = (255,0,0)
        self.blue = (0,0,255)
        self.yellow = (220,220,0)
        self.purple = (255,0,255)
        self.gray = (128,128,128)
        self.darkRed = (128,0,0)
        self.screenSize,self.screenSize2=pygame.display.get_surface().get_size()

        #splash screen
        self.titleText = self.titlefont.render("ProteinEDU",True,self.purple)
        self.descriptionText = self.gamefont.render("An interactive protein-folding educational game that teaches the user the basics of primary structure amino acid interactions.",True,self.blue)
        self.continueText = self.smallfont.render("Click anywhere to continue.",True,self.darkRed)


        #help screen
        self.helpTitleText = self.titlefont.render("HELP SCREEN",True,self.blue)
        self.objectiveText = self.gamefont.render("Objective: Solve the puzzle by moving each amino acid into a satisfactory conformation.",True,self.black)
        self.explanationText = self.gamefont.render("To be satisfied each amino acid on the screen must be EITHER in a satisfactory location OR interact with another amino acid.",True,self.black)
        self.categoriesText = self.gamefont.render("The 20 basic amino acids have been separated into five types.",True,self.black)
        self.nonpolarGreen = self.gamefont.render("GREEN- nonpolar amino acids",True,self.green)
        self.nonpolarYellow = self.gamefont.render("YELLOW- nonpolar amino acids, contains a sulfide",True,self.yellow)
        self.chargedRed = self.gamefont.render("RED - positive charged amino acids",True,self.red)
        self.chargedBlue = self.gamefont.render("BLUE - negatively charged amino acids",True,self.blue)
        self.polarPurple = self.gamefont.render("PURPLE - polar amino acids",True,self.purple)

        self.locationRules = self.gamefont.render("To satisfy via location, follow the below rules.", True, self.gray)
        self.nonpolarRules = self.gamefont.render("-Nonpolars prefer to be near the center of the protein due to hydrophobic interaction with water.",True,self.black)
        self.polarLocationRules = self.gamefont.render("-Polars prefer to be away from the center of the protein due to hydrophilic interaction with water",True,self.black)
        self.chargedLocationRules = self.gamefont.render("-Positive/negatively charged amino acids are also hydrophilic, but additionally must not be next to an like charge.",True,self.black)

        self.interactionRules = self.gamefont.render("To satisfy via interaction (satisfies both amino acids), follow the below rules.", True, self.gray)
        self.disulfideBridgeRules = self.gamefont.render("-Sulfide containing amino acids can form disulfide bridges.", True, self.black)
        self.electrostaticRules = self.gamefont.render("-Charged amino acids can experience electrostatic attraction with an opposite charge.", True, self.black)
        self.hbondingRules = self.gamefont.render("-All polar and charged amino acids can hydrogen bond with each other.", True, self.black)
        self.goalText = self.gamefont.render("In other words, each amino acid must meet one of the above six criteria to be considered biochemically satisified.", True, self.red)
        self.goal2Text = self.gamefont.render("Load an example FASTA file, then click and drag each amino acid to move it. Keep in mind that the peptide bonds are of a fixed length.", True, self.black)
        self.forMoreInfoText = self.smallfont.render("This project offers a very simplified view of proteins. For more information about amino acid biochemistry, visit: https://tinyurl.com/applicableinformation", True,self.darkRed)
 

        #game runtime screen
        #initialize heights/widths based on window size
        self.loadBarHeight = self.screenSize*(0.10)
        self.loadBarWidth = self.screenSize*(0.70)
        self.gameScreenWidth = self.loadBarWidth
        self.gameScreenHeight = self.screenSize

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

        self.solveButtonX1 = self.loadBarWidth*0.7
        self.solveButtonY1 = self.loadButtonY1


        
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

        self.resetSelection()
        self.infoImageBox = pygame.Rect(self.infoImageX,self.infoImageY,self.infoImageWidth,self.infoImageHeight)

        self.aminoAcidInfoText = self.gamefont.render(self.aminoAcidInfo,True,self.black)
        
        self.gameOverText = self.titlefont.render("GAMEOVER",True,self.black)
        self.clickAnywhereText = self.gamefont.render("Click anywhere to end the game.",True,self.black)

        self.solveButton =pygame.Rect(self.solveButtonX1,self.solveButtonY1,
                                        self.buttonWidth,self.buttonHeight)

        self.solveText = self.gamefont.render("AUTOSOLVE",True,self.black)

        self.issueText = self.gamefont.render(self.issue,True,self.black)



    def loadFASTA(self): #prompts the user to select a FASTA file, and loads it using the FASTA class
        dataDir="assets/data"
        
        fileName = filedialog.askopenfilename(initialdir=dataDir) #prompts user to specify file

        if ("fasta" in fileName.lower()):
            self.gameFASTA = FASTA(fileName)
            sequence = self.gameFASTA.getSequence()
            start = random.randint(0,len(sequence)-6) #takes a random fragment
            self.gameSequence = sequence[start:start+5]
            while("stop" in self.gameSequence):
                start = random.randint(0,len(sequence)-6) 
                self.gameSequence = sequence[start:start+5]
            self.fileLoaded=True
            self.unsolvedSequence = copy.copy(self.gameSequence)
            return True

      
        else:
            self.issue = "Wrong file type. Please load a FASTA file."
            return False
        

    def FASTAtest(self,fileName): #test of FASTA translation
        testTranslator = FASTA(fileName)
        final = testTranslator.getSequence()
        module_manager.review()
        return final

    def getGameSequence(self):
        return self.gameSequence
    
    def initializeGameScreen(self): #initializes amino acids on the game screen
        dx = self.gameScreenWidth/(len(self.gameSequence)+1)
        currentX = dx
        self.bondLength = dx * 2.5
        currentY = self.gameScreenHeight/3
        for aminoAcid in self.gameSequence:
            aminoAcid.x = currentX
            aminoAcid.y = currentY
            aminoAcid.particle = Particle(str(aminoAcid),aminoAcid.abbr,
                                            aminoAcid.color,aminoAcid.x,
                                            aminoAcid.y,aminoAcid.r)
            currentX+=dx
    

    
    def checkLegal(self,aminoAcid,otherAminoAcids=None): #checks if amino acid is biochemically satisfied relative to others
        if (otherAminoAcids==None):
            otherAminoAcids = copy.copy(self.getGameSequence())
            otherAminoAcids.remove(aminoAcid)

        totalX = 0
        totalY = 0
        for otherAminoAcid in otherAminoAcids:
            totalX+=otherAminoAcid.particle.x
            totalY+=otherAminoAcid.particle.y
            if (aminoAcid.particle.boundingBox.colliderect(otherAminoAcid.particle.boundingBox)):
                self.advice = "Collision with other amino acid!"
                return False
        centerX = totalX/len(otherAminoAcids)
        centerY = totalY/len(otherAminoAcids)

        sidechain = aminoAcid.sidechain
        threshold = (self.gameScreenWidth/7)
  

        if (sidechain.sulfide):
            for other in otherAminoAcids:
                if(other.sidechain.sulfide):
                    distance = abs(other.particle.x-aminoAcid.particle.x)
                    if(distance > threshold/2):
                        self.advice = "sulfide containing amino acids should be near other sulfide containing amino acids"
                  
                    else:
                        return self.satisfiedInteraction(aminoAcid, other)

        if(sidechain.charge==1):
            for other in otherAminoAcids:
                distance = abs(other.particle.x-aminoAcid.particle.x)
                if (other.sidechain.charge==1):
                    if (distance < threshold*1.5):
                        self.advice = "positive charged amino acid cannot be near other positive charges"
                  
         
                elif (other.sidechain.charge==-1):
                    if (distance < threshold*1.5):
                  
                        return self.satisfiedInteraction(aminoAcid, other)
     
        if(sidechain.charge==-1):
        
            for other in otherAminoAcids:
             
                distance = abs(other.particle.x-aminoAcid.particle.x)
                if (other.sidechain.charge==-1):
                    if (distance < threshold*1.5):
                        self.advice = "negative charged amino acid cannot be near other negative charges"
         
                elif (other.sidechain.charge==1):
                    if (distance < threshold*1.5):
           
                        return self.satisfiedInteraction(aminoAcid, other)
    

        if (sidechain.hbond):
            for other in otherAminoAcids:
                if(other.sidechain.hbond):
                    distance = abs(other.particle.x-aminoAcid.particle.x)
                    if(distance > threshold/2):
                        self.advice = "hbonding amino acids should be near other hbonding amino acids"
                     
                    else:
                        return self.satisfiedInteraction(aminoAcid, other)


        if (not sidechain.polar):
            if (abs(aminoAcid.particle.x-centerX)>threshold or \
                abs(aminoAcid.particle.y-centerY)>threshold):
                self.advice = "nonpolar amino acid should be near the center of the protein"
                
            else:
                return self.satisfiedLocation(aminoAcid)

        if (sidechain.polar or sidechain.charge!=0):
            if (abs(aminoAcid.particle.x-centerX)<threshold/4 or \
                abs(aminoAcid.particle.y-centerY)<threshold/4):
                self.advice = "polar or charged amino acid should be near the outside of the protein"
                
            else:
                return self.satisfiedLocation(aminoAcid)
        
        return False

    
    def satisfiedInteraction(self,aminoAcid, other): #if satisfied via interaction

        aminoAcid.associates.add(other)
        other.associates.add(aminoAcid)
        aminoAcid.satisfied = True
        other.satisfied = True
        if (aminoAcid in self.unsolvedSequence):
            self.unsolvedSequence.remove(aminoAcid)
            self.solvedSequence.append(aminoAcid)

        if (other in self.unsolvedSequence):
            self.unsolvedSequence.remove(other)
            self.solvedSequence.append(other)

        return True

    def satisfiedLocation(self,aminoAcid): #if satisfied via location
        aminoAcid.satisfied = True
        if (aminoAcid in self.unsolvedSequence):
            self.unsolvedSequence.remove(aminoAcid)
            self.solvedSequence.append(aminoAcid)

        return True


    def checkForWin(self): #check if all amino acids are satisfied
        for aminoAcid in self.gameSequence:
            if not aminoAcid.satisfied:
                return False
        return True

    def undoMove(self,aminoAcid,iterations=10): #reverts amino acid to previous position if not satisfied
        
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
    
    def drawSequence(self): #draws the amino acid sequence on the game screen
        if (self.fileLoaded):
            prevX = None
            prevY = None 
            for aminoAcid in self.getGameSequence(): #draw each amino acid
                if (prevX!=None):
                    pygame.draw.line(self.screen,self.black,(prevX,prevY),
                                            (aminoAcid.particle.x,
                                            aminoAcid.particle.y),10)
                prevX = aminoAcid.particle.x
                prevY = aminoAcid.particle.y
                aminoAcid.draw(self)
         
                for associate in aminoAcid.associates: #draw interactions as lines
                    if (aminoAcid.sidechain.sulfide and associate.sidechain.sulfide):
                        pygame.draw.line(self.screen,self.yellow,(associate.particle.x,associate.particle.y),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)
                    elif (aminoAcid.sidechain.hbond and associate.sidechain.hbond):
                        pygame.draw.line(self.screen,self.blue,(associate.particle.x,associate.particle.y),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)
                    elif ((abs(aminoAcid.sidechain.charge) - abs(associate.sidechain.charge))==0):
                        pygame.draw.line(self.screen,self.purple,(associate.particle.x,associate.particle.y),
                                                (aminoAcid.particle.x,
                                                aminoAcid.particle.y),10)                            

            pygame.display.update() 
    
    def displayInformation(self,aminoAcid): #updates/displays selected amino acid information
        self.infoText = self.namefont.render(aminoAcid.name,True,self.black)  
        path = "assets/images/" +aminoAcid.name.lower()+".png"
        self.infoImage = pygame.image.load(path)
        self.infoImage = pygame.transform.scale(self.infoImage,[200,200])
        self.aminoAcidInfo = ""
        if (aminoAcid.sidechain.polar):
            self.aminoAcidInfo+="polar, "
        elif (aminoAcid.sidechain.charge==1):
            self.aminoAcidInfo+="positively charged, "
        elif (aminoAcid.sidechain.charge==1):
            self.aminoAcidInfo+="negatively charged, "
        else:
            self.aminoAcidInfo+="nonpolar, "
        if (aminoAcid.sidechain.hbond):
            self.aminoAcidInfo+="hydrogen bonding, "
        else:
            self.aminoAcidInfo+="no hydrogen bonding, "
        if (aminoAcid.sidechain.sulfide):
            self.aminoAcidInfo+="contains sulfide. "
        else:
            self.aminoAcidInfo+="no sulfide. "

        self.aminoAcidInfoText = self.gamefont.render(self.aminoAcidInfo,True,self.black)
        

    def resetSelection(self): 
        self.infoText = self.namefont.render("No selection",True,self.black) 
        self.infoImage= pygame.image.load("assets/images/empty.png")

    def gameOver(self):
        self.done = False
        self.solutionDisplay = True
        advice = "Game Over."

        
def main(): #main game loop
    game = Game(10)
    intro(game)
    while not game.done:
        #draw main game elements
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
        
        pygame.draw.rect(game.screen,game.yellow,game.solveButton)
        game.screen.blit(game.solveText,game.solveButton)
        

        sideBarWidth = game.gameScreenWidth*1.02
        game.screen.blit(game.infoText,(sideBarWidth,game.gameScreenHeight/20))
        game.screen.blit(game.infoImage, game.infoImageBox.topleft)

        game.aminoAcidInfoText = game.gamefont.render(game.aminoAcidInfo,True,game.black)
        game.screen.blit(game.aminoAcidInfoText,(sideBarWidth, 15*game.gameScreenHeight/40))

        game.adviceText = game.smallfont.render(game.advice,True,game.black)
        game.screen.blit(game.adviceText,(sideBarWidth, game.gameScreenHeight/2))

        game.issueText = game.gamefont.render(game.issue,True,game.black)
        game.screen.blit(game.issueText,(sideBarWidth, 8*game.gameScreenHeight/15))

        game.satisfiedAminoAcidsText = game.smallfont.render("Solved: "+ str(game.solvedSequence),True,game.black)
        game.unsatisfiedAminoAcidsText = game.smallfont.render("Unsolved: "+ str(game.unsolvedSequence),True,game.black)

        game.screen.blit(game.satisfiedAminoAcidsText,(sideBarWidth, 18*game.gameScreenHeight/40))
        game.screen.blit(game.unsatisfiedAminoAcidsText,(sideBarWidth, 19*game.gameScreenHeight/40))
        

        if (game.solutionDisplay): #if puzzle has been solved and solution is being viewed

            game.screen.blit(game.gameOverText,(game.width/3, 4*game.height/14))
            game.screen.blit(game.clickAnywhereText,(game.width/3, 5*game.height/14))
        
        game.drawSequence()

        for event in pygame.event.get():
            if (event.type == pygame.QUIT):
                game.done = True
            elif (event.type == pygame.MOUSEBUTTONDOWN):
                if (game.solutionDisplay):
                    game.done = True
                else:
                    mousePos = pygame.mouse.get_pos()
                    mouseX,mouseY = mousePos
                    left,right,middle=pygame.mouse.get_pressed()
                    if (game.loadButton.collidepoint(mousePos) and left): #load button clicked
                        if (game.loadFASTA()):
                            game.initializeGameScreen()
                            game.fileLoaded = True
                    elif (game.clearButton.collidepoint(mousePos) and left): #clear button clicked
                        game.gameSequence = None
                        game.gameFASTA = None
                        game.fileLoaded = False
        
                    elif (game.helpButton.collidepoint(mousePos) and left): #help button clicked
                        help(game)

                    elif (game.solveButton.collidepoint(mousePos) and left): #solve button clicked
                        
                        if (game.fileLoaded):
                            game.advice = "solving...."
                            solve(game)

                    elif(game.screen.get_at(mousePos)!=game.white \
                        and mouseX<game.gameScreenWidth \
                        and mouseY>game.loadBarHeight \
                        and mouseY<game.gameScreenHeight):
                        #if clicked somewhere on the screen    
                        for aminoAcid in game.getGameSequence():
                            sqx = (mouseX - aminoAcid.particle.x)**2
                            sqy = (mouseY - aminoAcid.particle.y)**2

                            if (math.sqrt(sqx + sqy) < aminoAcid.r): 
                                #determines if click is within a circle

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
                        game.advice = "Bond length stretched too far."
                        game.undoMove(brokenLink,iterations=5)
                        brokenLink=None
                       
        pygame.display.update()
   
    
    pygame.quit()

def intro(game): #intro (splash screen)
    introDone = False
    while (not introDone):
        #load intro screen static elements
        game.screen.fill((255,255,255))
        game.screen.blit(game.titleText,(game.screenSize/3,game.screenSize/12))
        game.screen.blit(game.descriptionText,(game.screenSize/12,game.screenSize/6))
        game.screen.blit(game.continueText,(game.screenSize/3,game.screenSize/120))
        for event in pygame.event.get():
            if (event.type == pygame.MOUSEBUTTONDOWN):
                introDone=True
                help(game)
            if (event.type == pygame.QUIT):
                introDone=True
                game.done = True
                
        pygame.display.update() 

def help(game): #help screen
    helpDone = False
    while (not helpDone):
        #load help screen static elements
        game.screen.fill((255,255,255))
        textX1 = game.screenSize/15
        textX2 = game.screenSize/5
        game.screen.blit(game.continueText,(game.screenSize/3,game.screenSize/120))
        game.screen.blit(game.helpTitleText,(game.screenSize/3,game.screenSize/40))
        game.screen.blit(game.objectiveText,(textX1,3*game.screenSize/40))
        game.screen.blit(game.explanationText,(textX1,4*game.screenSize/40))
        game.screen.blit(game.categoriesText,(textX1,5*game.screenSize/40))
        game.screen.blit(game.nonpolarGreen,(textX1,6*game.screenSize/40))
        game.screen.blit(game.nonpolarYellow,(textX1,7*game.screenSize/40))
        game.screen.blit(game.chargedRed,(textX1,8*game.screenSize/40))
        game.screen.blit(game.chargedBlue,(textX1,9*game.screenSize/40))
        game.screen.blit(game.polarPurple,(textX1,10*game.screenSize/40))
        game.screen.blit(game.locationRules,(textX1,11*game.screenSize/40))
        game.screen.blit(game.nonpolarRules,(textX2,12*game.screenSize/40))
        game.screen.blit(game.polarLocationRules,(textX2,13*game.screenSize/40))
        game.screen.blit(game.chargedLocationRules,(textX2,14*game.screenSize/40))
        game.screen.blit(game.interactionRules,(textX1,15*game.screenSize/40))
        game.screen.blit(game.disulfideBridgeRules,(textX2,16*game.screenSize/40))
        game.screen.blit(game.electrostaticRules,(textX2,17*game.screenSize/40))
        game.screen.blit(game.hbondingRules,(textX2,18*game.screenSize/40))
        game.screen.blit(game.goalText,(textX1,19*game.screenSize/40))
        game.screen.blit(game.goal2Text,(textX1,20*game.screenSize/40))
        game.screen.blit(game.forMoreInfoText,(textX1,22*game.screenSize/40))
      
        for event in pygame.event.get():
            if (event.type == pygame.MOUSEBUTTONDOWN):
                helpDone=True
            if (event.type == pygame.QUIT):
                helpDone=True
                game.done = True
                
        pygame.display.update() 

def solve(game): #wrapper method for the autosolver
    solver = ProteinSolver(game)
    if(solver.solve()==None):
        game.issue="Autosolver could not find a solution."
    game.checkLegal(game.gameSequence[len(game.gameSequence)-1])
    game.unsolvedSequence = []
    game.solvedSequence = game.gameSequence
    if(game.checkForWin()):
        game.gameOver()

if __name__ == '__main__':
    main()

