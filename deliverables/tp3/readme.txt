ProteinEDU

An interactive protein-folding educational game that teaches the user the basics of primary structure amino acid interactions. The user can load one of several FASTA files (dna sequences), which is then automatically translated to a polypeptide and displayed in the game window. These amino acids can be interacted with until they are all in the biochemically-correct configuration relative to each other, at which point the puzzle is solved.

Inspiration for this project was taken primarily from Foldit, a similarly styled game that involves folding proteins, which actually crowdsources the "work" of determining the structure of real proteins by presenting them to the user in the form of puzzles to solve. The results can then actually be directly used by scientists to determine the correct folding conformation of these proteins, making a real world impact by helping them create models used to develop new medicines to treat disease. Another program called EteRNA, developed right here at Carnegie Mellon, is another "game with a purpose" that enlists the user to help discover the optimal folding of RNA molecules.

This game is purely educational, with the "proteins" being modeled solely at the lowest level of biological complexity, the primary structure (like beads on a string). It uses the same puzzle styled approach as the previous two, but instead of collecting data for the good of science, the idea is to teach anyone (regardless of background biology knowledge) how different amino acids interact with one another, because those basic interactions (hydrogen bonds, disulfide bridge, et. cetera) are what govern the folding of all proteins in every living being on earth. The project is also similar to PyMOL, a software (developed partly in Python) designed to allow users to visualize high-level structures of complicated proteins. However, that software is tailored primarily towards those with a fair deal of experience in biology, as it is difficult to use otherwise. Hence, my game is a combination of both of these types, aiming to teach anyone at all who plays the game the basics of protein folding, all the while never once touching a textbook.


To run the game, first download and unzip the tp3 file. In it, there are two python files and one folder which contains all the game assets.

Install the pygame and multi_key-dict libraries, which can be found in the libraries folder, under assets. The project also uses random, tkinter, math, os, sys, and copy.

The project itself runs out of a single python file (ProteinEDU.py), and uses module_manager.py to verify that all modules are installed correctly. To run, simply run the ProteinEDU.py file. No shortcuts are used, and the rest should be evident from the game.

Citations:

A comprehensive list of every external resource used in this project can be found at tinyurl.com/tp3citations. 