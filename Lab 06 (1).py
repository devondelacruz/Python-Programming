#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
class FastAreader :
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence     
##

class FindUniques:
    def __init__(self, seqs, header):
        """Initializes the sequences, getting the header and sequences to be used. It also creates the lists that are used accross
        the methods"""
        self.seqs = seqs #gets header from other class and makes it assessible to other parts of assignment
        self.powSetList = [] #Initializes powSetList to be used in other parts
        self.dupSetList = [] #Initializes dupSetList to be used in other methods
        self.uniqueList = [] #Initializes uniqueSetList to be used in other methods
        self.essSetList = [] #Initializes essSetList to be used in other methods
        self.header = header # gets header from other class and makes it assessible to other parts of assignment
    def powerSet(self):
        """Gets the total sequence, and parses through each sequence separated by header, and iterates through it, creating
        a powerset with every combination of nucleotides adjacent to each other."""
        for i in (self.seqs):#Outer nested loop, iterates through each sequences within the sequence list
            indSet = set() #Sets individual power set to nothing
            dupSet = set() #sets duplicate power set to nothing
            for start in range(len(i)): #middle nested loop,iterates through character in sequence
                end = len(i) #establishes end for the while loop, and allows power set to decrease sequentially
                while end > start: #innermost loop and adds the strings of a certain character and adds to indset and dupset
                    if i[start:end] not in indSet:#if statement to add to duplicate set list if already in individual set, otherwise just add to individual set
                        indSet.add(i[start:end])
                    else:
                        dupSet.add(i[start:end])
                    end -= 1#Makes the next iteration a little shorter so that the powerset can be created
            self.powSetList.append(indSet)#Adds individual powerset to a powerset list to be used in other methods
            self.dupSetList.append(dupSet)#Adds duplicate powerset to a duplicate powerset list to be used in other methods
        return
    def uniques(self):
        """Gets the total powersets, and iterates through each of them, creating one set containing the duplicates in the powerset
        from the other class and the other powersets not being iterated through. Then it takes the created set and subtracts it from
        the set of the power list being iterated through, removing all duplicates from the set."""
        dupSubtract = set()
        for i in range(len(self.powSetList)): #Outer portion of nested loop, loops through each sequence in powerset list
            dupSubtract = (self.powSetList[:i] + self.powSetList[i+1:]) #Establishes dupsubtract being every other sequence in powerset aside from the one being iterated over
            dupSubSet = set().union(self.dupSetList[i]) #Adds the duplicate set list of the same sequence being iterated over and adds to subtract set
            for j in dupSubtract:# inner loop that adds each of the other sequences to the subtract set instead of a list
                dupSubSet.union(j)
            uniqueSet = self.powSetList[i].difference(dupSubSet) #Unique Set is created that subtracts the unnecessary strings from subtract set from the wanted set, making a new list in the process
            self.uniqueList.append(uniqueSet)#Adds unique sets to a list to be used in later methods
        return

    def essential(self):
        """Takes each unique, and determines if when a nucleotide on either end is taken off, if there is another
        copy, at any point in the power set, then it will be removed from essential, as it is not unique nor is it essential"""
        self.essList = []
        for uniSet in self.uniqueList:#Outermost nested loop
            nonEssSet = set()
            for string in uniSet:#Inner loop that takes string and subtracts one from either end of all sequences to check if it has duplicates
                start = 1
                end = 1
                left = string[start:len(string)]#Left checks if taking one from the left makes it unique
                right = string[0:len(string)-end]#'' but with right ''
                if left in uniSet:
                    nonEssSet.add(string)#Adds longer copy strings to non essential set, same with right below
                    start += 1
                if right in uniSet:
                    nonEssSet.add(string)
                    end += 1
            essSet = uniSet.difference(nonEssSet) #Subtracts non essential sets to make only the wanted sets
            self.essList.append(essSet)#Adds sets to respective list
        return(self.essList)    
    
    def output(self):
        """This method formats and returns the desired output. It gets the header and sequence and prints them sequentially,
        along with printing portions of the string if they can be found in the sequence."""
        for head in range(0,len(self.header)):#Iterates through each header
            print(self.header[head])#Prints each header sequentially
            print(self.seqs[head])#Prints respective sequence sequentially
            strTest = ''
            strTestTot = ''
            for character in self.seqs[head]: #Loop to make a longer string to check positions of set strings
                strTestTot +=(character) #Allows for position of set to be correct
                for char in range(0,len(strTestTot)):
                    strTest = strTestTot[char:len(strTestTot)]#gets and creates new smaller strings that match the powerset strings
                    if strTest in self.essList[head]:# If statement to check if in essential list to then print it
                        print(('.'*(len(strTestTot)-len(strTest))),end ='') #formats the dots, by comparing where the string is found and multiplying dots by that amount
                        print(strTest)#Adds the string from powerset to the output

        
########################################################################
# Main
# Here is the main program
# 
########################################################################

def main(inCL=None):
    '''Gets file and puts it into fastaReader, formatting it correctly, then passes to FindUniques class, printing the output '''
    fasta = FastAreader(inCL)
    headers, sequences = [],[]
    for head, seq in fasta.readFasta():
        head = head.replace(' ','')#replaces unwanted characters in header
        headers.append(head)
        seq = seq.replace('-','')#Replaces unwanted characters in sequence, same as two below
        seq = seq.replace('.','')
        seq = seq.replace('_','')
        sequences.append(seq)

    unique = FindUniques(sequences, headers)# puts sequences and headers into finduniques class and puts it into unique to be used quicker
    (unique.powerSet())#Initializes powerSet class
    (unique.uniques()) #Initializes uniques class
    (unique.essential()) #Initializes essential class
    return(unique.output())

if __name__ == "__main__":
    main('bos-tRNA.fa') #For any specific dataset, the input can be modified, takes fastA files however. 


# In[ ]:




