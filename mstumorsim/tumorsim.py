from collections import deque, Counter
from uuid import uuid4
from scipy.stats import bernoulli, multinomial
import numpy as np
import networkx as nx
import copy
from mstumorsim.utils import get_spectrum, load_mutations
import pandas
import random
class Spectrum:
    def __init__(self,sigs):
        self.sigs = sigs
        self.spectrum = get_spectrum(self.sigs)



class Mutation:
    def __init__(self,spectrum):
        self.id = str(uuid4())
        # multinomial returns an array with a 1 in the bin that holds the returned value
        # and 0 everywhere else. 1-indexed for further processing.
        self.mutation_class = np.argmax(multinomial.rvs(1,spectrum)) + 1
    
    def __repr__(self):
        return("mutation {} of type {}".format(self.id,self.mutation_class))


class Cell:
       # TODO: add bozic parameters Should mutations be their own object. With ID and spectrum? '''
        def __init__(self,parent,mu = 3e-6, spectrum = None, mutations = None, seq_type = 'exome'):
            self.rep_rate = 2 # Cells just divide
            self.mu = mu
            self.seq_type = seq_type
            if mutations == None:
                self.mutations = []
            else:
                self.mutations = mutations
            self.id = str(uuid4())
            self.parent_id = parent
            self.spectrum = spectrum # Each cell has a spectrum and division generates mutations from that spectrum

        def reproduce(self,new_sigs=None):
         # This is where we reproduce TODO: add bozic parameters '''
            daughters = [] 
            if self.seq_type == "exome":
                num_muts = 3e7 * self.mu 
            else:
                num_muts = 3e9 * self.mu
                ## Not generating an empty tree
            tmp_spec = Spectrum(new_sigs)
            for i in range(self.rep_rate):
                current_mutations = copy.deepcopy(self.mutations)
                current_mutations.extend([Mutation(tmp_spec.spectrum) for x in list(range(int(num_muts)))])
                daughter = Cell(mu = self.mu,parent = self.id, spectrum = tmp_spec,mutations =current_mutations)
                daughters.append(daughter)
            return(daughters)
        
        def get_mutations(self):
            return([(m.id,m.mutation_class) for m in self.mutations])

        def get_mutation_ids(self):
            return(m.id for m in self.mutations)
            
        def get_sigs(self):
            return self.spectrum.sigs
        
        def __repr__(self):
            return("A {} with signatures {} and mutations {}".format(self.__class__.__name__,self.spectrum,[mutation for mutation in self.get_mutations()]))


class SNVtree:
   # A tree that does stuff. More documentation please '''

    def __init__(self, num_cells = 2, mu = 3e-6, sigs = None, timepoints = None, seq_type = 'exome'):
        '''
        num_cells: The total number of cells in the tree at completion
        mu: The point mutation rate
        sigs: List of integers (1-30) specifying the COSMIC mutational signatures to include in the tumor
        timepoints: List of floats the same length as "sigs" at which the signature should become active.
                    timepoint 0 for signatures present at initiation, .25 for sigs appearing after 25% of
                    cells have accumulated, and so on.
        seq_type: one of "exome" or "wgs". Determines the mutations selected and the number of mutations
                    per cell division
        '''
        self.mu = mu
        self.input_sigs = np.array(sigs)
        # Begin with current sigs = input sigs.
        # When signatures change they change for the whole tree, not just the cell.
        # Doing it for a single cell will likely be nearly invisible, but it is possible.
        
        self.input_timepoints = np.array(timepoints)
        self.init_sigs = self.input_sigs[np.where(self.input_timepoints == 0.)]
        self.current_sigs = self.init_sigs
        self.init_spectrum = Spectrum(self.init_sigs)
        self.add_sigs = self.input_sigs[np.where(self.input_timepoints > 0.)]
        self.add_timepoints = self.input_timepoints[np.where(self.input_timepoints > 0.)]
        self.next_timepoint_index = 0
        self.seq_type = seq_type
        self.n = num_cells
        if self.seq_type == "exome":
            num_muts = 3e7 * self.mu 
        else:
            num_muts = 3e9 * self.mu 
        self.queue = deque()
        initial_cell = [Cell(mu = self.mu, seq_type = self.seq_type, parent = 1,spectrum = self.init_spectrum, \
            mutations = [Mutation(self.init_spectrum.spectrum) for x in list(range(int(num_muts)))])] 
                            
        self.queue.extend(initial_cell)



    def run(self):
        while len(self.queue) < self.n:
            c = self.queue.popleft()
            if self._test_for_new_spectrum(len(self.queue)):
                # Awkward, we already incremented timepoint index
                self.current_sigs = np.append(c.get_sigs(),self.add_sigs[self.next_timepoint_index - 1])  
            new_cells = c.reproduce(new_sigs = self.current_sigs)
            self.queue.extend(new_cells)


    def get_cells(self):
        return list(self.queue)

    def get_mutations(self):
        # Should return a list of tuples
        # (mutation_id, mutation_class, population_AF)
        cells = self.get_cells()
        mutations = []
        for cell in cells:
            mutations.extend(cell.get_mutations())
        return(mutations)
        
    def get_mutation_ids(self):
        # Should return a list of tuples
        # (mutation_id, mutation_class, population_AF)
        cells = self.get_cells()
        mutations = []
        for cell in cells:
            mutations.extend(cell.get_mutation_ids())
        return(mutations)
    
    def get_mutation_vafs(self,min_vaf = 0.01):
        mc = Counter(self.get_mutations())
        mcount = []
        num_cells = len(self.get_cells())
        for k,v in mc.items():
            vaf = float(v)/(2*num_cells)
            if vaf < min_vaf:
                continue
            mcount.append((k,vaf))
        return mcount


    def make_bed(self, chr_style = "UCSC",min_vaf = 0.01, from_pickle = True):
        # Gather all of the mutations from this tree that
        # have a vaf exceeding min_vaf, and write a bed
        # file with vaf in the format required by bamsurgeon
        # where chromosome names are in UCSC or NCBI style
        mutation_vafs = self.get_mutation_vafs(min_vaf = min_vaf)
        print('got mutations')
        fr = pandas.DataFrame.from_records(mutation_vafs)
        print('made data frame')
        fr[['mut_id','mut_class']] = fr[0].apply(pandas.Series)
        fr['vaf'] = fr[1]
        fr = fr[['mut_id','mut_class','vaf']]
        if self.seq_type == "exome":
            mdict = load_mutations(sig_type = 'exome', from_pickle = from_pickle)
        else:
            mdict = load_mutations(sig_type = 'wgs', from_pickle = from_pickle) 
        
        # lambda function that takes in a mutation class
        # and returns a single random mutation from that class
        get_mutation = lambda x: mdict[str(x)][random.sample(range(len(mdict[str(x)])),1)[0]]
        print("pulling mutations")
        fr["mutation"] = fr["mut_class"].apply(get_mutation)
        print("expanding mutations")
        fr = fr.join(fr["mutation"].str.split('_', expand = True))
        fr.rename(columns = {0:"chr", 1:"start", 2:"end", 3:"ref", 4:"alt"}, inplace = True)
        fr = fr[["chr","start","end","vaf","ref","alt"]]
        if chr_style == "UCSC":
            fr["chr"] = "chr" + fr["chr"]
        # The sampling scheme is not guaranteed to produce unique mutations
        # Dropping is not the best answer
        #print("NOT! dropping duplicates, don't use this with bamsurgeon")
        print("dropping duplicates")
        fr.drop_duplicates(inplace = True)
        return(fr)
    
    def write_bed(self,filename,chr_style = "UCSC" , min_vaf = 0.01, from_pickle = True):
        fr = self.make_bed(chr_style = chr_style, min_vaf = min_vaf, from_pickle = from_pickle)
        print("got bed frame")
        fr.to_csv(filename, sep = "\t", index = False, header = False)
        return 0

    def _test_for_new_spectrum(self, ncells):
        # We reproduce at rate 2, so there is no guarantee
        # that we will hit the add index exactly
        # Incrementing the next_timepoint_index here prevents getting more
        # than one cell with the new spectrum
        if self.next_timepoint_index == self.add_timepoints.size:
            # We have used up all of the timepoints
            return False
        low = ncells - 5
        high = ncells + 5
        timepoint = int(self.add_timepoints[self.next_timepoint_index] * self.n)
        
        if low <= timepoint <= high:
            self.next_timepoint_index += 1
            return True
        else:
            return False 
    
    def __repr__(self):
        return(cell for cell in self.get_cells())
