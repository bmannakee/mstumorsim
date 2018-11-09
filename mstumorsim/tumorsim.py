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


class Cell:
       # TODO: add bozic parameters Should mutations be their own object. With ID and spectrum? '''
        def __init__(self,mut_rate,parent,is_empty=False,spectrum=None,mutations=None, seq_type = 'exome'):
            self.is_empty = is_empty
            self.rep_rate = 2 # Cells just divide
            self.dormant = False
            self.mr = mut_rate
            self.seq_type = seq_type
            if mutations == None:
                self.mutations = []
            else:
                self.mutations = mutations
            self.id = str(uuid4())
            self.parent_id = parent
            self.spectrum = spectrum # Each cell has a spectrum and division generates mutations from that spectrum

        def reproduce(self,new_sigs=None, force=False):
         # This is where we reproduce TODO: add bozic parameters '''
           # a d/b ratio of .72 implies 72 deaths per 100
           # births. bernoulli(.58) gets that. 
           # force = True forces reproduction, useful for starting new tree
            daughters = [] # Daughters has to be an iterable even if empty
            if self.seq_type == "exome":
                num_muts = 3
            else:
                num_muts = 300
            if self.is_empty:
                ## Generating an empty tree
                if force:
                    rep = 1
                else:
                    rep = bernoulli.rvs(size=1,p=0.51)

                if (rep==1 and not self.dormant):
                    for i in range(self.rep_rate):
                        daughter = Cell(mut_rate = self.mr,parent = self.id, is_empty = True, seq_type = self.seq_type)
                        daughters.append(daughter)
            else:
                ## Not generating an empty tree
                tmp_spec = Spectrum(new_sigs)
                if (np.array_equal(tmp_spec.sigs,self.spectrum.sigs) and not force):
                    rep = bernoulli.rvs(size=1,p=0.51) #  d/1-d=.72 gives, b=.58, d=.42
                else:
                    # if sigs change, force at least the initial reproduction
                    rep = 1
                if (rep==1 and not self.dormant):
                    for i in range(self.rep_rate):
                        current_mutations = copy.deepcopy(self.mutations)
                        current_mutations.extend([Mutation(tmp_spec.spectrum) for x in list(range(num_muts))])
                        daughter = Cell(mut_rate = self.mr,parent = self.id,is_empty = False, spectrum = tmp_spec,mutations =current_mutations)
                        daughters.append(daughter)
                else:
                    self.dormant = True # Safety issues to work out here. should be private?
            return(daughters)
        
        def get_mutations(self):
            return([(m.id,m.mutation_class) for m in self.mutations])

        def get_mutation_ids(self):
            return(m.id for m in self.mutations)
            
        def get_sigs(self):
            return self.spectrum.sigs

class SNVtree:
   # A tree that does stuff. More documentation please '''

    def __init__(self,num_cells,empty=False,sigs=None,timepoints=None,seq_type='exome'):
        '''
        num_cells: The total number of cells in the tree at completion
        empty: Bool - True to create a tree of empty cells, False to simulate cells with mutations.
        sigs: List of integers (1-30) specifying the COSMIC mutational signatures to include in the tumor
        timepoints: List of floats the same length as "sigs" at which the signature should become active.
                    timepoint 0 for signatures present at initiation, .25 for sigs appearing after 25% of
                    cells have accumulated, and so on.
        seq_type: one of "exome" or "wgs". Determines the mutations selected and the number of mutations
                    per cell division
        '''
        self.make_empty = empty
        if not self.make_empty:
            self.input_sigs = np.array(sigs)
            self.input_timepoints = np.array(timepoints)
            self.init_sigs = self.input_sigs[np.where(self.input_timepoints == 0.)]
            self.init_spectrum = Spectrum(self.init_sigs)
            self.add_sigs = self.input_sigs[np.where(self.input_timepoints > 0.)]
            self.add_timepoints = self.input_timepoints[np.where(self.input_timepoints > 0.)]
            self.next_timepoint_index = 0
            self.seq_type = seq_type
        else:
            self.input_sigs = None
            self.input_timepoints = None
            self.init_sigs = None
            self.init_spectrum = None
            self.add_sigs = None
            self.add_timepoints = None
            self.next_timepoint_index = None
        self.n = num_cells
        self.queue = deque()
        # Begin with a single cell with parent = 1. during the run force the first 100 cells to reproduce
        if self.make_empty:
            initial_cell = [Cell(mut_rate=2,parent=1,is_empty=True)]
        else:
            initial_cell = [Cell(seq_type = self.seq_type, mut_rate = 2, parent = 1,spectrum = self.init_spectrum, \
                mutations = [Mutation(self.init_spectrum.spectrum),Mutation(self.init_spectrum.spectrum)])] 
                            
        self.queue.extend(initial_cell)



    def run(self):
        if self.make_empty:
            while len(self.queue) < 10:
                # Get the tree started with 10 cells that are gauranteed to reproduce
                c = self.queue.popleft()
                if c.dormant:
                    self.queue.append(c)
                    continue
                else:
                    new_cells = c.reproduce(force=True)
                    self.queue.append(c)
                    self.queue.extend(new_cells)

            while len(self.queue) < self.n:
                c = self.queue.popleft()
                if c.dormant:
                    self.queue.append(c)   
                else:
                    new_cells = c.reproduce()
                    self.queue.append(c)
                    self.queue.extend(new_cells)
        else:
            while len(self.queue) < 4:
                # Get the tree started with 4 cells that are gauranteed to reproduce
                c = self.queue.popleft()
                if c.dormant:
                    self.queue.append(c)
                    continue
                else:
                    new_sigs = c.get_sigs()
                    new_cells = c.reproduce(new_sigs = new_sigs, force = True)
                    self.queue.append(c)
                    self.queue.extend(new_cells)
                # Get the tree started with 10 cells that are gauranteed to reproduce
            while len(self.queue) < self.n:
                if 98 < len(self.queue) < 100:
                    print('made it to 10')
                c = self.queue.popleft()
                # test for dormancy here to speed things up
                if c.dormant:   
                    self.queue.append(c)
                    continue
                if self._test_for_new_spectrum(len(self.queue)):
                    # Awkward, we already incremented timepoint index
                    new_sigs = np.append(c.get_sigs(),self.add_sigs[self.next_timepoint_index - 1])  
                else:
                    new_sigs = c.get_sigs()
                new_cells = c.reproduce(new_sigs = new_sigs, force = False)
                self.queue.append(c) # Push current cell back onto the end of the queue. If cell is marked dormant it will never reproduce
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
            vaf = float(v)/num_cells
            if vaf < min_vaf:
                continue
            mcount.append((k,vaf))
        return mcount

    def get_graph(self):
        g = nx.DiGraph()
        g.add_node(1,)
        nodes = self.get_cells()
        for node in nodes:
            g.add_node(node.id,cell=node)
        edges = [(node.parent_id,node.id) for node in nodes]
        g.add_edges_from(edges)
        return(g)

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
        fr = fr[["chr","start","end","vaf","alt"]]
        if chr_style == "UCSC":
            fr["chr"] = "chr" + fr["chr"]
        # The sampling scheme is not guaranteed to produce unique mutations
        # Dropping is not the best answer
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
