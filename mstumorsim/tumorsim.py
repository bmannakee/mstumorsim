from collections import deque
from uuid import uuid4
from scipy.stats import bernoulli, multinomial
import numpy as np
import networkx as nx
from mstumorsim.utils import get_spectrum

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
        def __init__(self,mut_rate,parent,spectrum,mutations):
           self.rep_rate = 2 # Cells just divide
           self.dormant = False
           self.mr = mut_rate
           self.mutations = mutations
           self.id = str(uuid4())
           self.parent_id = parent
           self.spectrum = spectrum # Each cell has a spectrum and division generates mutations from that spectrum

        def reproduce(self,new_sigs):
         # This is where we reproduce TODO: add bozic parameters '''
           # a d/b ratio of .72 implies 72 deaths per 100
           # births. bernoulli(.58) gets that. 
            daughters = []
            tmp_spec = Spectrum(new_sigs)
            if np.array_equal(tmp_spec.sigs,self.spectrum.sigs):
                rep = bernoulli.rvs(size=1,p=0.58) #  d/1-d=.72 gives, b=.58, d=.42
            else:
                # if sigs change, force at least the initial reproduction
                rep = 1
            if (rep==1 and not self.dormant):
                for i in range(self.rep_rate):
                    new_mutations = [Mutation(tmp_spec.spectrum) for x in range(self.mr)]
                    daughter = Cell(mut_rate = self.mr,parent = self.id,spectrum = tmp_spec,mutations = new_mutations)
                    daughters.append(daughter)
            else:
                self.dormant = True # Safety issues to work out here. should be private?
                return(daughters)
            return(daughters)
        
        def get_mutations(self):
            return([(m.id,m.mutation_class) for m in self.mutations])

        def get_sigs(self):
            return self.spectrum.sigs


class SNVtree:
   # A tree that does stuff. More documentation please '''

    def __init__(self,num_cells,sigs,timepoints):
        '''
        num_cells: The total number of cells in the tree at completion
        sigs: List of integers (1-30) specifying the COSMIC mutational signatures to include in the tumor
        timepoints: List of floats the same length as "sigs" at which the signature should become active.
                    timepoint 0 for signatures present at initiation, .25 for sigs appearing after 25% of
                    cells have accumulated, and so on.
        '''
        self.input_sigs = np.array(sigs)
        self.input_timepoints = np.array(timepoints)
        self.init_sigs = self.input_sigs[np.where(self.input_timepoints == 0.)]
        self.init_spectrum = Spectrum(self.init_sigs)
        self.add_sigs = self.input_sigs[np.where(self.input_timepoints > 0.)]
        self.add_timepoints = self.input_timepoints[np.where(self.input_timepoints > 0.)]
        self.next_timepoint_index = 0
        self.n = num_cells
        self.queue = deque()
        # seed with enough cells to survive. All have parent 1, and 2 different mutations.
        # here 1 is like a stem cell. TODO: more elegant way to do this?
        initial_cells = [Cell(mut_rate = 2, parent = 1,spectrum = self.init_spectrum, \
                         mutations = [Mutation(self.init_spectrum.spectrum),Mutation(self.init_spectrum.spectrum)]) \
                         for i in range(10)] 
        self.queue.extend(initial_cells)



    def run(self):
        while len(self.queue) < self.n:
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
            new_cells = c.reproduce(new_sigs)
            self.queue.append(c) # Push current cell back onto the end of the queue. If cell is marked dormant it will never reproduce
            self.queue.extend(new_cells)


    def get_cells(self):
        return list(self.queue)
    

    def get_graph(self):
        tree = nx.DiGraph()
        tree.add_node(1,)
        nodes = self.get_cells()
        for node in nodes:
            tree.add_node(node.id,cell=node)
        edges = [(node.parent_id,node.id) for node in nodes]
        tree.add_edges_from(edges)
        return(tree)
        
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
