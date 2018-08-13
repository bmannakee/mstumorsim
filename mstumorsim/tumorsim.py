from collections import deque
from uuid import uuid4
from scipy.stats import bernoulli, multinomial
import numpy as np
import networkx as nx


class Mutation:
    def __init__(self,spectrum):
        self.id = str(uuid4())
        # multinomial returns an array with a 1 in the bin that holds the returned value
        # and 0 everywhere else. TODO: 1-indexed for further processing?
        self.mutation_class = np.argmax(multinomial.rvs(1,spectrum))


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

        def reproduce(self):
         # This is where we reproduce TODO: add bozic parameters '''
           # a d/b ratio of .72 implies 72 deaths per 100
           # births. bernoulli(.58) gets that. 
            daughters = []
            rep = bernoulli.rvs(size=1,p=0.58) #  d/1-d=.72 gives, b=.58, d=.42
            if (rep==1 and not self.dormant):
                for i in range(self.rep_rate):
                    new_mutations = [Mutation(self.spectrum) for x in range(self.mr)]
                    daughter = Cell(mut_rate = self.mr,parent = self.id,spectrum = self.spectrum,mutations = new_mutations)
                    daughters.append(daughter)
            else:
                self.dormant = True # Safety issues to work out here. should be private?
                return(daughters)
            return(daughters)
        
        def get_mutations(self):
            return([(m.id,m.mutation_class) for m in self.mutations])

class SNVtree:
   # A tree that does stuff. More documentation please '''

    def __init__(self,num_cells,spectrum):
        # TODO: add Bozic parameters '''
        self.spectrum = spectrum
        # QA for spectrum input
        if not (len(spectrum)==96 and sum(spectrum)==1):
           print("ERROR: spectrum must have 96 entries summing to 1")
        self.n = num_cells
        self.queue = deque()
        # seed with enough cells to survive. All have parent 1, and 2 different mutations.
        # here 1 is like a stem cell. TODO: more elegant way to do this?
        initial_cells = [Cell(mut_rate = 2, parent = 1,spectrum = self.spectrum,mutations = [Mutation(spectrum),Mutation(spectrum)]) for i in range(10)] 
        self.queue.extend(initial_cells)



    def run(self):
        while len(self.queue) < self.n:
           c = self.queue.popleft()
           new_cells = c.reproduce()
           self.queue.append(c) # Push current cell back onto the end of the queue
           self.queue.extend(new_cells)


    def get_cells(self):
        return list(self.queue)
    

    def get_graph(self):
        tree = nx.Graph()
        tree.add_node(1,)
        nodes = self.get_cells()
        for node in nodes:
            tree.add_node(node.id,mutations = node.get_mutations())
        edges = [(node.parent_id,node.id) for node in nodes]
        tree.add_edges_from(edges)
        return(tree)

