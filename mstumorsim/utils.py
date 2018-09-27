import os
import pandas
import pickle

this_path,this_file = os.path.split(__file__)
SIG_FILE = os.path.join(this_path,'data','signatures_probabilities.txt')
MUT_PICKLE = os.path.join(this_path,'data','mutation_dict.pkl')
MUT_TSV = os.path.join(this_path,'data','mutations_with_contexts_GRCh38.tsv')

def load_signatures():
    sigs = pandas.read_table(SIG_FILE)
    sigs.columns = sigs.columns.str.replace(" ","_")
    return sigs

# Load signature frame ONCE, rather than several million times
SIGS_FR = load_signatures()

def get_spectrum(sigs):
    
    cols = [f'Signature_{sig}' for sig in sigs]
    spec_fr = SIGS_FR[cols]
    spec = spec_fr.sum(axis=1)
    spec = spec/sum(spec)
    return(spec)

# Load mutations into a dictionary where key is the mutation class
# and value is a list of mutations.
def load_mutations(from_pickle = True, mutation_file = MUT_TSV):
    if from_pickle:
        mdict = pickle.load(open(MUT_PICKLE, "rb"))
    else:
        # bring everything in as a string
        fr = pandas.read_table(mfile,  dtype = "object").dropna()
        # make a single field chr_start_stop_ref_alt
        fr['mutation'] = fr[['chr','start','end','ref','alt']].apply(lambda x: '_'.join(x), axis = 1)
        mdict = {k: list(v) for k,v in fr.groupby("context_id")["mutation"]}
        pickle.dump(mdict, open(MUT_PICKLE, "wb"))
    return(mdict)