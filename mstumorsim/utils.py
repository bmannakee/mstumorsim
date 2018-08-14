import os
import pandas

this_path,this_file = os.path.split(__file__)
SIG_FILE = os.path.join(this_path,'data','signatures_probabilities.txt')

def load_signatures():
    sigs = pandas.read_table(SIG_FILE)
    sigs.columns = sigs.columns.str.replace(" ","_")
    return sigs

def get_spectrum(sigs):
    sigs_fr = load_signatures()
    cols = [f'Signature_{sig}' for sig in sigs]
    spec_fr = sigs_fr[cols]
    spec = spec_fr.sum(axis=1)
    spec = spec/sum(spec)
    return(spec)

