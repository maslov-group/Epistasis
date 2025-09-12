import pandas as pd
import esm
import matplotlib.pyplot as plt
import torch

from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import itertools
from typing import List, Tuple
import numpy as np
import pickle

import math
import os



class ProteinModel:
    """
    To try a model other than ESM, make a similar object to ProteinModel. The important thing is that it has a get_scores function
    """
    def __init__(self, model_location):
        self.model, self.alphabet, self.batch_converter = self.load_model(model_location)
        
    def load_model(self, model_location):
        model, alphabet = pretrained.load_model_and_alphabet(model_location)
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
            print("Transferred model to GPU")
        batch_converter = alphabet.get_batch_converter()
        return model, alphabet, batch_converter

    def get_scores(self, sequence, aas):
    
        """
        This method calculates mutation scores for a given protein sequence and a set of possible amino acid mutations.

        Args:
            sequence (str): The wild-type protein sequence.
            aas (str): A string containing all possible mutant amino acids.

        Returns:
            pd.DataFrame: A DataFrame containing mutation identifiers and their corresponding scores.
        """
        
        mutation_scores = []
        mutations = []
        logpaa = []
        aaenc = []
        logpwt = []
        wtenc = []
        
        data = [("protein1", sequence)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        
        with torch.no_grad():
            token_probs = torch.log_softmax(self.model(batch_tokens.cuda())["logits"], dim=-1)
        
        # Iterate over each position in the sequence
        for idx, wt in enumerate(sequence):
            wt_encoded = self.alphabet.get_idx(wt)
            tokenprobswt = token_probs[0, 1 + idx, wt_encoded]

            # Iterate over all possible mutant amino acids in the alphabet
            for mt in aas:  # aas = all possible amino acids
                if mt == wt:
                    continue  # Skip the case where the mutant is the same as wild-type

                
                mt_encoded = self.alphabet.get_idx(mt)
                tokenprobsaa = token_probs[0, 1 + idx, mt_encoded]
                
                # Compute the score similar to label_row
                score = tokenprobsaa - tokenprobswt
                
                # Append log prob of AA and encoding
                logpaa.append(tokenprobsaa.item())
                logpwt.append(tokenprobswt.item())
                aaenc.append(mt_encoded)
                wtenc.append(wt_encoded)

                # Append the mutation and score to the lists
                mutations.append(f"{wt}{idx + 1}{mt}")  # Add 1 to idx to make it 1-based index
                mutation_scores.append(score.item())
                
        return pd.DataFrame({"mutations": mutations, "score": mutation_scores, "log(p_AA)": logpaa, "AA": aaenc, "log(p_WT)": logpwt, "WT": wtenc})

def double_mutant_fitness_llm(dataset, model, wt_sequence, output_csv):
    """
    Calculates double mutant fitness scores using a cached LLM model.

    Args:
        dataset (str): Path to mutation CSV.
        model (str): Name of the ESM model.
        output_csv (str): Output file to save scores (appended every 1000 rows).

    Returns:
        None (writes to file incrementally)
    """

    df = pd.read_csv(dataset)
    pm = ProteinModel(model)

    # Cache for get_scores
    score_cache = {}
    def cached_get_scores(sequence, mutset):
        key = (sequence, mutset)
        if key not in score_cache:
            score_cache[key] = pm.get_scores(sequence, mutset)
        return score_cache[key]

    all_aas = "ACDEFGHIKLMNPQRSTVWY"
    _ = cached_get_scores(wt_sequence, all_aas)

    if os.path.exists(output_csv):
        os.remove(output_csv)

    batch = []

    for idx, row in df.iterrows():
        pos1 = row['Position 1']
        pos2 = row['Position 2']
        wt1, wt2 = row['Original AA 1'], row['Original AA 2']
        mut1, mut2 = row['Mutated AA 1'], row['Mutated AA 2']

        string1 = f"{wt1}{pos1}{mut1}"
        string2 = f"{wt2}{pos2}{mut2}"

        scs1 = cached_get_scores(wt_sequence, all_aas)

        mut1_row = scs1.loc[scs1['mutations'] == string1]
        mut2_row = scs1.loc[scs1['mutations'] == string2]
        s1 = mut1_row['score'].values[0] if not mut1_row.empty else 0
        s2 = mut2_row['score'].values[0] if not mut2_row.empty else 0
        
        # Apply mut1 at pos1, get mut2 score
        wt_21 = wt_sequence[:pos1 - 1] + mut1 + wt_sequence[pos1:]
        scs21 = cached_get_scores(wt_21, all_aas)
        mut21_score = scs21.loc[scs21['mutations'] == string2, 'score'].values
        s21 = mut21_score[0] if mut21_score.size else 0

        # Apply mut2 at pos2, get mut1 score
        wt_12 = wt_sequence[:pos2 - 1] + mut2 + wt_sequence[pos2:]
        scs12 = cached_get_scores(wt_12, all_aas)
        mut12_score = scs12.loc[scs12['mutations'] == string1, 'score'].values
        s12 = mut12_score[0] if mut12_score.size else 0

        # Final fitness scores
        path1 = s1 + s21
        path2 = s2 + s12
        avg_score = 0.5 * (path1 + path2)

        batch.append({
            'mut1': s1,
            'mut21': s21,
            'mut2': s2,
            'mut12': s12,
            'double_fitness_llm': avg_score
        })

        if (idx + 1) % 1000 == 0 or idx == len(df) - 1:
            df_batch = pd.DataFrame(batch)
            df_batch.to_csv(output_csv, mode='a', index=False, header=not os.path.exists(output_csv))
            batch.clear()
            print(f"{idx + 1} done")