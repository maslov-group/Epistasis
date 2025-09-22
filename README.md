# Protein Language Models Capture Structural and Functional Epistasis in a Zero-Shot Setting
Protein language models (PLMs) learn from large collections of natural sequences and achieve striking success across prediction tasks. Yet, it remains unclear what biological principles underlie their representations. We use epistasis, the dependence of a mutation’s effect on its sequence context, as a lens to probe what PLMs capture about proteins. Comparing PLM-derived scores with deep mutational scanning data, we find that epistasis emerges naturally from pretrained models, without supervision on experimental fitness. Raw model scores align with residue–residue contacts, indicating that PLMs internalize structural proximity. Applying a nonlinear transformation to bring model outputs onto the experimental scale, however, shifts the signal toward functional couplings between distant sites. These findings show that PLMs capture both structural and functional dependencies from sequence data alone, and that epistasis provides a powerful window into the biological principles embedded in their representations.

## Data
Data used in this study can be found at this [link](https://drive.google.com/drive/folders/1P3x6xGwvQAecAfn79s7vDwCubmcCOuOZ?usp=sharing)


## Notebooks
To reproduce the figures in this paper, you can use the Jupyter notebooks in Figures.


## Obtaining PLM scores
Raw PLM scores can be obtained using the ProteinModel object in src/PLM/PLMscores.py file


## Nonlinear Transformation
Nonlinear tansformations performed in this study can be found in src/NonlinearTransform
