import random
from deepchem.utils import ScaffoldGenerator
from deepchem.utils.save import log
import torch as pt
from torch.autograd import Variable
import torch
import numpy as np
from rdkit import Chem
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from collections import OrderedDict
import deepchem as dc


def over_sample(train_set, train_label):
    pos_label = []
    neg_label = []
    for i in range(len(train_set)):
        if train_label[i] == 0:
            neg_label.append(i)
        else:
            pos_label.append(i)

    if len(pos_label) > len(neg_label):
        for i in range(len(pos_label) - len(neg_label)):
            id = random.sample(neg_label, 1)[0]
            train_set.append(train_set[id])
            train_label.append(0)
    else:
        for i in range(len(neg_label) - len(pos_label)):
            id = random.sample(pos_label, 1)[0]
            train_set.append(train_set[id])
            train_label.append(1)

    state = np.random.get_state()
    np.random.shuffle(train_set)
    np.random.set_state(state)
    np.random.shuffle(train_label)

    return train_set, train_label


def load_dataset(filepath, sample=True):

    train_features, train_labels = [], []
    val_features, val_labels = [], []
    test_features, test_labels = [], []

    df = pd.read_csv(filepath + 'train.csv')
    for i in range(len(df)):
        train_features.append(df.smiles[i])
        train_labels.append(df.activity[i])

    df = pd.read_csv(filepath + 'dev.csv')
    for i in range(len(df)):
        val_features.append(df.smiles[i])
        val_labels.append(df.activity[i])

    df = pd.read_csv(filepath + 'test.csv')
    for i in range(len(df)):
        test_features.append(df.smiles[i])
        test_labels.append(df.activity[i])

    if sample == True:
        over_sample(train_features, train_labels)

    return train_features, train_labels, val_features, val_labels,  test_features, test_labels

def load_data(filepath):
    train_features, train_labels, val_features, val_labels, test_features, test_labels = load_dataset(filepath, sample=True)
    train_labels = Variable(torch.IntTensor(train_labels), requires_grad=False)
    val_labels = Variable(torch.IntTensor(val_labels), requires_grad=False)
    test_labels = Variable(torch.IntTensor(test_labels), requires_grad=False)

    return train_features, train_labels, val_features, val_labels, test_features, test_labels

def construct_multigraph(smile):
    g = OrderedDict({})
    h = OrderedDict({})

    molecule = Chem.MolFromSmiles(smile)
    AllChem.ComputeGasteigerCharges(molecule)
    asa = rdMolDescriptors._CalcLabuteASAContribs(molecule)
    tpsa = rdMolDescriptors._CalcTPSAContribs(molecule)
    crippen = rdMolDescriptors._CalcCrippenContribs(molecule)
    for i in range(0, molecule.GetNumAtoms()):
        atom_i = molecule.GetAtomWithIdx(i)
        feature_i = dc.feat.graph_features.atom_features(atom_i)
        h[i] = Variable(torch.FloatTensor(feature_i)).view(1, 75)
        for j in range(0, molecule.GetNumAtoms()):
            e_ij = molecule.GetBondBetweenAtoms(i, j)
            if e_ij != None:
                e_ij = list(map(lambda x: 1 if x == True else 0, dc.feat.graph_features.bond_features(e_ij)))  # ADDED edge feat
                e_ij = Variable(torch.FloatTensor(e_ij).view(1, 6))
                if i not in g:
                    g[i] = []
                g[i].append((e_ij, j))
    return g, h


def construct_graph(smile):
    g = DGLGraph()
    molecule = Chem.MolFromSmiles(smile)
    g.add_nodes(molecule.GetNumAtoms())
    n_feature = pt.ones(molecule.GetNumAtoms(), 75)
    e_feature = pt.ones(molecule.GetNumBonds()*2, 6)
    edge_n = 0
    for i in range(0, molecule.GetNumAtoms()):
        atom_i = molecule.GetAtomWithIdx(i)
        n_feature[i] = Variable(pt.FloatTensor(dc.feat.graph_features.atom_features(atom_i))).view(1, 75)
        for j in range(0, molecule.GetNumAtoms()):
            e_ij = molecule.GetBondBetweenAtoms(i, j)
            if e_ij != None:
                e_ij = list(
                    map(lambda x: 1 if x == True else 0, dc.feat.graph_features.bond_features(e_ij)))  # ADDED edge feat
                e_ij = Variable(pt.FloatTensor(e_ij).view(1, 6))
                e_feature[edge_n] = e_ij
                edge_n += 1
                g.add_edge(i, j)
    return g, n_feature, e_feature

def construct_feat(smile, n, m):
    n_feature = pt.ones(n,  75)
    e_feature = pt.ones(m, 6)
    edge_n = 0
    molecule = Chem.MolFromSmiles(smile)
    for i in range(0, molecule.GetNumAtoms()):
        atom_i = molecule.GetAtomWithIdx(i)
        n_feature[i] = Variable(pt.FloatTensor(dc.feat.graph_features.atom_features(atom_i))).view(1, 75)
        for j in range(0, molecule.GetNumAtoms()):
            if j > i:
                e_ij = molecule.GetBondBetweenAtoms(i, j)
                if e_ij != None:
                    e_ij = list(map(lambda x: 1 if x == True else 0, dc.feat.graph_features.bond_features(e_ij)))  # ADDED edge feat
                    e_ij = Variable(pt.FloatTensor(e_ij).view(1, 6))
                    e_feature[edge_n] = e_ij
                    edge_n += 1
    return n_feature, e_feature