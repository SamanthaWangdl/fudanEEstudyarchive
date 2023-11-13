import matplotlib.pyplot as plt
import torch
import torch.nn as nn
from torch.autograd import Variable
import torch.optim as optim
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, roc_curve
import numpy as np
from torch.utils.tensorboard import SummaryWriter
import random
import preprocessing

random.seed(2)
torch.manual_seed(2)
np.random.seed(2)


T = 5
BATCH_SIZE = 50
MAXITER = 150
LIMIT = 0
LR = 1e-4

class MPNN(nn.Module):
    def __init__(self, n_input_dim = 75,
                 n_layer=3,
                 e_input_dim=6,
                 e_hidden_dim = 6,
                 n_hidden_dim_1 = 50,
                 n_hidden_dim_2 = 80):
        super(MPNN, self).__init__()
        self.n_layer = n_layer
        self.e_hidden_dim = e_hidden_dim
        self.n_hidden_dim_1 = n_hidden_dim_1
        self.n_hidden_dim_2 = n_hidden_dim_2
        self.R = nn.Linear(2 * n_input_dim, n_hidden_dim_2)
        self.U = {i: nn.Linear(n_input_dim + n_hidden_dim_1 + e_hidden_dim, n_input_dim) for i in range(self.n_layer)}
        self.V = {i: nn.Linear(n_input_dim, n_hidden_dim_1) for i in range(self.n_layer)}
        self.E = nn.Linear(e_input_dim, e_hidden_dim)
        self.lin0 = nn.Linear(n_hidden_dim_2, 1)

    def forward(self,g,  n_feature):
        h = n_feature
        for step in range(self.n_layer):
            for v in g.keys():
                neighbors = g[v]
                for neighbor in neighbors:
                    e_vw = neighbor[0]  # feature of edge
                    w = neighbor[1]
                    m_w = self.V[k](h[w])  # 输出为V dim
                    m_e_vw = self.E(e_vw)  # 输出为E dim
                    reshaped = torch.cat((h[v], m_w, m_e_vw), 1)
                    h[v] = F.selu(U[k](reshaped))
        catted_reads = map(lambda x: torch.cat([h[x[0]], n_feature[x[1]]], 1), zip(n_feature.keys(), h.keys()))
        activated_reads = map(lambda x: F.selu(self.R(x)), catted_reads)
        readout = Variable(torch.zeros(1, self.n_hidden_dim_2))
        for read in activated_reads:
            readout = readout + read
        readout = torch.tanh(readout)
        return self.lin0(readout)


def message_pass(g, h, k):
    for v in g.keys():
        neighbors = g[v]
        for neighbor in neighbors:
            e_vw = neighbor[0]  # feature variable
            w = neighbor[1]
            m_w = V[k](h[w]) # 输出为V dim
            m_e_vw = E(e_vw) # 输出为E dim
            reshaped = torch.cat((h[v], m_w, m_e_vw), 1)
            h[v] = F.selu(U[k](reshaped))


def adjust_learning_rate(optimizer, epoch):
    """Sets the learning rate to the initial LR decayed by .8 every 5 epochs"""
    lr = LR * (0.9 ** (epoch // 10))
    print('new lr [%.5f]' % lr)
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr



def train_mpnn(model, optimizer,  train_smiles, train_labels, val_smiles, val_labels, test_fold):
    for i in range(0, MAXITER):
        optimizer.zero_grad()
        train_loss = Variable(torch.zeros(1, 1))
        y_hats_train = []
        y_val = []
        for j in range(0, BATCH_SIZE):
            sample_index = random.randint(0, len(train_smiles) - 2)
            smile = train_smiles[sample_index]
            g, h = preprocessing.construct_multigraph(smile)
            y_hat = model(g, h)
            y = train_labels[sample_index]
            y_hats_train.append(y_hat)
            y_val.append(y)
            error = (y_hat - y) * (y_hat - y) / Variable(torch.FloatTensor([BATCH_SIZE])).view(1, 1)
            train_loss = train_loss + error

        train_loss.backward()
        optimizer.step()
        # writer.add_scalar('[%s/%d]train_loss' % ('mpnn', test_fold), np.mean(train_loss.data.numpy()), i)
        if i % 10 == 0:
            y_hats_train = np.array(list(map(lambda x: x.data.numpy(), y_hats_train)))
            y_val = np.array(list(map(lambda x: x.data.numpy(), y_val)))
            y_hats_train = y_hats_train.reshape(-1, 1)
            y_val = y_val.reshape(-1, 1)
            roc_auc = roc_auc_score(y_val, y_hats_train)
            p, r, thr = precision_recall_curve(y_val, y_hats_train)
            prc_auc = auc(r, p)
            train_loss_ = np.mean(train_loss.data.numpy())
            print('epoch [%i/%i] train_loss [%f] roc_auc [%.4f], prc_auc [%.4f]' \
                  % (i, MAXITER, train_loss_, roc_auc, prc_auc))
            fpr, tpr, threshold = roc_curve(y_val, y_hats_train)
            plt.clf()
            plt.plot(fpr, tpr, color='red', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.legend(loc="lower right")
            plt.savefig('./fig/%d_mpnn0_iter_%d.png' % (test_fold, i))

            plt.clf()
            plt.plot(r, p, color='red', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.legend(loc="lower right")
            plt.savefig('./fig/%d_prc_iter_%d.png' % (test_fold, i))

        else:
            train_loss_ = np.mean(train_loss.data.numpy())
            print('epoch [%i/%i] train_loss [%f]'%(i, MAXITER, train_loss_))


        val_loss = Variable(torch.zeros(1, 1), requires_grad=False)
        y_hats_val = []
        for j in range(0, len(val_smiles)):
            g, h = preprocessing.construct_multigraph(val_smiles[j])
            y_hat = model(g,h)
            y = val_labels[j]
            y_hats_val.append(y_hat)
            error = (y_hat - y) * (y_hat - y) / Variable(torch.FloatTensor([len(val_smiles)])).view(1, 1)
            val_loss = val_loss + error

    y_hats_val = np.array(list(map(lambda x: x.data.numpy(), y_hats_val)))
    y_val = np.array(list(map(lambda x: x.data.numpy(), val_labels)))
    y_hats_val = y_hats_val.reshape(-1, 1)
    y_val = y_val.reshape(-1, 1)
    roc_auc = roc_auc_score(y_val, y_hats_val)
    p, r, thr = precision_recall_curve(y_val, y_hats_val)
    prc_auc = auc(r, p)
    fpr, tpr, threshold = roc_curve(y_val, y_hats_val)
    plt.clf()
    plt.plot(fpr, tpr, color='red', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig('./fig/%d_mpnn0_val.png' % (test_fold))
    print('Val: loss [%.4f] roc_auc [%.4f], prc_auc [%.4f]' % (val_loss, roc_auc, prc_auc))



if __name__ == '__main__':

    for k in range(10):
        test_fold = k
        train_smiles, train_labels, val_smiles, val_labels, test_smiles, test_labels = preprocessing.load_data('./train_cv/fold_%d/'%k)

        linear = nn.Linear(80, 1)
        model = MPNN(n_input_dim= 75,
                     n_layer = 3,
                     e_input_dim=6,
                     e_hidden_dim=6,
                     n_hidden_dim_1=50,
                     n_hidden_dim_2=80)
        #
        params = [{'params': model.R.parameters()},
                  {'params': model.U[0].parameters()},
                  {'params': model.U[1].parameters()},
                  {'params': model.U[2].parameters()},
                  {'params': model.E.parameters()},
                  {'params': model.V[0].parameters()},
                  {'params': model.V[1].parameters()},
                  {'params': model.V[2].parameters()},
                  {'params': model.lin0.parameters()}]

        optimizer = optim.Adam(params, lr=LR, weight_decay=1e-4)
        train_mpnn(model, optimizer, train_smiles, train_labels, val_smiles, val_labels, test_fold)
