from rdkit import Chem
import numpy as np
from torch.autograd import Variable
import deepchem as dc
from torch.utils.data import Dataset, DataLoader
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
import preprocessing
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, roc_curve


class GConvRegressor(torch.nn.Module):
    def __init__(self, n_feature=128, n_layer=10):
        super(GConvRegressor, self).__init__()
        self.W = nn.ModuleList([nn.Linear(n_feature, n_feature) for _ in range(n_layer)])
        self.embedding = nn.Linear(75, n_feature)
        self.fc = nn.Linear(n_feature, 1)

    def forward(self, x, A):
        x = self.embedding(x)
        for l in self.W:
            x = l(x)
            x = torch.einsum('ijk,ikl->ijl', (A.clone(), x))
            x = F.relu(x)
        x = x.mean(1)

        retval = self.fc(x)

        return retval


class MolDataset(Dataset):
    def __init__(self, smiles, max_natoms, labels):
        self.smiles = smiles
        self.max_natoms = max_natoms
        self.label = labels

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        s = self.smiles[idx]

        m = Chem.MolFromSmiles(s)
        natoms = m.GetNumAtoms()

        # adjacency matrix
        A = GetAdjacencyMatrix(m) + np.eye(natoms)
        A_padding = np.zeros((self.max_natoms, self.max_natoms))
        A_padding[:natoms, :natoms] = A

        # atom feature
        X = [dc.feat.graph_features.atom_features(m.GetAtomWithIdx(i)) for i in range(natoms)]
        for i in range(natoms, max_natoms):
            X.append(np.zeros(75))
        X = np.array(X)

        sample = dict()
        sample['X'] = torch.from_numpy(X)
        sample['A'] = torch.from_numpy(A_padding)
        sample['Y'] = self.label[idx]
        return sample



def load_dataset(filepath):
    train_features, train_labels, val_features, val_labels, test_features, test_labels = preprocessing.load_dataset(filepath)
    train_labels = Variable(torch.IntTensor(train_labels), requires_grad=False)
    val_labels = Variable(torch.IntTensor(val_labels), requires_grad=False)
    test_labels = Variable(torch.IntTensor(test_labels), requires_grad=False)


    return train_features, train_labels, val_features, val_labels, test_features, test_labels



if __name__ == '__main__':
    max_natoms = 200
    for i in range(10):
        train_fold = i
        train_smiles, train_labels, val_smiles, val_labels, test_smiles, test_labels = load_dataset('./train_cv/fold_%d/'%i)
        lr = 1e-4
        model = GConvRegressor(50, 5)
        train_dataset = MolDataset(train_smiles, max_natoms,train_labels)
        val_dataset = MolDataset(test_smiles, max_natoms, test_labels)
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        train_dataloader = DataLoader(train_dataset, batch_size=48, num_workers=0)
        val_dataloader = DataLoader(val_dataset, batch_size=48, num_workers=0)
        loss_fn = nn.MSELoss()
        loss_list = []
        for epoch in range(100):
            epoch_loss = []
            pred_set = []
            y_set = []
            for i_batch, batch in enumerate(train_dataloader):
                x, y, A = batch['X'].float(), batch['Y'].float(), batch['A'].float()
                pred = model(x, A).squeeze(-1)
                loss = loss_fn(pred, y)
                for j in range(pred.shape[0]):
                    pred_set.append(pred[j])
                    y_set.append(y[j])
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)

                optimizer.step()
                loss_list.append(loss.data.cpu().numpy())
                epoch_loss.append(loss.data.cpu().numpy())
                # writer.add_scalar('[%s/%d] train_loss' % ('gcn',train_fold), np.mean(np.array(epoch_loss)),epoch)

            if epoch % 10 == 0:
                y_set = np.array(list(map(lambda x: x.data.numpy(), y_set)))
                pred_set = np.array(list(map(lambda x: x.data.numpy(), pred_set)))
                roc_auc = roc_auc_score(y_set, pred_set)
                p, r, thr = precision_recall_curve(y_set, pred_set)
                prc_auc = auc(r, p)
                # writer.add_scalar('[%s/%d] roc_auc' %  ('gcn',train_fold), roc_auc, epoch)
                # writer.add_scalar('[%s/%d] prc_auc' %  ('gcn',train_fold), prc_auc, epoch)
                print('epoch [%i/%i] train_loss [%.6f] roc_auc [%.4f], prc_auc [%.4f]'
                      % (epoch, 100, np.mean(np.array(epoch_loss)), roc_auc, prc_auc))
                fpr, tpr, threshold = roc_curve(y_set, pred_set)
                plt.clf()
                plt.plot(fpr, tpr, color='red', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.legend(loc="lower right")
                plt.savefig('./fig/%d_gcn_test_%d.png'%(train_fold, epoch))

                plt.clf()
                plt.plot(r, p, color='red', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('Recall')
                plt.ylabel('Precision')
                plt.legend(loc="lower right")
                plt.savefig('./fig/%d_gcn_prc_%d.png' % (train_fold, epoch))

            else:
                print('epoch [%i/%i] train_loss [%.6f]' % (epoch, 100, np.mean(np.array(epoch_loss))))

        loss_list = []
        pred_set = []
        y_set = []
        for i, batch in enumerate(val_dataloader):
            x, y, A = batch['X'].float(), batch['Y'].float(), batch['A'].float()
            pred = model(x, A).squeeze(-1)
            loss = loss_fn(pred, y)
            for j in range(pred.shape[0]):
                pred_set.append(pred[j])
                y_set.append(y[j])

            loss_list.append(loss.data.cpu().numpy())

        y_set = np.array(list(map(lambda x: x.data.numpy(), y_set)))
        pred_set = np.array(list(map(lambda x: x.data.numpy(), pred_set)))
        roc_auc = roc_auc_score(y_set, pred_set)
        p, r, thr = precision_recall_curve(y_set, pred_set)
        prc_auc = auc(r, p)
        plt.clf()
        fpr, tpr, threshold = roc_curve(y_set, pred_set)
        plt.plot(fpr, tpr, color='red', lw= 2, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic example')
        plt.legend(loc="lower right")
        plt.savefig('./fig/%d_gcn_test.png'%(train_fold))

        print('val_loss[%.6f] roc_auc [%.4f], prc_auc [%.4f]' % (np.mean(np.array(loss_list)), roc_auc, prc_auc))

