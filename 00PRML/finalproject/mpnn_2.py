from dgl.model_zoo.chem import MPNNModel
import numpy as np
import torch as pt
from torch.autograd import Variable
import torch.optim as optim
from torch.utils.tensorboard import SummaryWriter
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, roc_curve
import  matplotlib.pyplot as plt
import random
import preprocessing

random.seed(2)
pt.manual_seed(2)
np.random.seed(2)

writer = SummaryWriter('./dgl_model_1')


def train_mpnn(n_iter, n_batch, train_set, train_label, val_set, val_label, test_set, test_label, model, optimizer):
    for i in range(n_iter):
        y_hats_val = []
        y_val = []
        train_loss = Variable(pt.zeros(1, 1))
        for j in range(n_batch):
            sample_index = random.randint(0, len(train_set) - 1)
            smile = train_set[sample_index]
            g, n_feature, e_feature = preprocessing.construct_graph(smile)
            out = model(g, n_feature, e_feature)
            y_hats_val.append(out)
            y_val.append(train_label[sample_index])
            error = (out - train_label[sample_index]) * (out - train_label[sample_index]) / Variable(pt.FloatTensor([n_batch])).view(1, 1)
            train_loss += error

        train_loss.backward()
        optimizer.step()
        train_loss_ = train_loss.data.numpy()[0]
        writer.add_scalar('[%s] train_loss' % 'mpnn', train_loss_, i)

        if i % 10 == 0:
            y_hats_val = np.array(list(map(lambda x: x.data.numpy(), y_hats_val)))
            y_val = np.array(list(map(lambda x: x.data.numpy(), y_val)))
            y_hats_val = y_hats_val.reshape(-1, 1)
            y_val = y_val.reshape(-1, 1)
            roc_auc = roc_auc_score(y_val, y_hats_val)
            fpr, tpr, threshold = roc_curve(y_val, y_hats_val)
            plt.plot(fpr, tpr, color='red', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.legend(loc="lower right")
            plt.savefig('./fig/%d_mpnn_0_dgl_iter.png' % i)

        if i % int(len(train_set) / n_batch) == 0:
            val_loss = Variable(pt.zeros(1, 1), requires_grad=False)
            y_hats_val = []
            for j in range(len(val_set)):
                smile = train_set[j]
                g, n_feature, e_feature  = preprocessing.construct_graph(smile)
                out = model(g, n_feature, e_feature)
                y = val_label[j]
                y_hats_val.append(out)
                error = (out - y) * (out - y) / Variable(pt.FloatTensor([len(val_set)])).view(1, 1)
                val_loss = val_loss + error
            y_hats_val = np.array(list(map(lambda x: x.data.numpy(), y_hats_val)))
            y_val = np.array(list(map(lambda x: x.data.numpy(), val_label)))
            y_hats_val = y_hats_val.reshape(-1, 1)
            y_val = y_val.reshape(-1, 1)


            p, r, thr = precision_recall_curve(y_val, y_hats_val)
            prc_auc = auc(r, p)
            train_loss_ = train_loss.data.numpy()[0]
            val_loss_ = val_loss.data.numpy()[0]
            print('epoch [%i/%i] train_loss [%f] val_loss [%f] roc_auc [%.4f], prc_auc [%.4f]' \
            % (i, n_iter, train_loss_, val_loss_, roc_auc, prc_auc))
            writer.add_scalar('[%s] roc_auc' % 'mpnn', roc_auc, i)
            writer.add_scalar('[%s] prc_auc' % 'mpnn', prc_auc, i)
            # n_iter += 1
        else:
            print('epoch [%i/%i] train_loss [%f] ]' \
                  % (i, n_iter, train_loss_))

    y_hats_val = []
    for j in range(len(test_set)):
        smile = test_set[j]
        g, n_feature, e_feature = preprocessing.construct_graph(smile)
        out = model(g, n_feature, e_feature)
        y_hats_val.append(out)

    y_hats_val = np.array(list(map(lambda x: x.data.numpy(), y_hats_val)))
    y_val = np.array(list(map(lambda x: x.data.numpy(), test_label)))
    y_hats_val = y_hats_val.reshape(-1, 1)
    y_val = y_val.reshape(-1, 1)

    roc_auc = roc_auc_score(y_val, y_hats_val)
    p, r, thr = precision_recall_curve(y_val, y_hats_val)
    prc_auc = auc(r, p)
    print('roc_auc [%.4f], prc_auc [%.4f]' % (roc_auc, prc_auc))



if __name__ == '__main__':
    for i in range(1, 10):
        train_fold = i
        train_features, train_labels, val_features, val_labels, test_features, test_labels = donkey.load_dataset('./train_cv/fold_%d/'%train_fold)
        train_labels = Variable(pt.IntTensor(train_labels), requires_grad=False)
        val_labels = Variable(pt.IntTensor(val_labels), requires_grad=False)
        test_labels = Variable(pt.IntTensor(test_labels), requires_grad=False)
        model = MPNNModel(
                node_input_dim=75,
                edge_input_dim=6,
                output_dim=1,
                node_hidden_dim = 50,
                edge_hidden_dim = 6,
                num_step_message_passing = 3,
                num_step_set2set = 6,
                num_layer_set2set = 3)
        optimizer = optim.Adam(model.parameters(), lr=1e-4)
        train_mpnn(100, 50, train_features, train_labels, val_features, val_labels, test_features, test_labels, model, optimizer)