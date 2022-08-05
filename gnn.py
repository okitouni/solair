import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path as osp
import torch
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import NeighborLoader
import torch_geometric.transforms as T
from src.dataset.dataset_processing import CFD_dataset
from src.models.models import GCN
import torch.nn as nn
from tqdm import trange
import wandb

wandb.init(project="gnn", entity="imdea_dolo")

######### create dataset
dataset = CFD_dataset('./data/')
data = dataset[0]

# %%
train_loader = NeighborLoader(
    data,
    num_neighbors=[4]*6,
    batch_size=128,
    input_nodes=data.train_mask,
)
test_loader = NeighborLoader(
    data,
    num_neighbors=[4]*6,
    batch_size=128,
    shuffle=False,
    input_nodes=~data.train_mask,
)

######### create model 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = GCN(256, 4).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)#, weight_decay=5e-4)
loss_mse = nn.MSELoss()

######### train model
model.train()
losses_train = []
losses_val = []

EPOCHS = 100
pbar = trange(len(train_loader) * EPOCHS, desc='Training')
for epoch in range(EPOCHS):
    for i, data in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data)
        loss = loss_mse(out, data.y)
        wandb.log({"loss": loss_val})
        losses_train.append(loss.item())
        loss.backward()
        optimizer.step()
        # validation step
        if i % 100 == 0:
            model.eval()
            with torch.no_grad():
                for data in test_loader:
                    data = data.to(device)
                    out = model(data)
                    loss_val = loss_mse(out, data.y)
                    wandb.log({"loss_val": loss_val})
                    losses_val.append(loss_val.item())
            model.train()
    
        pbar.update(1)
        pbar.set_description(f'Epoch {epoch + 1:03d}| Loss train: {loss.item():.4f}, Loss val: {loss_val.item():.4f}')

np.savez('losses', train=losses_train, val=losses_val)
torch.save(model.state_dict(), 'model.pt')
