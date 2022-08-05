import torch.nn.functional as F
from torch_geometric.nn import GCNConv
import torch

class GCN(torch.nn.Module):
    def __init__(self, hidden_channels=128, hidden_layers=3):
        super().__init__()
        self.layers = [GCNConv(6, hidden_channels, cached=False)]
        for _ in range(hidden_layers - 1):
            self.layers.append(GCNConv(hidden_channels, hidden_channels, cached=False))
        self.layers = torch.nn.ModuleList(self.layers)

        self.regressor = torch.nn.Sequential(
            torch.nn.Linear(hidden_channels, hidden_channels),
            torch.nn.LeakyReLU(),
            torch.nn.Linear(hidden_channels, hidden_channels),
            torch.nn.LeakyReLU(),
            torch.nn.Linear(hidden_channels, 1),
        )

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        for layer in self.layers:
            x = layer(x, edge_index)
            x = F.relu(x)
        x = self.regressor(x)
        return x