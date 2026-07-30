import gc
import random
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
DEFAULT_DTYPE = torch.float32


def get_device(requested_device=None):
    if requested_device is None or requested_device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(requested_device)


def get_device_info(requested_device=None):
    device = get_device(requested_device)
    return {
        "device": str(device),
        "cuda_available": torch.cuda.is_available(),
        "gpu_name": torch.cuda.get_device_name(0) if device.type == "cuda" else None,
        "dtype": str(DEFAULT_DTYPE).replace("torch.", "")
    }


def set_all_seeds(seed):
    seed = int(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    return seed


def clear_all():
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True


class HybridPCFModel(nn.Module):
    def __init__(self, pcf_channels, seq_len, out_dim, n_aux):
        super(HybridPCFModel, self).__init__()

        self.pcf_branch = nn.Sequential(
            nn.Conv1d(pcf_channels, 64, kernel_size=7, padding=0),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=5),
            nn.Conv1d(64, 64, kernel_size=7, padding=0),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=5),
            nn.Conv1d(64, 64, kernel_size=7, padding=0),
            nn.ReLU()
        )

        self.n_aux = int(n_aux)

        with torch.no_grad():
            dummy_pcf = torch.zeros(1, pcf_channels, seq_len, dtype=DEFAULT_DTYPE)
            pcf_dim = int(self.pcf_branch(dummy_pcf).view(1, -1).shape[1])
            flat_feats = pcf_dim + self.n_aux

        self.head = nn.Sequential(
            nn.Linear(flat_feats, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, out_dim)
        )

    def forward(self, x_pcf, x_aux):
        z_pcf = self.pcf_branch(x_pcf)
        z_pcf = z_pcf.view(z_pcf.size(0), -1)
        if x_aux.dim() == 1:
            x_aux = x_aux.unsqueeze(1)
        z = torch.cat([z_pcf, x_aux], dim=1)
        return self.head(z)


def from_batch_to_tensor_pcf(x, device):
    x = np.asarray(x, dtype=np.float32)
    x = np.ascontiguousarray(x)
    return torch.from_numpy(x).permute(0, 2, 1).contiguous().to(device=device, dtype=DEFAULT_DTYPE)


def from_batch_to_tensor_aux(x, device):
    x = np.asarray(x, dtype=np.float32)
    x = np.ascontiguousarray(x)
    return torch.from_numpy(x).to(device=device, dtype=DEFAULT_DTYPE)


def train_and_save_hybrid_NN_model(L_train, Y_train, AUX_train, model_path, batch_size=64, epochs=30, lr=1e-3, requested_device="auto", seed=None, verbose=True):
    batch_size = int(batch_size)
    epochs = int(epochs)
    lr = float(lr)
    device = get_device(requested_device)
    if seed is not None:
        set_all_seeds(seed)
    L_train = np.asarray(L_train)
    Y_train = np.asarray(Y_train)
    AUX_train = np.asarray(AUX_train)
    N = int(L_train.shape[0])

    pcf_channels = int(L_train.shape[2])
    seq_len = int(L_train.shape[1])
    n_out = int(Y_train.shape[1])
    n_aux = int(AUX_train.shape[1])
    model = HybridPCFModel(pcf_channels=pcf_channels, seq_len=seq_len, out_dim=n_out, n_aux=n_aux).to(device=device, dtype=DEFAULT_DTYPE)
    opt = optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.MSELoss()

    def run_epoch(train=True):
        model.train(train)
        total = 0

        for start in range(0, N, batch_size):
            end = min(start + batch_size, N)
            Lb = from_batch_to_tensor_pcf(L_train[start:end], device)
            Ab = from_batch_to_tensor_aux(AUX_train[start:end], device)
            Yb = from_batch_to_tensor_aux(Y_train[start:end], device)
            pred = model(Lb, Ab)
            pred_sigma = pred[:, [0, 2, 3]]
            true_sigma = Yb[:, [0, 2, 3]]
            pred_scale = pred[:, [1, 4, 5]]
            true_scale = Yb[:, [1, 4, 5]]
            loss_sigma = loss_fn(pred_sigma, true_sigma)
            loss_scale = loss_fn(pred_scale, true_scale)
            loss = loss_sigma + loss_scale

            if train:
                opt.zero_grad()
                loss.backward()
                opt.step()

            total += float(loss.item()) * int(end - start)

            del Lb, Ab, Yb, pred, pred_sigma, true_sigma, pred_scale, true_scale, loss_sigma, loss_scale, loss

        clear_all()
        return total / N

    train_list = []
    if verbose:
        info = get_device_info(requested_device)
        print(f"device={info['device']} dtype={info['dtype']} gpu_name={info['gpu_name']}")

    for ep in range(epochs):
        train_loss = run_epoch(train=True)
        train_list.append(train_loss)
        if verbose:
            print(f"epoch {ep}  train_mse={train_loss:.4f}")

    save_all = {
        "state_dict": model.state_dict(),
        "pcf_channels": pcf_channels,
        "seq_len": seq_len,
        "n_out": n_out,
        "n_aux": n_aux,
        "dtype": "float32",
        "requested_device": requested_device,
    }
    torch.save(save_all, model_path)

    return {
        "device": str(device)
    }


def load_hybrid_NN_model_predict(L_test, AUX_test, model_path, requested_device="auto", predict_batch_size=2000):
    predict_batch_size = int(predict_batch_size)
    device = get_device(requested_device)
    load_model = torch.load(model_path, map_location=device)
    model = HybridPCFModel(pcf_channels=load_model["pcf_channels"], seq_len=load_model["seq_len"], out_dim=load_model["n_out"], n_aux=load_model["n_aux"]).to(device=device, dtype=DEFAULT_DTYPE)
    model.load_state_dict(load_model["state_dict"])
    model.eval()
    L_test = np.asarray(L_test)
    AUX_test = np.asarray(AUX_test)
    N = int(L_test.shape[0])
    n_out = int(load_model["n_out"])
    preds_out = np.empty((N, n_out), dtype=np.float32)

    with torch.no_grad():
        for start in range(0, N, predict_batch_size):
            end = min(start + predict_batch_size, N)
            Lb = from_batch_to_tensor_pcf(L_test[start:end], device)
            Ab = from_batch_to_tensor_aux(AUX_test[start:end], device)
            pred = model(Lb, Ab).detach().cpu().numpy()
            preds_out[start:end, :] = pred

            del Lb, Ab, pred

    del model, load_model
    clear_all()
    return preds_out
