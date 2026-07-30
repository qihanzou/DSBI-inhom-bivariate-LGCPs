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


class ImageBranch(nn.Module):
    def __init__(self, in_channels):
        super(ImageBranch, self).__init__()

        self.branch = nn.Sequential(
            nn.Conv2d(in_channels, 32, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),
            nn.Conv2d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),
            nn.Conv2d(64, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool2d((4, 4))
        )

    def forward(self, x):
        z = self.branch(x)
        return z.view(z.size(0), -1)


class HybridImageModel(nn.Module):
    def __init__(self, img_channels, img_h, img_w, out_dim):
        super(HybridImageModel, self).__init__()

        self.img_branch = ImageBranch(img_channels)

        with torch.no_grad():
            dummy_img = torch.zeros(1, img_channels, img_h, img_w, dtype=DEFAULT_DTYPE)
            img_dim = int(self.img_branch(dummy_img).view(1, -1).shape[1])
            flat_feats = img_dim

        self.head = nn.Sequential(
            nn.Linear(flat_feats, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, out_dim)
        )

    def forward(self, x_img):
        z_img = self.img_branch(x_img)
        return self.head(z_img)


def get_img_norm(x_img):
    x_img = np.asarray(x_img)
    mean = x_img.mean(axis=(0, 1, 2), keepdims=True).astype(np.float32)
    std = x_img.std(axis=(0, 1, 2), keepdims=True).astype(np.float32)
    std[std == 0] = 1
    return mean, std


def from_batch_to_tensor_img(x, device, img_mean=None, img_std=None, standardize_img=True):
    x = np.asarray(x, dtype=np.float32)
    if standardize_img:
        x = (x - img_mean) / img_std
    x = np.ascontiguousarray(x)
    return torch.from_numpy(x).permute(0, 3, 1, 2).contiguous().to(device=device, dtype=DEFAULT_DTYPE)


def from_batch_to_tensor_aux(x, device):
    x = np.asarray(x, dtype=np.float32)
    x = np.ascontiguousarray(x)
    return torch.from_numpy(x).to(device=device, dtype=DEFAULT_DTYPE)


def train_and_save_hybrid_NN_model(IMG_train, Y_train, model_path, batch_size=64, epochs=30, lr=1e-3, requested_device="auto", seed = None, verbose=True, standardize_img=True):
    batch_size = int(batch_size)
    epochs = int(epochs)
    lr = float(lr)
    standardize_img = bool(standardize_img)
    device = get_device(requested_device)
    if seed is not None:
       set_all_seeds(seed)
    IMG_train = np.asarray(IMG_train)
    Y_train = np.asarray(Y_train)
    N = int(IMG_train.shape[0])
  
    if standardize_img:
        img_mean, img_std = get_img_norm(IMG_train)
    else:
        img_mean = np.zeros((1, 1, 1, int(IMG_train.shape[3])), dtype=np.float32)
        img_std = np.ones((1, 1, 1, int(IMG_train.shape[3])), dtype=np.float32)

    img_h = int(IMG_train.shape[1])
    img_w = int(IMG_train.shape[2])
    img_channels = int(IMG_train.shape[3])
    n_out = int(Y_train.shape[1])
    model = HybridImageModel(img_channels=img_channels, img_h=img_h, img_w=img_w, out_dim=n_out).to(device=device, dtype=DEFAULT_DTYPE)
    opt = optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.MSELoss()

    def run_epoch(train=True):
        model.train(train)
        total = 0

        for start in range(0, N, batch_size):
            end = min(start + batch_size, N)
            Ib = from_batch_to_tensor_img(IMG_train[start:end], device, img_mean=img_mean, img_std=img_std, standardize_img=standardize_img)
            Yb = from_batch_to_tensor_aux(Y_train[start:end], device)
            pred = model(Ib)
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

            del Ib, Yb, pred, pred_sigma, true_sigma, pred_scale, true_scale, loss_sigma, loss_scale, loss

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
        "img_channels": img_channels,
        "img_h": img_h,
        "img_w": img_w,
        "n_out": n_out,
        "img_mean": img_mean,
        "img_std": img_std,
        "standardize_img": standardize_img,
        "dtype": "float32",
        "requested_device": requested_device,
    }
    torch.save(save_all, model_path)

    return {
        "device": str(device)
    }


def load_hybrid_NN_model_predict(IMG_test, model_path, requested_device="auto", predict_batch_size=2000):
    predict_batch_size = int(predict_batch_size)
    device = get_device(requested_device)
    load_model = torch.load(model_path, map_location=device)
    model = HybridImageModel(img_channels=load_model["img_channels"], img_h=load_model["img_h"], img_w=load_model["img_w"], out_dim=load_model["n_out"]).to(device=device, dtype=DEFAULT_DTYPE)
    model.load_state_dict(load_model["state_dict"])
    model.eval()
    IMG_test = np.asarray(IMG_test)
    N = int(IMG_test.shape[0])
    n_out = int(load_model["n_out"])
    preds_out = np.empty((N, n_out), dtype=np.float32)
    standardize_img = bool(load_model.get("standardize_img", True))
    img_mean = load_model["img_mean"]
    img_std = load_model["img_std"]

    with torch.no_grad():
        for start in range(0, N, predict_batch_size):
            end = min(start + predict_batch_size, N)
            Ib = from_batch_to_tensor_img(IMG_test[start:end], device, img_mean=img_mean, img_std=img_std, standardize_img=standardize_img)
            pred = model(Ib).detach().cpu().numpy()
            preds_out[start:end, :] = pred

            del Ib, pred

    del model, load_model
    clear_all()
    return preds_out
