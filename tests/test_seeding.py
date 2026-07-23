import os
import random

import numpy as np
import torch

from signflow.config import seed_all


def test_seed_all_is_deterministic():
    seed_all(7)
    first = (random.random(), float(np.random.rand()), torch.rand(3))
    seed_all(7)
    second = (random.random(), float(np.random.rand()), torch.rand(3))
    assert first[0] == second[0]
    assert first[1] == second[1]
    assert torch.equal(first[2], second[2])
    assert os.environ["PYTHONHASHSEED"] == "7"
    assert torch.are_deterministic_algorithms_enabled() is True


def test_seed_all_changes_torch_stream_for_different_seeds():
    seed_all(7)
    first = torch.rand(5)
    seed_all(8)
    second = torch.rand(5)
    assert torch.equal(first, second) is False
