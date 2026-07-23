import random

import numpy as np

from signflow.config import seed_all


def test_seed_all_is_deterministic():
    seed_all(7)
    first = (random.random(), float(np.random.rand()))
    seed_all(7)
    second = (random.random(), float(np.random.rand()))
    assert first == second
